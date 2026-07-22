function out = variables_organizer(x, k_vec, f_nodes, S_halo, target_position, epoch_comet_flyby, c_const)
% variables_organizer  Extract and organize all physical quantities of the
% optimal multiple-shooting solution (output of optimization_ephe with MS).
%
% Inputs:
%   x                 - design variable vector (result.bestfeasible.x)
%   k_vec             - [k1 k2 k3 k4] nodes per segment
%   f_nodes           - 1x4 cell, fixed node time fractions
%                       (epochs reconstructed: T = Ts + f*(Te-Ts))
%   S_halo            - halo state (synodic non-dim, 6x1)
%   target_position   - comet position [km] in ECLIPJ2000
%   epoch_comet_flyby - comet arrival epoch [ET s]
%   c_const           - constants struct (mu, Lstar, Tstar, G, mEarth, mSun)
%
% Outputs:
%   out - struct with fields:
%     .departure  : epoch [ET s], state (.syn 1x6, .J2000 1x6 [km,km/s])
%     .injection  : dv (.J2000_kms 1x3, .syn_adim 1x3, .norm_ms double)
%     .dsm1       : epoch [ET s], state_pre (.syn, .J2000),
%                   dv (.J2000_kms, .syn_adim, .norm_ms)
%     .flyby      : epoch [ET s], moon_state [km,km/s],
%                   state_pre (.syn, .J2000), state_post (.syn, .J2000),
%                   vinf_in [km/s], vinf_out [km/s]
%     .dsm2       : epoch [ET s], state_pre (.syn, .J2000),
%                   dv (.J2000_kms, .syn_adim, .norm_ms)
%     .comet      : epoch [ET s], state (.syn 1x6, .J2000 1x6)

mu    = c_const.mu;
Tstar = c_const.Tstar;

%% ===== UNPACK BASE DESIGN VARIABLES =====================================
dv_inj_J2000  = x(1:3).'   ./ 1e3;    % [km/s] column
dv_dsm1_J2000 = x(4:6).'   ./ 10;     % [km/s] column
epoch_flyby   = epoch_comet_flyby - x(7)  * 86400;   % [ET s]  (x in days before comet)
vel_comet_arr = x(8:10).';             % [km/s] column
dv_dsm2_J2000 = x(11:13).' ./ 10;     % [km/s] column
t_halo_days   = x(14);                % [days]
epoch_dep     = epoch_comet_flyby - x(15) * 86400;   % [ET s]  (x in days before comet)
epoch_dsm1    = epoch_comet_flyby - x(16) * 86400;   % [ET s]  DSM1 epoch (explicit)
epoch_dsm2    = epoch_comet_flyby - x(17) * 86400;   % [ET s]  DSM2 epoch (explicit)

%% ===== UNPACK MULTIPLE SHOOTING NODES ===================================
k1 = k_vec(1);  k2 = k_vec(2);
k3 = k_vec(3);  k4 = k_vec(4);

% Node states (node epochs are NOT in x: reconstructed from f_nodes)
n_base = 17;
off1 = n_base;          off2 = off1 + 6*k1;
off3 = off2 + 6*k2;

S1 = reshape(x(off1+1 : off1+6*k1), 6, k1);
S2 = reshape(x(off2+1 : off2+6*k2), 6, k2);
S3 = reshape(x(off3+1 : off3+6*k3), 6, k3);

% Node epochs reconstructed from fixed fractions: T = Ts + f*(Te - Ts).
% Only T2 (last forward node -> flyby) and T3 (backward segment) are needed here.
T2 = epoch_dsm1  + f_nodes{2}(:) .* (epoch_flyby - epoch_dsm1);
T3 = epoch_flyby + f_nodes{3}(:) .* (epoch_dsm2  - epoch_flyby);

opt = odeset('AbsTol',1e-13,'RelTol',1e-13);

%% ===== 1. DEPARTURE FROM THE HALO =======================================
tof_halo_ad = t_halo_days * 86400 / Tstar;   % non-dim CR3BP

if abs(tof_halo_ad) > 3600/Tstar
    [~, S_prop] = ode45(@(t,S) CR3BP(t,S,mu), [0 tof_halo_ad], S_halo(:), opt);
    S_dep_syn = S_prop(end,:);   % 1×6
else
    S_dep_syn = S_halo(:).';
end

S_dep_J2000 = synodic2sun_J2000(S_dep_syn, epoch_dep, mu);   % 1×6

out.departure.epoch = epoch_dep;
out.departure.state.syn   = S_dep_syn;
out.departure.state.J2000 = S_dep_J2000;

%% ===== 2. INJECTION MANEUVER ============================================
out.injection.J2000_kms  = dv_inj_J2000.';                       % 1×3 [km/s]
out.injection.syn_adim   = dv_J2000_to_syn(dv_inj_J2000, epoch_dep, mu).'; % 1×3
out.injection.norm_ms    = norm(dv_inj_J2000) * 1e3;             % [m/s]

%% ===== 3. PRE-DSM1 STATE ================================================
S_preDSM1_syn   = S1(:,k1).';                                          % 1×6
S_preDSM1_J2000 = synodic2sun_J2000(S_preDSM1_syn, epoch_dsm1, mu);   % 1×6

out.dsm1.epoch          = epoch_dsm1;
out.dsm1.state_pre.syn   = S_preDSM1_syn;
out.dsm1.state_pre.J2000 = S_preDSM1_J2000;
out.dsm1.dv.J2000_kms    = dv_dsm1_J2000.';
out.dsm1.dv.syn_adim     = dv_J2000_to_syn(dv_dsm1_J2000, epoch_dsm1, mu).';
out.dsm1.dv.norm_ms      = norm(dv_dsm1_J2000) * 1e3;

%% ===== 4. PRE-FLYBY STATE (forward, last arc -> flyby) ================
% Propagate from the last forward node to epoch_flyby.
% Start point: post-DSM1 (if k2=0) or the last T2 node (if k2>0).
if k2 > 0
    S_fwd_start_syn = S2(:,k2);
    T_fwd_start     = T2(k2);
else
    % Apply DSM1 to the last T1 node and start from there
    S_postDSM1_J2000 = synodic2sun_J2000(S_preDSM1_syn, epoch_dsm1, mu);
    S_postDSM1_J2000 = S_postDSM1_J2000(:);
    S_postDSM1_J2000(4:6) = S_postDSM1_J2000(4:6) + dv_dsm1_J2000;
    S_fwd_start_syn = [];   % flag: start in J2000
    T_fwd_start     = epoch_dsm1;
end

dt_fwd = epoch_flyby - T_fwd_start;

if ~isempty(S_fwd_start_syn)
    S0_fwd = synodic2sun_J2000(S_fwd_start_syn(:).', T_fwd_start, mu);
    S0_fwd = S0_fwd(:);
else
    S0_fwd = S_postDSM1_J2000;
end

[~, S_fwd] = ode45(@(t,S) NBODY_J2000(t,S,T_fwd_start,c_const), [0 dt_fwd], S0_fwd, opt);

S_preFlyby_J2000 = S_fwd(end,:);                                    % 1×6
S_preFlyby_syn   = sun_J2000_to_synodic(S_preFlyby_J2000, epoch_flyby, mu); % 1×6

out.flyby.epoch              = epoch_flyby;
out.flyby.state_pre.syn      = S_preFlyby_syn;
out.flyby.state_pre.J2000    = S_preFlyby_J2000;

%% ===== 5. MOON STATE AT FLYBY ===========================================
moon_state = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN').'; % 1×6
out.flyby.moon_state = moon_state;

%% ===== 6. POST-FLYBY STATE (backward, last T3 node -> flyby) ===========
S_bwd_start_syn = S3(:,k3);
T_bwd_start     = T3(k3);
dt_bwd          = epoch_flyby - T_bwd_start;   % negative (backward)

S0_bwd = synodic2sun_J2000(S_bwd_start_syn(:).', T_bwd_start, mu);
S0_bwd = S0_bwd(:);

[~, S_bwd] = ode45(@(t,S) NBODY_J2000(t,S,T_bwd_start,c_const), [0 dt_bwd], S0_bwd, opt);

S_postFlyby_J2000 = S_bwd(end,:);                                     % 1×6
S_postFlyby_syn   = sun_J2000_to_synodic(S_postFlyby_J2000, epoch_flyby, mu);

out.flyby.state_post.syn     = S_postFlyby_syn;
out.flyby.state_post.J2000   = S_postFlyby_J2000;

%% ===== 7. VINF IN / OUT =================================================
out.flyby.vinf_in  = S_preFlyby_J2000(4:6)  - moon_state(4:6);   % [km/s]
out.flyby.vinf_out = S_postFlyby_J2000(4:6) - moon_state(4:6);   % [km/s]
out.flyby.vinf_in_norm_kms  = norm(out.flyby.vinf_in);
out.flyby.vinf_out_norm_kms = norm(out.flyby.vinf_out);

%% ===== 8. PRE-DSM2 STATE ================================================
S_preDSM2_syn   = S3(:,1).';
S_preDSM2_J2000 = synodic2sun_J2000(S_preDSM2_syn, epoch_dsm2, mu);

out.dsm2.epoch           = epoch_dsm2;
out.dsm2.state_pre.syn   = S_preDSM2_syn;
out.dsm2.state_pre.J2000 = S_preDSM2_J2000;
out.dsm2.dv.J2000_kms    = dv_dsm2_J2000.';
out.dsm2.dv.syn_adim     = dv_J2000_to_syn(dv_dsm2_J2000, epoch_dsm2, mu).';
out.dsm2.dv.norm_ms      = norm(dv_dsm2_J2000) * 1e3;

%% ===== 9. COMET ENCOUNTER ===============================================
S_comet_J2000 = [target_position(:); vel_comet_arr(:)].';           % 1×6
S_comet_syn   = sun_J2000_to_synodic(S_comet_J2000, epoch_comet_flyby, mu);

out.comet.epoch       = epoch_comet_flyby;
out.comet.state.J2000 = S_comet_J2000;
out.comet.state.syn   = S_comet_syn;

%% ===== Print summary ====================================================
fprintf('\n========== VARIABLES ORGANIZER ==========\n');
fprintf('  Departure epoch   : %.6e ET s  (%s)\n', epoch_dep, cspice_et2utc(epoch_dep,'C',0));
fprintf('  DSM1 epoch        : %.6e ET s  (%s)\n', epoch_dsm1, cspice_et2utc(epoch_dsm1,'C',0));
fprintf('  Flyby epoch       : %.6e ET s  (%s)\n', epoch_flyby, cspice_et2utc(epoch_flyby,'C',0));
fprintf('  DSM2 epoch        : %.6e ET s  (%s)\n', epoch_dsm2, cspice_et2utc(epoch_dsm2,'C',0));
fprintf('  Comet epoch       : %.6e ET s  (%s)\n', epoch_comet_flyby, cspice_et2utc(epoch_comet_flyby,'C',0));
fprintf('\n  Total TOF          : %.1f days\n', (epoch_comet_flyby - epoch_dep)/86400);
fprintf('  TOF dep->DSM1      : %.1f days\n', (epoch_dsm1 - epoch_dep)/86400);
fprintf('  TOF DSM1->flyby    : %.1f days\n', (epoch_flyby - epoch_dsm1)/86400);
fprintf('  TOF flyby->DSM2    : %.1f days\n', (epoch_dsm2 - epoch_flyby)/86400);
fprintf('  TOF DSM2->comet    : %.1f days\n', (epoch_comet_flyby - epoch_dsm2)/86400);
fprintf('\n  dV injection      : [%.4f %.4f %.4f] km/s  |  %.1f m/s\n', ...
    dv_inj_J2000, out.injection.norm_ms);
fprintf('  dV DSM1           : [%.4f %.4f %.4f] km/s  |  %.1f m/s\n', ...
    dv_dsm1_J2000, out.dsm1.dv.norm_ms);
fprintf('  dV DSM2           : [%.4f %.4f %.4f] km/s  |  %.1f m/s\n', ...
    dv_dsm2_J2000, out.dsm2.dv.norm_ms);
fprintf('  dV total           : %.1f m/s\n', ...
    out.injection.norm_ms + out.dsm1.dv.norm_ms + out.dsm2.dv.norm_ms);
fprintf('\n  Vinf in  (Moon)   : [%.4f %.4f %.4f] km/s  |  %.4f km/s\n', ...
    out.flyby.vinf_in, out.flyby.vinf_in_norm_kms);
fprintf('  Vinf out (Moon)   : [%.4f %.4f %.4f] km/s  |  %.4f km/s\n', ...
    out.flyby.vinf_out, out.flyby.vinf_out_norm_kms);
fprintf('==========================================\n\n');

end


%% =========================================================================
%  HELPER: convert a J2000 delta-v [km/s] to the non-dimensional synodic frame
%  For an impulse (position unchanged): dv_syn = R'(epoch) * dv_J / V(epoch)
%% =========================================================================
function dv_syn = dv_J2000_to_syn(dv_J2000, epoch, mu)
% dv_J2000_to_syn  Convert a J2000 delta-v vector to the non-dimensional synodic frame.
    ES = cspice_spkezr('EARTH', epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    rE = ES(1:3);   vE = ES(4:6);
    L  = norm(rE);
    h  = cross(rE, vE);
    e3 = h / norm(h);
    e1 = rE / L;
    e2 = cross(e3, e1);
    R  = [e1, e2, e3];
    omega = norm(h) / L^2;
    V     = L * omega;
    dv_syn = R' * dv_J2000(:) / V;   % 3×1
end


%% =========================================================================
%  HELPER: inverse conversion J2000 [km,km/s] -> non-dimensional synodic
%  Exact inverse of synodic2sun_J2000 (includes the Ldot term)
%  S_J2000: N×6 or 1×6   epochs: N×1 or scalar [ET s]
%% =========================================================================
function S_syn = sun_J2000_to_synodic(S_J2000, epochs, mu)
% sun_J2000_to_synodic  Local copy: convert Sun-centred J2000 states to the
% non-dimensional synodic frame (see Auxiliar/sun_J2000_to_synodic.m).
    epochs = epochs(:).';
    N      = numel(epochs);

    ES       = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');
    rE       = ES(1:3,:);
    vE       = ES(4:6,:);

    L_all    = sqrt(sum(rE.^2, 1));
    h_all    = cross(rE, vE);
    h_nrm    = sqrt(sum(h_all.^2, 1));
    om_all   = h_nrm ./ L_all.^2;
    V_all    = L_all .* om_all;
    Ldot_all = dot(rE, vE, 1) ./ L_all;

    e1 = rE ./ L_all;
    e3 = h_all ./ h_nrm;
    e2 = cross(e3, e1);

    S_J2000 = reshape(S_J2000, N, 6);
    S_syn   = zeros(N, 6);

    for i = 1:N
        R   = [e1(:,i), e2(:,i), e3(:,i)];
        r_J = S_J2000(i, 1:3).';
        v_J = S_J2000(i, 4:6).';

        r_syn_adim_mu = R' * r_J / L_all(i);
        r_syn         = r_syn_adim_mu - [mu; 0; 0];

        omega_vec = e3(:,i) * om_all(i);
        v_rot_km  = R' * (v_J - cross(omega_vec, r_J) - r_J * Ldot_all(i)/L_all(i));
        v_syn     = v_rot_km / V_all(i);

        S_syn(i,:) = [r_syn.' v_syn.'];
    end
end
