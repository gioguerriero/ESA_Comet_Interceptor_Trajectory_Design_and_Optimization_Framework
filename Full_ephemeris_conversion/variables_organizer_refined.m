function out2 = variables_organizer_refined(out, x_pre, x_post, ...
                                             state_tcm, epoch_tcm, t_tcm2, ...
                                             epoch_flyby, comet_pos, epoch_comet, c_const, dv_tcm_in)
% variables_organizer_refined  Reassemble all physical quantities after the
% pre- and post-flyby refinements in the full ephemeris model.
%
% Accepts two formats for x_post (auto-detected from its length):
%
%   x_post [7x1]  - from refinement_post_flyby:
%                     x(1:3) TCM2 DV [km/s]
%                     x(4:6) DSM2 DV [km/s]
%                     x(7)   epoch_dsm2/1e8
%                   In this case dv_tcm_in [3x1] must be provided.
%
%   x_post [12x1] - from refinement_global_post_flyby:
%                     x(1)   BT/1e3  [km/1e3]
%                     x(2)   BR/1e3  [km/1e3]
%                     x(3:5) DSM2 DV [km/s]
%                     x(6)   epoch_dsm2/1e8
%                     x(7:9)   post-TCM2 position / 1e8  [km/1e8]
%                     x(10:12) post-TCM2 velocity / 10   [km/s /10]
%                   dv_tcm_in not required (computed internally from BT, BR).
%
% Inputs:
%   out         - struct from variables_organizer (provides the departure state)
%   x_pre [7x1] - pre-flyby refinement solution:
%                   x(1:3) inj DV [km/s], x(4:6) DSM1 DV [km/s], x(7) epoch_dsm1/1e8
%   x_post      - post-flyby solution vector (7 or 12 elements, see above)
%   state_tcm   - J2000 state at TCM1 epoch [km, km/s] (pre-maneuver)
%   epoch_tcm   - SPICE ET of TCM1 [s]
%   t_tcm2      - days after flyby when TCM2 is applied [days]
%   epoch_flyby - SPICE ET of flyby [s]
%   comet_pos   - target comet position [km] (3x1 or 1x3)
%   epoch_comet - SPICE ET of comet rendezvous [s]
%   c_const     - constants struct (mu, ...)
%   dv_tcm_in   - [3x1] TCM1 DV [km/s], required if x_post is [7x1]
%
% Outputs:
%   out2 - struct with fields:
%     .departure   epoch, state (.syn, .J2000)
%     .injection   dv (.J2000_kms, .syn_adim, .norm_ms)
%     .dsm1        epoch, state_pre (.syn, .J2000), dv (...)
%     .tcm1        epoch, state_pre (.syn, .J2000), state_post (.syn, .J2000), dv (...)
%     .flyby       epoch, state (.syn, .J2000), moon_state, vinf_kms
%     .tcm2        epoch, state_pre (.syn, .J2000), state_post (.syn, .J2000), dv (...)
%     .dsm2        epoch, state_pre (.syn, .J2000), dv (...)
%     .comet       epoch, state (.syn, .J2000)

if nargin < 11
    dv_tcm_in = [];
end

mu = c_const.mu;

x_pre  = x_pre(:);
x_post = x_post(:);

%% --- Unpack design variables -----------------------------------------------

% Pre-flyby refinement (common to both cases)
%   x_pre(7) = DSM1 epoch variation [days] relative to out.dsm1.epoch
%   (consistent with NC_refinement_pre_flyby: epoch_dsm1 = x(7)*86400 + ref)
dv_inj     = x_pre(1:3);            % [km/s]
dv_dsm1    = x_pre(4:6);            % [km/s]
epoch_dsm1 = x_pre(7) * 86400 + out.dsm1.epoch;   % [ET s]

epoch_tcm2 = epoch_flyby + t_tcm2 * 86400;       % [ET s]

if length(x_post) == 7
    % --- refinement_post_flyby format ---
    dv_tcm2_vec     = x_post(1:3);               % explicit TCM2 DV [km/s]
    dv_dsm2         = x_post(4:6);               % [km/s]
    epoch_dsm2      = x_post(7) * 1e8;           % [ET s]
    use_control_pt  = false;

    if isempty(dv_tcm_in)
        error('variables_organizer_refined: dv_tcm_in is required when x_post is [7x1]');
    end
    dv_tcm = dv_tcm_in(:);

else
    % --- refinement_global_post_flyby format (12 elements, match-point) ---
    %   x(1)    BT / 5000                         [km/5000]
    %   x(2)    BR / 5000                         [km/5000]
    %   x(3:5)  DSM2 DV / 0.1                      [km/s /0.1]
    %   x(6)    DSM2 epoch variation [days] relative to out.dsm2.epoch
    %   x(7:9)  comet arrival velocity / 10        [km/s /10]  (used only by
    %           the optimizer's backward branch; not needed here because the
    %           reconstruction is forward and the comet is the target).
    %   x(10:12) TCM2 DV / 0.1                     [km/s /0.1]
    BT              = x_post(1)    * 5000;        % [km]
    BR              = x_post(2)    * 5000;        % [km]
    dv_dsm2         = x_post(3:5)  * 0.1;         % [km/s]
    epoch_dsm2      = x_post(6) * 86400 + out.dsm2.epoch;   % [ET s]
    dv_tcm2_vec     = x_post(10:12) * 0.1;        % explicit TCM2 DV [km/s]
    use_control_pt  = false;                      % TCM2 applied as a delta

    opt_bplane.verbose = 0;
    [dv_tcm, ~, ~] = bplane_tcm(state_tcm(:), epoch_tcm, epoch_flyby, ...
                                 'MOON', BT, BR, c_const.muMoon, c_const, opt_bplane);
    dv_tcm = dv_tcm(:);
end

%% --- Propagation setup -----------------------------------------------------
opt = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

state_dep = out.departure.state.J2000(:);   % [km; km/s]
epoch_dep = out.departure.epoch;

%% --- 1. DEPARTURE ----------------------------------------------------------
out2.departure.epoch      = epoch_dep;
out2.departure.state.J2000 = out.departure.state.J2000;
out2.departure.state.syn   = out.departure.state.syn;

%% --- 2. INJECTION ----------------------------------------------------------
out2.injection.J2000_kms = dv_inj.';
out2.injection.syn_adim  = dv_J2000_to_syn(dv_inj, epoch_dep, mu).';
out2.injection.norm_ms   = norm(dv_inj) * 1e3;

%% --- 3. ARC: departure → DSM1 (post-injection) -----------------------------
s0 = state_dep;
s0(4:6) = s0(4:6) + dv_inj;

tof_dep2dsm1 = epoch_dsm1 - epoch_dep;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_dep, c_const), ...
               [0, tof_dep2dsm1], s0, opt);

S_pre_dsm1  = S(end, :);
S_post_dsm1 = S_pre_dsm1(:);
S_post_dsm1(4:6) = S_post_dsm1(4:6) + dv_dsm1;

out2.dsm1.epoch            = epoch_dsm1;
out2.dsm1.state_pre.J2000  = S_pre_dsm1;
out2.dsm1.state_pre.syn    = sun_J2000_to_synodic(S_pre_dsm1, epoch_dsm1, mu);
out2.dsm1.dv.J2000_kms     = dv_dsm1.';
out2.dsm1.dv.syn_adim      = dv_J2000_to_syn(dv_dsm1, epoch_dsm1, mu).';
out2.dsm1.dv.norm_ms       = norm(dv_dsm1) * 1e3;

%% --- 4. ARC: DSM1 → TCM1 (pre-flyby chain, for state reporting only) -------
tof_dsm12tcm1 = epoch_tcm - epoch_dsm1;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_dsm1, c_const), ...
               [0, tof_dsm12tcm1], S_post_dsm1, opt);

S_pre_tcm1_propagated = S(end, :);   % state at TCM1 from pre-flyby chain

% TCM1: use state_tcm (the exact input to refinement_global_post_flyby) as reference.
% state_tcm was computed by back-propagating out.flyby.state_pre from epoch_flyby;
% it is the authoritative pre-TCM1 state for the post-flyby chain.
s0_post = state_tcm(:);
s0_post(4:6) = s0_post(4:6) + dv_tcm;   % apply TCM1

out2.tcm1.epoch            = epoch_tcm;
out2.tcm1.state_pre.J2000  = state_tcm.';
out2.tcm1.state_pre.syn    = sun_J2000_to_synodic(state_tcm.', epoch_tcm, mu);
out2.tcm1.state_post.J2000 = s0_post.';
out2.tcm1.state_post.syn   = sun_J2000_to_synodic(s0_post.', epoch_tcm, mu);
out2.tcm1.dv.J2000_kms     = dv_tcm.';
out2.tcm1.dv.syn_adim      = dv_J2000_to_syn(dv_tcm, epoch_tcm, mu).';
out2.tcm1.dv.norm_ms       = norm(dv_tcm) * 1e3;
% Mismatch between pre-flyby propagation and state_tcm (diagnostic)
out2.tcm1.pre_chain_pos_err_km = norm(S_pre_tcm1_propagated(1:3) - state_tcm(1:3).');

%% --- 5. ARC: TCM1 → flyby --------------------------------------------------
% Start from state_tcm + dv_tcm (same as NC_refinement_global_post_flyby).
tof_tcm12flyby = epoch_flyby - epoch_tcm;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm, c_const), ...
               [0, tof_tcm12flyby], s0_post, opt);

S_flyby = S(end, :);
moon_state = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN').';

out2.flyby.epoch       = epoch_flyby;
out2.flyby.state.J2000 = S_flyby;
out2.flyby.state.syn   = sun_J2000_to_synodic(S_flyby, epoch_flyby, mu);
out2.flyby.moon_state  = moon_state;
out2.flyby.vinf_kms    = S_flyby(4:6) - moon_state(4:6);
out2.flyby.vinf_norm_kms = norm(out2.flyby.vinf_kms);

%% --- 6. ARC: TCM1 → TCM2 (through flyby) -----------------------------------
tof_tcm12tcm2 = epoch_tcm2 - epoch_tcm;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm, c_const), ...
               [0, tof_tcm12tcm2], s0_post, opt);

S_pre_tcm2 = S(end, :);   % 1×6, pre-TCM2 state from the continuous propagation

if use_control_pt
    % x_post [12]: post-TCM2 velocity from the design variable (control point)
    dv_tcm2     = state_post_tcm2(4:6) - S_pre_tcm2(4:6)';   % 3×1 [km/s]
    S_post_tcm2 = [S_pre_tcm2(1:3)'; state_post_tcm2(4:6)];  % 6×1
else
    % x_post [7]: explicit TCM2 DV
    dv_tcm2     = dv_tcm2_vec(:);                              % 3×1 [km/s]
    S_post_tcm2 = [S_pre_tcm2(1:3)'; S_pre_tcm2(4:6)' + dv_tcm2];  % 6×1
end

out2.tcm2.epoch            = epoch_tcm2;
out2.tcm2.state_pre.J2000  = S_pre_tcm2;
out2.tcm2.state_pre.syn    = sun_J2000_to_synodic(S_pre_tcm2, epoch_tcm2, mu);
out2.tcm2.state_post.J2000 = S_post_tcm2.';
out2.tcm2.state_post.syn   = sun_J2000_to_synodic(S_post_tcm2.', epoch_tcm2, mu);
out2.tcm2.dv.J2000_kms     = dv_tcm2.';
out2.tcm2.dv.syn_adim      = dv_J2000_to_syn(dv_tcm2, epoch_tcm2, mu).';
out2.tcm2.dv.norm_ms       = norm(dv_tcm2) * 1e3;

%% --- 7. ARC: TCM2 → DSM2 ---------------------------------------------------
tof_tcm22dsm2 = epoch_dsm2 - epoch_tcm2;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm2, c_const), ...
               [0, tof_tcm22dsm2], S_post_tcm2, opt);

S_pre_dsm2  = S(end, :);
S_post_dsm2 = S_pre_dsm2(:);
S_post_dsm2(4:6) = S_post_dsm2(4:6) + dv_dsm2;

out2.dsm2.epoch            = epoch_dsm2;
out2.dsm2.state_pre.J2000  = S_pre_dsm2;
out2.dsm2.state_pre.syn    = sun_J2000_to_synodic(S_pre_dsm2, epoch_dsm2, mu);
out2.dsm2.dv.J2000_kms     = dv_dsm2.';
out2.dsm2.dv.syn_adim      = dv_J2000_to_syn(dv_dsm2, epoch_dsm2, mu).';
out2.dsm2.dv.norm_ms       = norm(dv_dsm2) * 1e3;

%% --- 8. ARC: DSM2 → comet --------------------------------------------------
tof_dsm22comet = epoch_comet - epoch_dsm2;
[~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_dsm2, c_const), ...
               [0, tof_dsm22comet], S_post_dsm2, opt);

S_comet = S(end, :);

out2.comet.epoch       = epoch_comet;
out2.comet.state.J2000 = S_comet;
out2.comet.state.syn   = sun_J2000_to_synodic(S_comet, epoch_comet, mu);
out2.comet.target_pos  = comet_pos(:).';
out2.comet.pos_err_km  = norm(S_comet(1:3) - comet_pos(:).');

%% --- Summary ---------------------------------------------------------------
dv_total_ms = out2.injection.norm_ms + out2.dsm1.dv.norm_ms + out2.tcm1.dv.norm_ms ...
            + out2.tcm2.dv.norm_ms   + out2.dsm2.dv.norm_ms;

fprintf('\n========== VARIABLES ORGANIZER REFINED ==========\n');
fprintf('  Departure epoch : %s\n', cspice_et2utc(epoch_dep,    'C', 0));
fprintf('  DSM1 epoch      : %s  (TOF dep→DSM1  = %.1f d)\n', ...
        cspice_et2utc(epoch_dsm1,  'C', 0), tof_dep2dsm1/86400);
fprintf('  TCM1 epoch      : %s  (TOF DSM1→TCM1 = %.1f d)\n', ...
        cspice_et2utc(epoch_tcm,   'C', 0), tof_dsm12tcm1/86400);
fprintf('  Flyby epoch     : %s  (TOF TCM1→flyby= %.1f d)\n', ...
        cspice_et2utc(epoch_flyby, 'C', 0), tof_tcm12flyby/86400);
fprintf('  TCM2 epoch      : %s  (TOF flyby→TCM2= %.1f d)\n', ...
        cspice_et2utc(epoch_tcm2,  'C', 0), t_tcm2);
fprintf('  DSM2 epoch      : %s  (TOF TCM2→DSM2 = %.1f d)\n', ...
        cspice_et2utc(epoch_dsm2,  'C', 0), tof_tcm22dsm2/86400);
fprintf('  Comet epoch     : %s  (TOF DSM2→comet= %.1f d)\n', ...
        cspice_et2utc(epoch_comet, 'C', 0), tof_dsm22comet/86400);
fprintf('\n  DV injection : %.2f m/s\n', out2.injection.norm_ms);
fprintf('  DV DSM1      : %.2f m/s\n', out2.dsm1.dv.norm_ms);
fprintf('  DV TCM1      : %.2f m/s\n', out2.tcm1.dv.norm_ms);
fprintf('  DV TCM2      : %.2f m/s\n', out2.tcm2.dv.norm_ms);
fprintf('  DV DSM2      : %.2f m/s\n', out2.dsm2.dv.norm_ms);
fprintf('  DV TOTAL     : %.2f m/s\n', dv_total_ms);
fprintf('\n  Vinf @ flyby : %.4f km/s\n', out2.flyby.vinf_norm_kms);
fprintf('  Comet pos err: %.4e km\n', out2.comet.pos_err_km);
fprintf('==================================================\n\n');

end


%% =========================================================================
function dv_syn = dv_J2000_to_syn(dv_J2000, epoch, mu)
% dv_J2000_to_syn  Convert a J2000 delta-v vector to the non-dimensional synodic frame.
    ES = cspice_spkezr('EARTH', epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    rE = ES(1:3);  vE = ES(4:6);
    L  = norm(rE);
    h  = cross(rE, vE);
    e3 = h / norm(h);
    e1 = rE / L;
    e2 = cross(e3, e1);
    R  = [e1, e2, e3];
    V  = L * norm(h) / L^2;
    dv_syn = R' * dv_J2000(:) / V;
end


%% =========================================================================
function S_syn = sun_J2000_to_synodic(S_J2000, epochs, mu)
% sun_J2000_to_synodic  Local copy: convert Sun-centred J2000 states to the
% non-dimensional synodic frame (see Auxiliar/sun_J2000_to_synodic.m).
    epochs = epochs(:).';
    N      = numel(epochs);

    ES       = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');
    rE       = ES(1:3, :);
    vE       = ES(4:6, :);
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

        r_rot_dim = R' * r_J;
        r_syn     = r_rot_dim / L_all(i) - [mu; 0; 0];

        omega_vec = e3(:,i) * om_all(i);
        v_rot_dim = R' * (v_J - cross(omega_vec, r_J) - r_J * Ldot_all(i) / L_all(i));
        v_syn     = v_rot_dim / V_all(i);

        S_syn(i, :) = [r_syn.', v_syn.'];
    end
end
