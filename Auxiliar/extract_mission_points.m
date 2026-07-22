function out = extract_mission_points(x_opt, S_halo_syn, target_position, ...
                                      epoch_comet_flyby, c)
% extract_mission_points  Extract key states (synodic non-dim) and TOFs from the
% ephemeris-optimized solution.
%
% Inputs:
%   x_opt             - 17x1 optimized design vector (with scalings)
%                       x = [dv_inj*1e3; tof_h2d1; dv_dsm1*10; tof_d1f;
%                            epoch_flyby/1e7; vel_comet; tof_d2c/100; dv_dsm2*10; t_halo]
%   S_halo_syn        - 1x6 or 6x1 halo state in SYNODIC NON-DIM
%   target_position   - 3x1 comet position [km, J2000]
%   epoch_comet_flyby - ET [s] of comet arrival
%   c                 - constants struct (c.mu, c.Lstar, c.Tstar, c.Vstar)
%
% Outputs:
%   out - struct with fields:
%     .S_inj_pre      pre-injection state (post-Halo wait, pre-dv_inj)   [synodic non-dim]
%     .S_inj          injection state (post-dv_inj)                      [synodic non-dim]
%     .S_dsm1         state after DSM1                                   [synodic non-dim]
%     .S_flyby_m1d    state 1 day before the flyby                       [synodic non-dim]
%     .S_flyby        state at the flyby (forward)                       [synodic non-dim]
%     .S_flyby_p16d   state 16 days after the flyby                      [synodic non-dim]
%     .S_dsm2         state after DSM2                                   [synodic non-dim]
%     .dv_inj_km, .dv_dsm1_km, .dv_dsm2_km   [km/s] delta-v magnitudes
%     .tof_total_d   [days] total mission time
%     .tof_inj2dsm1_ad
%     .tof_dsm12flybyM1d_ad   (DSM1 -> 1 day before the flyby)
%     .tof_flyby2dsm2_ad
%     .tof_dsm22comet_ad
%     .date_comet_flyby       UTC string (YYYY MON DD HH:MM:SS.sss)

%% ---- Unpack x_opt (17 variables with scalings) ---------------------------
dv_inj           = reshape(x_opt(1:3),   3, 1) / 1e+3;   % [km/s]
tof_inj2dsm1_d   = x_opt(4);                              % [days]
dv_dsm1          = reshape(x_opt(5:7),   3, 1) / 10;     % [km/s]
tof_dsm12flyby_d = x_opt(8);                              % [days]
epoch_flyby      = x_opt(9) * 1e+8;                       % [ET s]
vel_comet_arr    = reshape(x_opt(10:12), 3, 1);           % [km/s]
tof_dsm22comet_d = x_opt(13) * 100;                       % [days]
dv_dsm2          = reshape(x_opt(14:16), 3, 1) / 10;     % [km/s]
tof_onHalo_d     = x_opt(17);                             % [days]

day = 86400;
tof_inj2dsm1_s   = tof_inj2dsm1_d   * day;
tof_dsm12flyby_s = tof_dsm12flyby_d * day;
tof_dsm22comet_s = tof_dsm22comet_d * day;
tof_onHalo_s     = tof_onHalo_d     * day;
tof_flyby2dsm2_s = epoch_comet_flyby - epoch_flyby - tof_dsm22comet_s;
initial_epoch    = epoch_flyby - tof_dsm12flyby_s - tof_inj2dsm1_s;

%% ---- Arc propagation (J2000) ----------------------------------------
opt           = odeset('AbsTol',1e-10,'RelTol',1e-10);
S_halo_syn_0  = S_halo_syn(:)';   % force 1×6

% Shift along the Halo in CR3BP (synodic non-dim) if tof_onHalo is significant
tof_onHalo_ad = tof_onHalo_s / c.Tstar;
if abs(tof_onHalo_ad) > 3600 / c.Tstar
    opt_cr3bp = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [~, S_shift_syn] = ode45(@(t,S) CR3BP(t,S, c.mu), [0 tof_onHalo_ad], S_halo_syn_0', opt_cr3bp);
    S_halo_syn_end = S_shift_syn(end,:);   % 1×6
else
    S_halo_syn_end = S_halo_syn_0;         % 1×6
end

% Convert the Halo point (synodic non-dim -> J2000) at the initial_epoch
S_inj_pre = synodic2sun_J2000(S_halo_syn_end, initial_epoch, c.mu)';  % 6×1

% Injection (after the optional Halo shift)
S_inj_post = S_inj_pre;   S_inj_post(4:6) = S_inj_post(4:6) + dv_inj;
ep_inj = initial_epoch;

% Injection -> DSM1
[~, S2] = ode45(@(t,S) NBODY_J2000(t,S,ep_inj,c), [0 tof_inj2dsm1_s], S_inj_post, opt);
S_dsm1_pre = S2(end,:)';

% DSM1
S_dsm1_post = S_dsm1_pre;  S_dsm1_post(4:6) = S_dsm1_post(4:6) + dv_dsm1;
ep_dsm1 = ep_inj + tof_inj2dsm1_s;

% DSM1 -> (1 day before the flyby)
tof_dsm12flybyM1d_s = tof_dsm12flyby_s - day;
[~, S3a] = ode45(@(t,S) NBODY_J2000(t,S,ep_dsm1,c), [0 tof_dsm12flybyM1d_s], S_dsm1_post, opt);
S_flybyM1d = S3a(end,:)';

% (1 day before the flyby) -> flyby
[~, S3b] = ode45(@(t,S) NBODY_J2000(t,S,ep_dsm1+tof_dsm12flybyM1d_s,c), ...
                 [0 day], S_flybyM1d, opt);
S_flyby = S3b(end,:)';     % forward state at the flyby

% flyby -> (flyby + 16 days)   [forward, same post-flyby forward velocity]
[~, S4]  = ode45(@(t,S) NBODY_J2000(t,S,epoch_flyby,c), [0 16*day], S_flyby, opt);
S_flybyP16d = S4(end,:)';

% Backward from the comet -> DSM2
S_comet_ic = [target_position(:); vel_comet_arr];
[~, S5]   = ode45(@(t,S) NBODY_J2000(t,S,epoch_comet_flyby,c), [0 -tof_dsm22comet_s], S_comet_ic, opt);
S_dsm2_post = S5(end,:)';  % state AFTER DSM2 (forward direction toward the comet)

%% ---- Epoch of each point (needed for the synodic conversion) ---
ep_inj_post    = ep_inj;
ep_dsm1_post   = ep_dsm1;
ep_flybyM1d    = ep_dsm1 + tof_dsm12flybyM1d_s;
ep_flyby       = epoch_flyby;
ep_flybyP16d   = epoch_flyby + 16*day;
ep_dsm2_post   = epoch_comet_flyby - tof_dsm22comet_s;

%% ---- Conversion J2000 -> synodic non-dim -------------------------------
out.S_inj_pre    = j2000_to_synodic_adim(S_inj_pre,   ep_inj_post,  c);
out.S_inj        = j2000_to_synodic_adim(S_inj_post,  ep_inj_post,  c);
out.S_dsm1       = j2000_to_synodic_adim(S_dsm1_post, ep_dsm1_post, c);
out.S_flyby_m1d  = j2000_to_synodic_adim(S_flybyM1d,  ep_flybyM1d,  c);
out.S_flyby      = j2000_to_synodic_adim(S_flyby,     ep_flyby,     c);
out.S_flyby_p16d = j2000_to_synodic_adim(S_flybyP16d, ep_flybyP16d, c);
out.S_dsm2       = j2000_to_synodic_adim(S_dsm2_post, ep_dsm2_post, c);

%% ---- TOF [days] --------------------------------------------------------
tof_flyby2dsm2_d = tof_flyby2dsm2_s / day;

out.tof_inj2dsm1_d      = tof_inj2dsm1_d;
out.tof_dsm12flyby_d    = tof_dsm12flyby_d;
out.tof_flyby2dsm2_d    = tof_flyby2dsm2_d;
out.tof_dsm22comet_d    = tof_dsm22comet_d;
out.tof_total_d         = tof_inj2dsm1_d + tof_dsm12flyby_d + tof_flyby2dsm2_d + tof_dsm22comet_d;

% Non-dimensional TOFs (for downstream use)
out.tof_inj2dsm1_ad      = tof_inj2dsm1_s      / c.Tstar;
out.tof_dsm12flybyM1d_ad = (tof_dsm12flyby_s - day) / c.Tstar;
out.tof_flyby2dsm2_ad    = tof_flyby2dsm2_s    / c.Tstar;
out.tof_dsm22comet_ad    = tof_dsm22comet_s    / c.Tstar;

%% ---- ΔV [m/s] ---------------------------------------------------------
out.dv_inj_ms   = norm(dv_inj)  * 1e3;
out.dv_dsm1_ms  = norm(dv_dsm1) * 1e3;
out.dv_dsm2_ms  = norm(dv_dsm2) * 1e3;
out.dv_total_ms = out.dv_inj_ms + out.dv_dsm1_ms + out.dv_dsm2_ms;

%% ---- Comet arrival date -----------------------------------------------
out.date_comet_flyby = cspice_et2utc(epoch_comet_flyby, 'C', 3);

%% ---- High-precision print -----------------------------------------
pretty_print(out);

%% ---- Plot key points in the synodic frame -------------------------------
plot_mission_points(out);

end

% =========================================================================
function pretty_print(o)
fmt_s  = '  %-16s  [% .15f  % .15f  % .15f   % .15f  % .15f  % .15f]\n';
fmt_dv = '  %-16s  %10.4f  m/s\n';
fmt_t  = '  %-16s  %10.4f  days\n';
fprintf('\n========== Mission points (synodic non-dim) ==========\n');
fprintf(fmt_s, 'S_inj_pre',    o.S_inj_pre);
fprintf(fmt_s, 'S_inj',        o.S_inj);
fprintf(fmt_s, 'S_dsm1',       o.S_dsm1);
fprintf(fmt_s, 'S_flyby_m1d',  o.S_flyby_m1d);
fprintf(fmt_s, 'S_flyby',      o.S_flyby);
fprintf(fmt_s, 'S_flyby_p16d', o.S_flyby_p16d);
fprintf(fmt_s, 'S_dsm2',       o.S_dsm2);
fprintf('\n========== ΔV [m/s] ==========\n');
fprintf(fmt_dv, 'dv_inj',   o.dv_inj_ms);
fprintf(fmt_dv, 'dv_dsm1',  o.dv_dsm1_ms);
fprintf(fmt_dv, 'dv_dsm2',  o.dv_dsm2_ms);
fprintf(fmt_dv, 'dv_total', o.dv_total_ms);
fprintf('\n========== TOF [days] ==========\n');
fprintf(fmt_t, 'inj → DSM1',   o.tof_inj2dsm1_d);
fprintf(fmt_t, 'DSM1 → Flyby', o.tof_dsm12flyby_d);
fprintf(fmt_t, 'Flyby → DSM2', o.tof_flyby2dsm2_d);
fprintf(fmt_t, 'DSM2 → Comet', o.tof_dsm22comet_d);
fprintf(fmt_t, 'TOT',          o.tof_total_d);
fprintf('\n========== Date ==========\n');
fprintf('  date_comet_flyby   %s\n\n', o.date_comet_flyby);
end

% =========================================================================
function plot_mission_points(o)
% plot_mission_points  Synodic non-dim figure with all key points,
%                     including the pre-injection state (post-Halo wait).

labels  = {'S\_inj\_pre','S\_inj','S\_dsm1','S\_flyby\_m1d', ...
           'S\_flyby','S\_flyby\_p16d','S\_dsm2'};
states  = {o.S_inj_pre, o.S_inj, o.S_dsm1, o.S_flyby_m1d, ...
           o.S_flyby,   o.S_flyby_p16d, o.S_dsm2};
colors  = {'k','b','g','m','r','c','y'};
markers = {'o','s','d','^','p','h','v'};

figure('Name','Mission points - synodic non-dim','NumberTitle','off');
hold on; grid on; axis equal;

for k = 1:numel(states)
    s = states{k};
    plot3(s(1), s(2), s(3), [colors{k}, markers{k}], ...
          'MarkerSize', 8, 'MarkerFaceColor', colors{k}, 'DisplayName', labels{k});
    text(s(1), s(2), s(3), ['  ' strrep(labels{k},'\_','_')], 'FontSize', 7);
end

% L2 point (1-mu, 0, 0) as a reference
mu_plot = 3.0034e-6;   % Sun-Earth
xL2 = 1 - mu_plot + (mu_plot/3)^(1/3);
plot3(xL2, 0, 0, 'k+', 'MarkerSize', 10, 'DisplayName', 'L2');
text(xL2, 0, 0, '  L2', 'FontSize', 7);

xlabel('x [adim]'); ylabel('y [adim]'); zlabel('z [adim]');
title('Mission points - synodic non-dim frame');
legend('Location','best','FontSize',7);
view(2);
hold off;
end

% =========================================================================
function S_syn = j2000_to_synodic_adim(S_J, et, c)
% j2000_to_synodic_adim  Inverse of synodic2sun_J2000 for a single epoch (with Ldot correction).
% S_J: 6×1 [km; km/s] J2000 Sun-centred. Returns 1×6 synodic non-dim.

ES = cspice_spkezr('EARTH', et, 'ECLIPJ2000', 'NONE', 'SUN');
rE = ES(1:3);   vE = ES(4:6);

L      = norm(rE);
h      = cross(rE, vE);
om     = norm(h) / L^2;
V      = L * om;
Ldot   = dot(rE, vE) / L;

e1 = rE / L;
e3 = h  / norm(h);
e2 = cross(e3, e1);
R  = [e1, e2, e3];

r_J = S_J(1:3);
v_J = S_J(4:6);

r_rot      = R' * r_J;
r_body_dim = r_rot * (Ldot / L);
omega_vec  = e3 * om;
v_rot      = R' * (v_J - cross(omega_vec, r_J)) - r_body_dim;

S_syn = [r_rot(1)/L - c.mu, r_rot(2)/L, r_rot(3)/L, ...
         v_rot(1)/V,         v_rot(2)/V, v_rot(3)/V];
end
