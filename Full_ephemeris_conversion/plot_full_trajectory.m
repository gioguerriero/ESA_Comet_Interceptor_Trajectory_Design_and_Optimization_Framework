function plot_full_trajectory(out2, S_halo, c, comet_name)
% PLOT_FULL_TRAJECTORY  Full mission trajectory in Sun-Earth synodic frame.
%
% Plots the complete trajectory from injection maneuver to comet arrival,
% using full-ephemeris propagation (NBODY_J2000_full_ephe) for all arcs.
% Includes the Halo parking orbit, Moon orbit, and a side data panel.
%
% INPUTS
%   out2        – struct from variables_organizer_refined
%   S_halo      – Halo orbit initial condition [6x1, synodic adim]
%   c           – constants struct (c.mu, c.Tstar)
%   comet_name  – comet name string for the title (e.g. 'C2008A1')
%
% USAGE (from main, after variables_organizer_refined):
%   plot_full_trajectory(out2, S_halo, c, 'C2008A1')

mu = c.mu;

%% --- Propagation setup -----------------------------------------------------
opt_plot = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);
N_pts    = 400;

%% --- Arc 1: Departure (post-inj) → DSM1 ------------------------------------
s0 = out2.departure.state.J2000(:);
s0(4:6) = s0(4:6) + out2.injection.J2000_kms(:);
[t1, S_arc1] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, out2.departure.epoch, c), ...
    linspace(0, out2.dsm1.epoch - out2.departure.epoch, N_pts), s0, opt_plot);
ep_arc1 = out2.departure.epoch + t1';

%% --- Arc 2: DSM1 (post) → TCM1 ---------------------------------------------
s0 = out2.dsm1.state_pre.J2000(:);
s0(4:6) = s0(4:6) + out2.dsm1.dv.J2000_kms(:);
[t2, S_arc2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, out2.dsm1.epoch, c), ...
    linspace(0, out2.tcm1.epoch - out2.dsm1.epoch, N_pts), s0, opt_plot);
ep_arc2 = out2.dsm1.epoch + t2';

%% --- Arc 3: TCM1 (post) → TCM2  [through flyby] ----------------------------
s0 = out2.tcm1.state_post.J2000(:);
[t3, S_arc3] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, out2.tcm1.epoch, c), ...
    linspace(0, out2.tcm2.epoch - out2.tcm1.epoch, N_pts), s0, opt_plot);
ep_arc3 = out2.tcm1.epoch + t3';

%% --- Arc 4: TCM2 (post) → DSM2 ---------------------------------------------
s0 = out2.tcm2.state_post.J2000(:);
[t4, S_arc4] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, out2.tcm2.epoch, c), ...
    linspace(0, out2.dsm2.epoch - out2.tcm2.epoch, N_pts), s0, opt_plot);
ep_arc4 = out2.tcm2.epoch + t4';

%% --- Arc 5: DSM2 (post) → Comet --------------------------------------------
s0 = out2.dsm2.state_pre.J2000(:);
s0(4:6) = s0(4:6) + out2.dsm2.dv.J2000_kms(:);
[t5, S_arc5] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, out2.dsm2.epoch, c), ...
    linspace(0, out2.comet.epoch - out2.dsm2.epoch, N_pts), s0, opt_plot);
ep_arc5 = out2.dsm2.epoch + t5';

%% --- Convert all arcs to synodic adimensional ------------------------------
SYN1 = j2000_to_synodic(S_arc1, ep_arc1, mu);
SYN2 = j2000_to_synodic(S_arc2, ep_arc2, mu);
SYN3 = j2000_to_synodic(S_arc3, ep_arc3, mu);
SYN4 = j2000_to_synodic(S_arc4, ep_arc4, mu);
SYN5 = j2000_to_synodic(S_arc5, ep_arc5, mu);

% Key nodes in synodic
syn_dep      = out2.departure.state.syn;
syn_dsm1_pre = out2.dsm1.state_pre.syn;
syn_tcm1_pre = out2.tcm1.state_pre.syn;
syn_flyby    = out2.flyby.state.syn;
syn_tcm2_pre = out2.tcm2.state_pre.syn;
syn_dsm2_pre = out2.dsm2.state_pre.syn;
syn_comet    = out2.comet.state.syn;

%% --- Halo parking orbit (use pre-computed points directly) -----------------
% S_halo is an N×6 (or 6×N) matrix of synodic adim states on the Halo orbit.
if size(S_halo, 2) == 6
    S_halo_orbit = S_halo;          % N×6
else
    S_halo_orbit = S_halo';         % transpose to N×6
end

%% --- Moon orbit (last ~35 days before flyby) --------------------------------
t_moon_before = 35 * 86400;   % look back 35 days from flyby
N_moon = 200;
t_moon_vec = linspace(out2.flyby.epoch - t_moon_before, out2.flyby.epoch, N_moon);
r_moon_syn = zeros(N_moon, 6);
r_moon_J2000 = zeros(N_moon, 6);   % also store J2000 for distance calc
for i = 1:N_moon
    st = cspice_spkezr('MOON', t_moon_vec(i), 'ECLIPJ2000', 'NONE', 'SUN');
    r_moon_syn(i, :) = j2000_to_synodic(st', t_moon_vec(i), mu);
    r_moon_J2000(i, :) = st';
end

% Moon at flyby epoch (correct ephemeris position)
st_moon_flyby  = cspice_spkezr('MOON', out2.flyby.epoch, 'ECLIPJ2000', 'NONE', 'SUN');
syn_moon_flyby = j2000_to_synodic(st_moon_flyby', out2.flyby.epoch, mu);

%% --- Lagrange points -------------------------------------------------------
rL = (mu/3)^(1/3);
L1 = [1-mu-rL, 0];
L2 = [1-mu+rL, 0];

%% --- TOF and DV data for text panel ----------------------------------------
tof1 = (out2.dsm1.epoch  - out2.departure.epoch) / 86400;
tof2 = (out2.tcm1.epoch  - out2.dsm1.epoch)      / 86400;
tof3 = (out2.flyby.epoch - out2.tcm1.epoch)       / 86400;
tof4 = (out2.tcm2.epoch  - out2.flyby.epoch)      / 86400;
tof5 = (out2.dsm2.epoch  - out2.tcm2.epoch)       / 86400;
tof6 = (out2.comet.epoch - out2.dsm2.epoch)       / 86400;
tof_total = (out2.comet.epoch - out2.departure.epoch) / 86400;

dv_inj  = out2.injection.norm_ms;
dv_dsm1 = out2.dsm1.dv.norm_ms;
dv_tcm1 = out2.tcm1.dv.norm_ms;
dv_tcm2 = out2.tcm2.dv.norm_ms;
dv_dsm2 = out2.dsm2.dv.norm_ms;
dv_total = dv_inj + dv_dsm1 + dv_tcm1 + dv_tcm2 + dv_dsm2;

dep_str   = cspice_et2utc(out2.departure.epoch, 'C', 0);
comet_str = cspice_et2utc(out2.comet.epoch,     'C', 0);

%% --- Compute minimum distances to Earth and Moon ----------------------------
% Collect all trajectory points and epochs
traj_J2000_all = [S_arc1; S_arc2; S_arc3; S_arc4; S_arc5];
traj_syn_all   = [SYN1; SYN2; SYN3; SYN4; SYN5];
epochs_all     = [ep_arc1(:); ep_arc2(:); ep_arc3(:); ep_arc4(:); ep_arc5(:)];

% Min distance to Earth center (in synodic frame, Earth at (1-mu, 0, 0))
r_earth_syn = [1-mu, 0, 0];
dist_to_earth = sqrt(sum((traj_syn_all(:,1:3) - r_earth_syn).^2, 2));
[min_dist_earth_km, idx_min_earth] = min(dist_to_earth);
min_dist_earth_km = min_dist_earth_km * c.Lstar;  % convert from adim to km

% Min distance to Moon center (compute for all trajectory points)
dist_to_moon = zeros(length(epochs_all), 1);
for i = 1:length(epochs_all)
    st_moon = cspice_spkezr('MOON', epochs_all(i), 'ECLIPJ2000', 'NONE', 'SUN');
    r_moon_J2000_pos = st_moon(1:3)';
    r_sc_J2000_pos = traj_J2000_all(i, 1:3);
    dist_to_moon(i) = norm(r_sc_J2000_pos - r_moon_J2000_pos);
end
[min_dist_moon_km, idx_min_moon] = min(dist_to_moon);

% Altitudes above surface (assuming Earth radius 6371 km, Moon radius 1737 km)
r_earth_surface = 6371;
r_moon_surface = 1737;
alt_earth_km = min_dist_earth_km - r_earth_surface;
alt_moon_km  = min_dist_moon_km - r_moon_surface;

txt = { ...
    '— MISSION DATA —', ...
    '', ...
    sprintf('Departure:  %s', dep_str(1:11)), ...
    sprintf('Arrival:    %s', comet_str(1:11)), ...
    '', ...
    '— TIME OF FLIGHT —', ...
    sprintf('Dep → DSM1 :  %6.1f d', tof1), ...
    sprintf('DSM1 → TCM1:  %6.1f d', tof2), ...
    sprintf('TCM1 → flyby: %6.1f d', tof3), ...
    sprintf('Flyby → TCM2: %6.1f d', tof4), ...
    sprintf('TCM2 → DSM2:  %6.1f d', tof5), ...
    sprintf('DSM2 → comet: %6.1f d', tof6), ...
    sprintf('Total TOF:    %6.1f d', tof_total), ...
    '', ...
    '— DELTA-V —', ...
    sprintf('Injection:  %7.1f m/s', dv_inj), ...
    sprintf('DSM1:       %7.1f m/s', dv_dsm1), ...
    sprintf('TCM1:       %7.1f m/s', dv_tcm1), ...
    sprintf('TCM2:       %7.1f m/s', dv_tcm2), ...
    sprintf('DSM2:       %7.1f m/s', dv_dsm2), ...
    sprintf('Total DV:   %7.1f m/s', dv_total), ...
    '', ...
    '— FLYBY —', ...
    sprintf('v_{inf}: %.3f km/s', out2.flyby.vinf_norm_kms), ...
    '', ...
    '— MINIMUM DISTANCES —', ...
    sprintf('Min dist Earth:  %.0f km', min_dist_earth_km), ...
    sprintf('Min alt Earth:   %.0f km', alt_earth_km), ...
    sprintf('Min dist Moon:   %.0f km', min_dist_moon_km), ...
    sprintf('Min alt Moon:    %.0f km', alt_moon_km), ...
    '', ...
    '— COMET —', ...
    sprintf('Pos err: %.2e km', out2.comet.pos_err_km), ...
};

%% --- Figure layout ---------------------------------------------------------
hfig = figure('Name', sprintf('Full Trajectory — %s', comet_name), ...
              'Color', 'w', 'Position', [80, 80, 1350, 720]);

% Main axes: left 72% of figure
ax = axes('Parent', hfig, 'Position', [0.04, 0.07, 0.68, 0.88]);
hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

%% --- Plot arcs -------------------------------------------------------------
plot3(ax, SYN1(:,1), SYN1(:,2), SYN1(:,3), ...
      '-',  'Color', [0.13 0.47 0.71], 'LineWidth', 2.0, ...
      'DisplayName', 'Arc 1: Dep \rightarrow DSM1');
plot3(ax, SYN2(:,1), SYN2(:,2), SYN2(:,3), ...
      '-',  'Color', [0.12 0.63 0.33], 'LineWidth', 2.0, ...
      'DisplayName', 'Arc 2: DSM1 \rightarrow TCM1');
plot3(ax, SYN3(:,1), SYN3(:,2), SYN3(:,3), ...
      '-',  'Color', [0.89 0.47 0.10], 'LineWidth', 2.0, ...
      'DisplayName', 'Arc 3: TCM1 \rightarrow TCM2  (flyby)');
plot3(ax, SYN4(:,1), SYN4(:,2), SYN4(:,3), ...
      '-',  'Color', [0.58 0.22 0.73], 'LineWidth', 2.0, ...
      'DisplayName', 'Arc 4: TCM2 \rightarrow DSM2');
plot3(ax, SYN5(:,1), SYN5(:,2), SYN5(:,3), ...
      '-',  'Color', [0.84 0.15 0.16], 'LineWidth', 2.0, ...
      'DisplayName', 'Arc 5: DSM2 \rightarrow Comet');

% Halo orbit
plot3(ax, S_halo_orbit(:,1), S_halo_orbit(:,2), S_halo_orbit(:,3), ...
      'k-', 'LineWidth', 1.2, 'DisplayName', 'Halo orbit');

% Moon orbit (dashed, last ~35 days before flyby)
plot3(ax, r_moon_syn(:,1), r_moon_syn(:,2), r_moon_syn(:,3), ...
      '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9, 'DisplayName', 'Moon orbit (last 35 d)');

% Moon at flyby epoch
plot3(ax, syn_moon_flyby(1), syn_moon_flyby(2), syn_moon_flyby(3), ...
      'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
      'MarkerSize', 8, 'DisplayName', 'Moon @ flyby');

% Fixed bodies
plot3(ax, -mu, 0, 0, 'o',  'Color', [1 0.8 0], 'MarkerFaceColor', [1 0.8 0], ...
      'MarkerSize', 12, 'DisplayName', 'Sun');
plot3(ax, 1-mu, 0, 0, 'o', 'Color', [0.2 0.6 1], 'MarkerFaceColor', [0.2 0.6 1], ...
      'MarkerSize',  8, 'DisplayName', 'Earth');
plot3(ax, L1(1), 0, 0, 'k+', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'L1');
plot3(ax, L2(1), 0, 0, 'kx', 'MarkerSize', 9, 'LineWidth', 1.5, 'DisplayName', 'L2');

% Maneuver nodes
plot3(ax, syn_dep(1),      syn_dep(2),      syn_dep(3),      ...
      'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Departure');
plot3(ax, syn_dsm1_pre(1), syn_dsm1_pre(2), syn_dsm1_pre(3), ...
      '^', 'Color', [0.12 0.63 0.33], 'MarkerFaceColor', [0.12 0.63 0.33], ...
      'MarkerSize', 8, 'DisplayName', 'DSM1');
plot3(ax, syn_tcm1_pre(1), syn_tcm1_pre(2), syn_tcm1_pre(3), ...
      '^', 'Color', [0.89 0.47 0.10], 'MarkerFaceColor', [0.89 0.47 0.10], ...
      'MarkerSize', 8, 'DisplayName', 'TCM1');
plot3(ax, syn_tcm2_pre(1), syn_tcm2_pre(2), syn_tcm2_pre(3), ...
      '^', 'Color', [0.58 0.22 0.73], 'MarkerFaceColor', [0.58 0.22 0.73], ...
      'MarkerSize', 8, 'DisplayName', 'TCM2');
plot3(ax, syn_dsm2_pre(1), syn_dsm2_pre(2), syn_dsm2_pre(3), ...
      '^', 'Color', [0.84 0.15 0.16], 'MarkerFaceColor', [0.84 0.15 0.16], ...
      'MarkerSize', 8, 'DisplayName', 'DSM2');
plot3(ax, syn_comet(1),    syn_comet(2),    syn_comet(3),    ...
      'p', 'Color', [0.84 0.15 0.16], 'MarkerFaceColor', [0.84 0.15 0.16], ...
      'MarkerSize', 14, 'DisplayName', sprintf('Comet %s', comet_name));

xlabel(ax, 'x [adim]'); ylabel(ax, 'y [adim]'); zlabel(ax, 'z [adim]');
title(ax, sprintf('Full Mission Trajectory — %s  |  Sun-Earth Synodic Frame', comet_name), ...
      'FontSize', 12);
legend(ax, 'show', 'Location', 'best', 'FontSize', 8);

% zlim([-0.004, 0.0055])
view(0, 90)
%xlim([0.5, 1.5])       % centrato su x=1
%ylim([-0.5, 0.5])      % centrato su y=0


%% --- Side text panel -------------------------------------------------------
ax_txt = axes('Parent', hfig, 'Position', [0.74, 0.05, 0.25, 0.92]);
axis(ax_txt, 'off');
text(ax_txt, 0.05, 0.97, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'FontName', 'Courier New', ...
    'FontSize', 9.5, ...
    'Interpreter', 'none');

end


%% =========================================================================
function syn = j2000_to_synodic(S_J2000, epochs, mu)
% Convert N×6 heliocentric J2000 states to synodic adimensional.

    epochs = epochs(:)';
    N = numel(epochs);
    if size(S_J2000, 1) ~= N
        S_J2000 = S_J2000';
    end

    ES   = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');
    rE   = ES(1:3, :);
    vE   = ES(4:6, :);
    L    = sqrt(sum(rE.^2, 1));
    h    = cross(rE, vE);
    hnrm = sqrt(sum(h.^2, 1));
    om   = hnrm ./ L.^2;
    V    = L .* om;
    Ldot = dot(rE, vE, 1) ./ L;

    e1 = rE ./ L;
    e3 = h  ./ hnrm;
    e2 = cross(e3, e1);

    syn = zeros(N, 6);
    for i = 1:N
        R = [e1(:,i), e2(:,i), e3(:,i)];
        r_J = S_J2000(i, 1:3)';
        v_J = S_J2000(i, 4:6)';

        r_rot_dim = R' * r_J;
        r_syn     = r_rot_dim / L(i) - [mu; 0; 0];

        omega_vec = om(i) * e3(:, i);
        v_rot_dim = R' * (v_J - cross(omega_vec, r_J) - r_rot_dim * Ldot(i) / L(i));
        v_syn     = v_rot_dim / V(i);

        syn(i, :) = [r_syn', v_syn'];
    end
end
