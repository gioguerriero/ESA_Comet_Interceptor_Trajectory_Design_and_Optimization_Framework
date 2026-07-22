function plot_initial_guess_ephe(x0, c, S_halo, target_position, epoch_comet_flyby, S_halo_orbit_J2000)
% plot_initial_guess_ephe  Visualize the initial-guess trajectory using the
% N-body model. Propagates all arcs exactly as in nonlinear_constraints_ephe
% and produces a 3D Sun-centred ECLIPJ2000 [AU] plot.
%
% Inputs:
%   x0                 - 16-element initial guess
%   c                  - constants struct
%   S_halo             - halo state in synodic non-dim (6x1 or 1x6)
%   target_position    - comet position at arrival [km] (3x1 or 1x3)
%   epoch_comet_flyby  - comet arrival epoch [ET s]
%   S_halo_orbit_J2000 - (optional) full halo orbit Nx6, Sun-centred ECLIPJ2000 [km; km/s]
%
% Outputs:
%   (none) - produces a 3D figure
%
%   x layout (17 variables, with scalings):
%     x(1:3)   – dv_inj   * 1e+3     [m/s scaled]   → / 1e3 = km/s
%     x(4)     – tof_halo2dsm1_d     [days]
%     x(5:7)   – dv_dsm1  * 10       [km/s scaled]  → / 10  = km/s
%     x(8)     – tof_dsm12flyby_d    [days]
%     x(9)     – epoch_flyby / 1e+8  [ET s scaled]  → * 1e7 = ET s
%     x(10:12) – vel_comet_arr       [km/s]
%     x(13)    – tof_dsm22comet_d / 100  [days scaled] → * 100 = days
%     x(14:16) – dv_dsm2  * 10       [km/s scaled]  → / 10  = km/s
%     x(17)    - t_halo              [days]         (shift along the Halo, CR3BP)

AU = 1.495978707e8;   % km/AU

has_halo_orbit = nargin >= 6 && ~isempty(S_halo_orbit_J2000);

%% ===== NORMALIZE S_HALO (synodic non-dim) -> 6×1 column =================
S_halo_syn = S_halo(:);   % synodic non-dim, 6×1

%% ===== UNPACK x0 ========================================================
% All delta-v as 3×1 columns, all TOFs as scalars

dv_inj           = reshape(x0(1:3),   3, 1) / 1e+3;   % [km/s]
tof_inj2dsm1_d   = x0(4);                              % [days]
dv_dsm1          = reshape(x0(5:7),   3, 1) / 10;     % [km/s]
tof_dsm12flyby_d = x0(8);                              % [days]
epoch_flyby      = x0(9) * 1e+8;                       % [ET s]
vel_comet_arr    = reshape(x0(10:12), 3, 1);           % [km/s]
tof_dsm22comet_d = x0(13) * 100;                       % [days]
dv_dsm2          = reshape(x0(14:16), 3, 1) / 10;     % [km/s]
tof_onHalo_d     = x0(17);                             % [days]

tof_inj2dsm1_s   = tof_inj2dsm1_d   * 86400;
tof_dsm12flyby_s = tof_dsm12flyby_d * 86400;
tof_dsm22comet_s = tof_dsm22comet_d * 86400;

tof_flyby2dsm2_s = epoch_comet_flyby - epoch_flyby - tof_dsm22comet_s;
% epoch_flyby fixed; Halo departure forced to initial_epoch
initial_epoch    = epoch_flyby - tof_dsm12flyby_s - tof_inj2dsm1_s;

% --- Shift along the Halo in CR3BP (synodic non-dim) --------------------------
tof_onHalo_ad = tof_onHalo_d * 86400 / c.Tstar;
if abs(tof_onHalo_ad) > 3600 / c.Tstar
    opt_cr3bp = odeset('AbsTol',1e-10,'RelTol',1e-10);
    [~, S_shift_syn] = ode45(@(t,S) CR3BP(t,S, c.mu), ...
                              [0 tof_onHalo_ad], S_halo_syn, opt_cr3bp);
    S_halo_syn_end = S_shift_syn(end,:)';
else
    S_shift_syn    = [];
    S_halo_syn_end = S_halo_syn;
end

% Convert the (post-shift) Halo state from synodic non-dim to J2000
S_halo_ic = synodic2sun_J2000(S_halo_syn_end', initial_epoch, c.mu)';   % 6×1

%% ===== PROPAGATION (identical to nonlinear_constraints_ephe) =============
% Convention: ode45 initial conditions are always 6×1 columns.
% States extracted from the ode45 output (N×6) as S(end,:)' -> 6×1 column.

opt = odeset('AbsTol',1e-9,'RelTol',1e-9);

% ---- Injection (fixed departure from the Halo) -----------------------------
S_inj_pre  = S_halo_ic;                        % 6×1
S_inj_post = S_inj_pre;
S_inj_post(4:6) = S_inj_post(4:6) + dv_inj;    % 3×1 + 3×1 ✓

% ---- Arc 1: injection -> DSM1 ------------------------------------------
[~, S_arc2] = ode45(@(t,S) NBODY_J2000(t,S, initial_epoch, c), ...
                     [0 tof_inj2dsm1_s], S_inj_post, opt);
S_dsm1_pre = S_arc2(end,:)';      % 6×1

% ---- DSM1 ---------------------------------------------------------------
S_dsm1_post = S_dsm1_pre;
S_dsm1_post(4:6) = S_dsm1_post(4:6) + dv_dsm1;

% ---- Arc 2: DSM1 -> flyby ----------------------------------------------
[~, S_arc3] = ode45(@(t,S) NBODY_J2000(t,S, initial_epoch+tof_inj2dsm1_s, c), ...
                     [0 tof_dsm12flyby_s], S_dsm1_post, opt);
S_flyby_pre = S_arc3(end,:)';     % 6×1

% ---- Arc 4 (backward): comet -> DSM2 ----------------------------------
S_comet_ic = [target_position(:); vel_comet_arr];   % 6×1
[~, S_arc4] = ode45(@(t,S) NBODY_J2000(t,S, epoch_comet_flyby, c), ...
                     [0 -tof_dsm22comet_s], S_comet_ic, opt);
S_dsm2_post = S_arc4(end,:)';     % 6×1

% ---- DSM2 ---------------------------------------------------------------
S_dsm2_pre = S_dsm2_post;
S_dsm2_pre(4:6) = S_dsm2_pre(4:6) - dv_dsm2;

% ---- Arc 5 (backward): DSM2 -> flyby -----------------------------------
[~, S_arc5] = ode45(@(t,S) NBODY_J2000(t,S, epoch_comet_flyby-tof_dsm22comet_s, c), ...
                     [0 -tof_flyby2dsm2_s], S_dsm2_pre, opt);
S_flyby_post = S_arc5(end,:)';    % 6×1

%% ===== BODY EPHEMERIDES (for the plot) ===================================

moon_state  = cspice_spkezr('MOON',  epoch_flyby,   'ECLIPJ2000','NONE','SUN');
earth_dep   = cspice_spkezr('EARTH', initial_epoch, 'ECLIPJ2000','NONE','SUN');
earth_flyby = cspice_spkezr('EARTH', epoch_flyby,   'ECLIPJ2000','NONE','SUN');

%% ===== PLOT =============================================================

figure('Name','Initial Guess - NBODY propagation','Color','w');
hold on; grid on; axis equal;

% --- Halo orbit (if provided) -------------------------------------------
if has_halo_orbit
    plot3(S_halo_orbit_J2000(:,1)/AU, S_halo_orbit_J2000(:,2)/AU, S_halo_orbit_J2000(:,3)/AU, ...
          'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'LineStyle',':', ...
          'DisplayName','Halo orbit');
end

% --- Forward arcs (blue) ------------------------------------------------
plot3(S_arc2(:,1)/AU, S_arc2(:,2)/AU, S_arc2(:,3)/AU, ...
      'b-', 'LineWidth',1.5, 'DisplayName','Inj \rightarrow DSM1');
plot3(S_arc3(:,1)/AU, S_arc3(:,2)/AU, S_arc3(:,3)/AU, ...
      'Color',[0 0.4 0.8], 'LineWidth',1.5, 'DisplayName','DSM1 \rightarrow Flyby');

% --- Backward arcs (red, reversed to show flyby->comet) --------
plot3(flipud(S_arc5(:,1))/AU, flipud(S_arc5(:,2))/AU, flipud(S_arc5(:,3))/AU, ...
      'r-', 'LineWidth',1.5, 'DisplayName','Flyby \rightarrow DSM2');
plot3(flipud(S_arc4(:,1))/AU, flipud(S_arc4(:,2))/AU, flipud(S_arc4(:,3))/AU, ...
      'Color',[0.8 0 0], 'LineWidth',1.5, 'DisplayName','DSM2 \rightarrow Comet');

% --- Injection delta-v arrow -----------------------------------------------
dv_dir       = dv_inj / norm(dv_inj);          % direction unit vector [-]
arrow_len_AU = 0.01;                            % visual length [AU]
quiver3(S_inj_pre(1)/AU, S_inj_pre(2)/AU, S_inj_pre(3)/AU, ...
        dv_dir(1)*arrow_len_AU, dv_dir(2)*arrow_len_AU, dv_dir(3)*arrow_len_AU, ...
        0, 'Color',[0 0.75 0], 'LineWidth',2.5, 'MaxHeadSize',0.8, ...
        'DisplayName', sprintf('\\DeltaV_{inj}  (%.1f m/s)', norm(dv_inj)*1e3));

% --- Key points --------------------------------------------------------
plot3(S_halo_ic(1)/AU,  S_halo_ic(2)/AU,  S_halo_ic(3)/AU, ...
      'cs', 'MarkerSize',8,  'MarkerFaceColor','c',          'DisplayName','Halo (departure)');
plot3(S_inj_pre(1)/AU,  S_inj_pre(2)/AU,  S_inj_pre(3)/AU, ...
      'b^', 'MarkerSize',7,  'MarkerFaceColor','b',          'DisplayName','Injection');
plot3(S_dsm1_pre(1)/AU, S_dsm1_pre(2)/AU, S_dsm1_pre(3)/AU, ...
      'bd', 'MarkerSize',8,  'MarkerFaceColor',[0 0.6 1],   'DisplayName','DSM1');
plot3(S_flyby_pre(1)/AU,  S_flyby_pre(2)/AU,  S_flyby_pre(3)/AU, ...
      'mo', 'MarkerSize',9,  'MarkerFaceColor','m',          'DisplayName','Flyby (forward)');
plot3(S_flyby_post(1)/AU, S_flyby_post(2)/AU, S_flyby_post(3)/AU, ...
      'mx', 'MarkerSize',9,  'LineWidth',2,                  'DisplayName','Flyby (backward)');
plot3(S_dsm2_post(1)/AU, S_dsm2_post(2)/AU, S_dsm2_post(3)/AU, ...
      'rd', 'MarkerSize',8,  'MarkerFaceColor',[1 0.4 0],   'DisplayName','DSM2');
plot3(target_position(1)/AU, target_position(2)/AU, target_position(3)/AU, ...
      'rp', 'MarkerSize',12, 'MarkerFaceColor','r',          'DisplayName','Comet (arrival)');

% --- Corpi celesti -------------------------------------------------------
plot3(0, 0, 0, ...
      'yo', 'MarkerSize',14, 'MarkerFaceColor','y',          'DisplayName','Sole');
plot3(earth_dep(1)/AU,   earth_dep(2)/AU,   earth_dep(3)/AU, ...
      'g^', 'MarkerSize',8,  'MarkerFaceColor','g',          'DisplayName','Earth (departure)');
plot3(earth_flyby(1)/AU, earth_flyby(2)/AU, earth_flyby(3)/AU, ...
      'gv', 'MarkerSize',8,  'MarkerFaceColor',[0 0.6 0],   'DisplayName','Terra (flyby)');
plot3(moon_state(1)/AU,  moon_state(2)/AU,  moon_state(3)/AU, ...
      'ko', 'MarkerSize',7,  'MarkerFaceColor',[0.5 0.5 0.5],'DisplayName','Luna (flyby)');

% --- Decorazioni ---------------------------------------------------------
xlabel('X [AU]'); ylabel('Y [AU]'); zlabel('Z [AU]');
title('Initial Guess - NBODY propagation (ECLIPJ2000)');
legend('Location','best','FontSize',8);
view(2);

% --- Flyby gap --------------------------------------------------------
gap_km = norm(S_flyby_pre(1:3) - S_flyby_post(1:3));
fprintf('  Gap posizione al flyby (forward vs backward): %.1f km\n', gap_km);

hold off;
end
