%% propagate_S_pre_inj.m
%
% Script: propagate S_pre_inj_manual for n_days days with NBODY_J2000, convert
% the trajectory back to the non-dimensional synodic frame, and plot it with
% the Halo orbit, the Sun and the Earth. view(2).
%
% Requires: S_pre_inj_manual in the workspace (1x6, synodic non-dim), or define
%           it in the parameters section below.

%% ───────── PARAMETERS TO EDIT ────────────────────────────────
n_days = 100;         % propagation duration [days]
N_pts  = 3000;        % number of points in the trajectory

% Initial epoch [ET seconds]  - same as used in main
epoch0 = 752530462;

% Uncomment to override S_pre_inj_manual from the workspace:
% S_pre_inj_manual = [x, y, z, vx, vy, vz];
%────────────────────────────────────────────────────────────────────

%% Setup
addpath("Auxiliar/");
startup
cspice_furnsh({'kernels/sckernel.tm'});

%% CR3BP Sun-Earth constants
c.G        = 6.67430e-20;
c.Lstar    = 1.496e+8;
c.Tstar    = 5021870.4424055;
c.Vstar    = 29.806;
c.mSun     = getAstroConstants('Sun',  'mass');
c.mEarth   = getAstroConstants('Earth','mass');
c.mMoon    = 7.34767309e22;
c.rSun     = 695700;
c.rEarth   = 6378;
c.rMoon    = 1737;
c.rMoon_ad = 384400 / c.Lstar;
c.mu       = c.mEarth / (c.mSun + c.mEarth);

%% Check that the initial state exists
if ~exist('S_pre_inj_manual','var')
    error('S_pre_inj_manual not found in the workspace. Define it in the parameters section.');
end
fprintf('Initial state (synodic non-dim):\n');
fprintf('  r = [%+.9f  %+.9f  %+.9f]\n', S_pre_inj_manual(1:3));
fprintf('  v = [%+.9f  %+.9f  %+.9f]\n', S_pre_inj_manual(4:6));

%% 1) Synodic non-dim -> J2000 [km, km/s]
[S_J2000_0, ~] = synodic2sun_J2000(S_pre_inj_manual', epoch0, c.mu);
fprintf('\nInitial J2000 state:\n');
fprintf('  r = [%+.6e  %+.6e  %+.6e] km\n',   S_J2000_0(1:3));
fprintf('  v = [%+.6e  %+.6e  %+.6e] km/s\n', S_J2000_0(4:6));

%% 2) NBODY_J2000 propagation
tof_s = n_days * 86400;
t_vec = linspace(0, tof_s, N_pts);
opt   = odeset('AbsTol',1e-10, 'RelTol',1e-10);

fprintf('\nPropagation in progress (%d days)...\n', n_days);
[t_traj, S_J2000_traj] = ode45(@(t,S) NBODY_J2000(t, S, epoch0, c), ...
    t_vec, S_J2000_0', opt);
fprintf('Done: %d integrated points.\n', length(t_traj));

%% 3) J2000 -> synodic non-dim  (exact inverse of synodic2sun_J2000, including the Ldot correction)
N_traj      = length(t_traj);
epochs_traj = (epoch0 + t_traj)';        % 1×N [ET s]

ES      = cspice_spkezr('EARTH', epochs_traj, 'ECLIPJ2000', 'NONE', 'SUN');
rE      = ES(1:3,:);   vE = ES(4:6,:);

L_all    = sqrt(sum(rE.^2, 1));           % 1×N [km]
h_all    = cross(rE, vE);                 % 3×N
h_nrm    = sqrt(sum(h_all.^2, 1));        % 1×N
om_all   = h_nrm ./ L_all.^2;            % 1×N [rad/s]
V_all    = L_all .* om_all;              % 1×N [km/s]
Ldot_all = dot(rE, vE, 1) ./ L_all;     % 1×N [km/s]  - correction from synodic2sun_J2000

e1 = rE ./ L_all;
e3 = h_all ./ h_nrm;
e2 = cross(e3, e1);

S_syn_traj = zeros(N_traj, 6);
for i = 1:N_traj
    R  = [e1(:,i), e2(:,i), e3(:,i)];
    L  = L_all(i);   V = V_all(i);   om = om_all(i);

    r_J = S_J2000_traj(i,1:3)';
    v_J = S_J2000_traj(i,4:6)';

    % Position
    r_rot_dim = R' * r_J;
    x_syn = r_rot_dim(1)/L - c.mu;
    y_syn = r_rot_dim(2)/L;
    z_syn = r_rot_dim(3)/L;

    % Velocity  (inverse of: v_J = R*v_rot + omega×r_J + R*r_body_dim)
    omega_vec  = e3(:,i) * om;
    r_body_dim = r_rot_dim * (Ldot_all(i) / L);      % Ldot term
    v_rot_dim  = R' * (v_J - cross(omega_vec, r_J)) - r_body_dim;
    vx_syn = v_rot_dim(1)/V;
    vy_syn = v_rot_dim(2)/V;
    vz_syn = v_rot_dim(3)/V;

    S_syn_traj(i,:) = [x_syn, y_syn, z_syn, vx_syn, vy_syn, vz_syn];
end

%% 4) Plot
load('S_halo.mat');   % N×6, synodic non-dim

r_sun   = [-c.mu,   0, 0];
r_earth = [1-c.mu,  0, 0];

figure('Name', sprintf('NBODY - %d days | synodic non-dim', n_days), ...
       'Color','w', 'Units','normalized', 'Position',[0.1 0.1 0.7 0.75]);
hold on; grid on; axis equal;

% Halo orbit
plot3(S_halo(:,1), S_halo(:,2), S_halo(:,3), ...
    'k--', 'LineWidth', 1.2, 'DisplayName', 'Halo orbit (CR3BP)');

% NBODY trajectory (colour-coded by time)
t_days = t_traj / 86400;
patch([S_syn_traj(:,1); NaN], [S_syn_traj(:,2); NaN], [S_syn_traj(:,3); NaN], ...
      [t_days; NaN], 'EdgeColor','interp', 'FaceColor','none', ...
      'LineWidth', 1.8, 'DisplayName', sprintf('NBODY (%d d)', n_days));
cb = colorbar;
cb.Label.String = 'Time [days]';
colormap(gca, jet);

% Initial and final points
scatter3(S_pre_inj_manual(1), S_pre_inj_manual(2), S_pre_inj_manual(3), ...
    120, 'g', 'filled', 'MarkerEdgeColor','k', 'DisplayName', 'Initial point');
scatter3(S_syn_traj(end,1), S_syn_traj(end,2), S_syn_traj(end,3), ...
    120, 'r', 'filled', 'MarkerEdgeColor','k', ...
    'DisplayName', sprintf('Final point (%d d)', n_days));

% Sun
scatter3(r_sun(1), r_sun(2), r_sun(3), 300, [1 0.85 0], 'filled', ...
    'MarkerEdgeColor','k', 'DisplayName', 'Sun');

% Earth
scatter3(r_earth(1), r_earth(2), r_earth(3), 180, [0.2 0.5 1], 'filled', ...
    'MarkerEdgeColor','k', 'DisplayName', 'Earth');

xlabel('x [adim]'); ylabel('y [adim]'); zlabel('z [adim]');
title(sprintf('NBODY\\_J2000 - %d days | Synodic non-dim frame\nEpoch: %s', ...
    n_days, cspice_et2utc(epoch0,'C',0)), 'FontSize', 12);
legend('Location','best', 'FontSize', 9);
view(2);

fprintf('\nFinal synodic state (%d d):\n', n_days);
fprintf('  r = [%+.9f  %+.9f  %+.9f] adim\n', S_syn_traj(end,1:3));
fprintf('  v = [%+.9f  %+.9f  %+.9f] adim\n', S_syn_traj(end,4:6));
