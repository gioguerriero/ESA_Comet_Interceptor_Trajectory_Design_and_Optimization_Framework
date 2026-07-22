function plot_trajectory_heliocentric( ...
    S_traj, t_traj_s, ...
    epoch_dep, epoch_dsm1, epoch_flyby, epoch_dsm2, ...
    S_pre_inj_J2000, S_pre_dsm1_J2000, S_pre_flyby_J2000, S_pre_dsm2_J2000, ...
    selected_comet)
% plot_trajectory_heliocentric  Plot the full mission trajectory in the
% Sun-centred ECLIPJ2000 frame, colour-coded by phase, with event markers.
%
% Inputs:
%   S_traj            - Nx3(+) trajectory states [km]
%   t_traj_s          - Nx1 time from departure for each sample [s]
%   epoch_dep/dsm1/flyby/dsm2 - event epochs [ET s]
%   S_pre_inj_J2000, S_pre_dsm1_J2000, S_pre_flyby_J2000, S_pre_dsm2_J2000 - event states [km,km/s]
%   selected_comet    - comet struct (uses comet_pos for the arrival marker)
%
% Outputs:
%   (none) - produces a 3D figure

figure; hold on; grid on; axis equal
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]')

% =========================
% Trajectory arcs (mission phases)
% =========================
i1 = t_traj_s <= epoch_dsm1 - epoch_dep;
i2 = t_traj_s >  epoch_dsm1 - epoch_dep & t_traj_s <= epoch_flyby - epoch_dep;
i3 = t_traj_s >  epoch_flyby - epoch_dep & t_traj_s <= epoch_dsm2 - epoch_dep;
i4 = t_traj_s >  epoch_dsm2 - epoch_dep;

h1 = plot3(S_traj(i1,1), S_traj(i1,2), S_traj(i1,3), 'Color',[0 0.45 0.74],'LineWidth',1.5);
h2 = plot3(S_traj(i2,1), S_traj(i2,2), S_traj(i2,3), 'Color',[0.85 0.33 0.10],'LineWidth',1.5);
h3 = plot3(S_traj(i3,1), S_traj(i3,2), S_traj(i3,3), 'Color',[0.47 0.67 0.19],'LineWidth',1.5);
h4 = plot3(S_traj(i4,1), S_traj(i4,2), S_traj(i4,3), 'Color',[0.49 0.18 0.56],'LineWidth',1.5);

% =========================
% Events
% =========================
h_inj   = plot3(S_pre_inj_J2000(1),   S_pre_inj_J2000(2),   S_pre_inj_J2000(3),   'o','Color','k','MarkerSize',6,'LineWidth',2);
h_dsm1  = plot3(S_pre_dsm1_J2000(1),  S_pre_dsm1_J2000(2),  S_pre_dsm1_J2000(3),  's','Color','k','MarkerSize',6,'LineWidth',2);
h_flyby = plot3(S_pre_flyby_J2000(1), S_pre_flyby_J2000(2), S_pre_flyby_J2000(3), 'd','Color','k','MarkerSize',6,'LineWidth',2);
h_dsm2  = plot3(S_pre_dsm2_J2000(1),  S_pre_dsm2_J2000(2),  S_pre_dsm2_J2000(3),  '^','Color','k','MarkerSize',6,'LineWidth',2);

% =========================
% Earth and comet
% =========================
state_earth_dep = cspice_spkezr('EARTH', epoch_dep, 'ECLIPJ2000', 'NONE', 'SUN');
h_earth = plot3(state_earth_dep(1), state_earth_dep(2), state_earth_dep(3), ...
    'o','Color',[0 0.45 0.74],'MarkerSize',8,'LineWidth',2);

comet_pos = selected_comet.comet_pos;
h_comet = plot3(comet_pos(1), comet_pos(2), comet_pos(3), ...
    'o','Color',[0.85 0.33 0.10],'MarkerSize',8,'LineWidth',2);

% =========================
% Earth orbit
% =========================
t_orb = linspace(epoch_dep, epoch_dep + 365*86400, 1000);
earth_orb = zeros(length(t_orb),3);
for k = 1:length(t_orb)
    st = cspice_spkezr('EARTH', t_orb(k), 'ECLIPJ2000', 'NONE', 'SUN');
    earth_orb(k,:) = st(1:3)';
end
h_earth_orb = plot3(earth_orb(:,1), earth_orb(:,2), earth_orb(:,3), '--','Color',[0.3 0.3 0.3]);

% =========================
% Moon orbit
% =========================
moon_orb = zeros(length(t_orb),3);
for k = 1:length(t_orb)
    st = cspice_spkezr('MOON', t_orb(k), 'ECLIPJ2000', 'NONE', 'SUN');
    moon_orb(k,:) = st(1:3)';
end
h_moon_orb = plot3(moon_orb(:,1), moon_orb(:,2), moon_orb(:,3), '--','Color',[0 0.75 0.75]);

% =========================
% Sun
% =========================
h_sun = plot3(0,0,0,'o','Color',[1 0.8 0],'MarkerSize',10,'LineWidth',2);

% =========================
% Legend
% =========================
legend([h1 h2 h3 h4 ...
        h_inj h_dsm1 h_flyby h_dsm2 ...
        h_earth h_comet ...
        h_earth_orb h_moon_orb h_sun], ...
       {'Halo departure → DSM1', ...
        'DSM1 → Lunar flyby', ...
        'Flyby → DSM2', ...
        'DSM2 → Comet arrival', ...
        'Injection burn', ...
        'DSM1', ...
        'Lunar flyby', ...
        'DSM2', ...
        'Earth at departure', ...
        'Comet at arrival', ...
        'Earth orbit', ...
        'Moon orbit', ...
        'Sun'}, ...
        'Location','best')

view(2)
end