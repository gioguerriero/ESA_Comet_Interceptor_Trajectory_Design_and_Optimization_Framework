function plot_from_state_vector_MS(x0_ms, k_vec, f_nodes, S_traj_syn, c_const, ...
                               target_position, epoch_comet_flyby, S_halo_orbit_syn)
% plot_from_state_vector_MS  Visualize the multiple-shooting solution/guess in
%   the non-dimensional CR3BP Sun-Earth synodic frame. Produces a figure.
%
%   Inputs:
%     x0_ms             - full design-variable vector (17 base + 6*sum(k) node states)
%     k_vec             - [k1 k2 k3 k4] nodes per segment
%     f_nodes           - 1x4 cell, fixed node time fractions
%     S_traj_syn        - (Nx6) reference CR3BP trajectory (synodic non-dim)
%     c_const           - CR3BP constants struct (c_const.mu, c_const.Lstar, ...)
%     target_position   - comet position [km] in J2000
%     epoch_comet_flyby - comet arrival epoch [ET s]
%     S_halo_orbit_syn  - (optional) full halo orbit Nx6 (synodic non-dim)
%
%   x0_ms layout:
%     x(1:3)   dv_inj*1e3      x(4:6)  dv_dsm1*10
%     x(7)     epoch_flyby [days before comet]  x(8:10) vel_comet_arr [km/s]
%     x(11:13) dv_dsm2*10      x(14)   t_halo
%     x(15)    epoch_dep   [days before comet]
%     x(16)    epoch_dsm1  [days before comet]
%     x(17)    epoch_dsm2  [days before comet]
%     then: [states 6×k1] segment 1, [states 6×k2] segment 2,
%          [states 6×k3] segment 3, [states 6×k4] segment 4.
%     Node epochs are reconstructed from f_nodes (T = Ts + f*(Te-Ts)).

has_halo_orb = nargin >= 8 && ~isempty(S_halo_orbit_syn);

mu    = c_const.mu;
Lstar = c_const.Lstar;

%% ===== UNPACK X0 ========================================================
k1 = k_vec(1);  k2 = k_vec(2);
k3 = k_vec(3);  k4 = k_vec(4);

epoch_flyby = epoch_comet_flyby - x0_ms(7)  * 86400;   % [ET s]  (x in days before comet)
epoch_dep   = epoch_comet_flyby - x0_ms(15) * 86400;   % [ET s]  (x in days before comet)
epoch_dsm1  = epoch_comet_flyby - x0_ms(16) * 86400;   % [ET s]
epoch_dsm2  = epoch_comet_flyby - x0_ms(17) * 86400;   % [ET s]

dv_dsm1_km    = x0_ms(4:6).'  ./ 10;   % [km/s] column
vel_comet_arr = x0_ms(8:10).';          % [km/s] column

% Node states (node epochs are NOT in x: reconstructed from f_nodes)
n_base = 17;
off1 = n_base;          off2 = off1 + 6*k1;
off3 = off2 + 6*k2;     off4 = off3 + 6*k3;

S1 = reshape(x0_ms(off1+1 : off1+6*k1),  6, k1);
S2 = reshape(x0_ms(off2+1 : off2+6*k2),  6, k2);
S3 = reshape(x0_ms(off3+1 : off3+6*k3),  6, k3);
S4 = reshape(x0_ms(off4+1 : off4+6*k4),  6, k4);

% Node epochs reconstructed from fixed fractions: T = Ts + f*(Te - Ts)
T1 = epoch_dep   + f_nodes{1}(:) .* (epoch_dsm1        - epoch_dep);
T2 = epoch_dsm1  + f_nodes{2}(:) .* (epoch_flyby       - epoch_dsm1);
T3 = epoch_flyby + f_nodes{3}(:) .* (epoch_dsm2        - epoch_flyby);
T4 = epoch_dsm2  + f_nodes{4}(:) .* (epoch_comet_flyby - epoch_dsm2);

%% ===== CELESTIAL BODY POSITIONS ==========================================
gamma = (mu/3)^(1/3);
L1_x  = 1 - mu - gamma;
L2_x  = 1 - mu + gamma;

r_moon_ad      = 384400 / Lstar;
theta_circ     = linspace(0, 2*pi, 120);
moon_orbit_x   = r_moon_ad * cos(theta_circ) + (1 - mu);
moon_orbit_y   = r_moon_ad * sin(theta_circ);

moon_syn   = j2000_to_synodic_pos(cspice_spkezr('MOON', epoch_flyby, ...
                'ECLIPJ2000','NONE','SUN'), epoch_flyby, mu);
target_syn = j2000_to_synodic_pos(target_position(:), epoch_comet_flyby, mu);

%% ===== ΔV info ===========================================================
dv_inj  = norm(x0_ms(1:3))  / 1e+3 * 1e+3;   % [m/s]
dv_dsm1 = norm(x0_ms(4:6))  / 10   * 1e+3;   % [m/s]
dv_dsm2 = norm(x0_ms(11:13))/ 10   * 1e+3;   % [m/s]
dv_tot  = dv_inj + dv_dsm1 + dv_dsm2;

tof_d = (epoch_comet_flyby - epoch_dep) / 86400;
title_str = sprintf('MS Initial Guess | k=[%d %d %d %d] | TOF=%.0f d | dV_{tot}=%.0f m/s', ...
            k1, k2, k3, k4, tof_d, dv_tot);

%% ===== NBODY PROPAGATION OF THE MS SEGMENTS ===============================
%
% Segment 1 (blue):  interior arcs T1(i)->T1(i+1)
%                  [arc dep->T1(1) unavailable: injection state missing]
%
% Segment 2 (green): bridge post-DSM1 (S1(k1)+dv_dsm1 -> T2(1) or flyby)
%                   + interior arcs T2
%                   + final arc T2(k2)->flyby
%
% Segment 3 (orange): interior arcs T3(i)->T3(i+1)
%                     + final arc T3(k3)->flyby (backward)
%
% Segment 4 (purple): arc comet->T4(1) (backward)
%                    + interior arcs T4
%                    + bridge T4(k4)->T3(1) (backward, visual gap at DSM2)

fprintf('Propagating MS segments with NBODY_J2000...\n');

% --- Segment 1: interior arcs -------------------------------------------
arcs1 = propagate_ms_arc(S1, T1(:), k1, c_const);

% --- Segment 2: post-DSM1 bridge + interior + arc to flyby ---------------
arcs2 = {};

% Bridge: S1(:,k1) + dv_dsm1 -> T2(1) if k2>0, else -> flyby
if k2 > 0
    t_bridge_end = T2(1);
else
    t_bridge_end = epoch_flyby;
end
seg = arc_syn_with_dv(S1(:,k1), T1(k1), dv_dsm1_km, t_bridge_end, c_const);
if ~isempty(seg), arcs2{end+1} = seg; end

% Interior arcs T2
arcs2 = [arcs2, propagate_ms_arc(S2, T2(:), k2, c_const)];

% Final arc T2(k2) -> flyby (only if k2>0)
if k2 > 0
    seg = arc_syn(S2(:,k2), T2(k2), epoch_flyby, c_const);
    if ~isempty(seg), arcs2{end+1} = seg; end
end

% --- Segment 3: interior arcs + final arc T3(k3) -> flyby (backward) --
arcs3 = propagate_ms_arc(S3, T3(:), k3, c_const);

% Final arc T3(k3) -> flyby (backward, dt < 0)
seg = arc_syn(S3(:,k3), T3(k3), epoch_flyby, c_const);
if ~isempty(seg), arcs3{end+1} = seg; end

% --- Segment 4: comet->T4(1) + interior arcs + bridge T4(k4)->T3(1) -----
arcs4 = {};

% Arc comet -> T4(1) if k4>0, else -> T3(1) (backward)
S_comet_J2000 = [target_position(:); vel_comet_arr(:)];
if k4 > 0
    t_back_end = T4(1);
else
    t_back_end = T3(1);
end
seg = arc_J2000(S_comet_J2000, epoch_comet_flyby, t_back_end, c_const);
if ~isempty(seg), arcs4{end+1} = seg; end

% Interior arcs T4
arcs4 = [arcs4, propagate_ms_arc(S4, T4(:), k4, c_const)];

% Bridge T4(k4) -> T3(1) [pre-DSM2, backward] (only if k4>0)
% The visual gap between this arc and S3(:,1) represents dv_dsm2.
if k4 > 0
    seg = arc_syn(S4(:,k4), T4(k4), T3(1), c_const);
    if ~isempty(seg), arcs4{end+1} = seg; end
end

prop_arcs = {arcs1, arcs2, arcs3, arcs4};

fprintf('Propagation complete.\n');

%% ===== FIGURA 1: piano x-y (view 2D) ====================================
figure('Name','MS Initial Guess - Synodic XY','Color','w','Position',[80 80 1100 750]);
hold on; grid on; axis equal;
title(title_str, 'FontSize',10, 'FontWeight','bold');
xlabel('x_{syn} [-]');  ylabel('y_{syn} [-]');

plot_bodies_and_orbit(mu, L1_x, L2_x, moon_orbit_x, moon_orbit_y, ...
                      moon_syn, target_syn, has_halo_orb, S_halo_orbit_syn);
plot_propagated_arcs(prop_arcs);
plot_traj_and_nodes(S_traj_syn, S1, S2, S3, S4, k1, k2, k3, k4);

view(2);
legend('Location','eastoutside','FontSize',8,'NumColumns',1);

%% ===== FIGURA 2: 3D =====================================================
% figure('Name','MS Initial Guess – Sinodico 3D','Color','w','Position',[200 80 1100 750]);
% hold on; grid on; axis equal;
% title(title_str, 'FontSize',10, 'FontWeight','bold');
% xlabel('x_{syn} [-]');  ylabel('y_{syn} [-]');  zlabel('z_{syn} [-]');
% 
% plot_bodies_and_orbit(mu, L1_x, L2_x, moon_orbit_x, moon_orbit_y, ...
%                       moon_syn, target_syn, has_halo_orb, S_halo_orbit_syn);
% plot_propagated_arcs(prop_arcs);
% plot_traj_and_nodes(S_traj_syn, S1, S2, S3, S4, k1, k2, k3, k4);
% 
% view(3);
% legend('Location','eastoutside','FontSize',8,'NumColumns',1);

%% ===== STAMPA INFO EPOCHE ===============================================
fprintf('\n--- Multiple Shooting Initial Guess ---\n');
fprintf('  epoch_dep    : %.6e ET s\n', epoch_dep);
fprintf('  T1 epochs    : '); fprintf('%.4e  ', T1); fprintf('\n');
if k2>0; fprintf('  T2 epochs    : '); fprintf('%.4e  ', T2); fprintf('\n'); end
fprintf('  epoch_flyby  : %.6e ET s\n', epoch_flyby);
fprintf('  T3 epochs    : '); fprintf('%.4e  ', T3); fprintf('\n');
if k4>0; fprintf('  T4 epochs    : '); fprintf('%.4e  ', T4); fprintf('\n'); end
fprintf('  epoch_comet  : %.6e ET s\n', epoch_comet_flyby);
fprintf('  TOF totale   : %.1f giorni\n', tof_d);
fprintf('  ΔV inj/DSM1/DSM2 = %.1f / %.1f / %.1f m/s\n', dv_inj, dv_dsm1, dv_dsm2);
fprintf('  ΔV totale    : %.1f m/s\n\n', dv_tot);
fprintf('  NOTE: the arc epoch_dep -> T1(1) is not shown\n');
fprintf('        (stato injection non disponibile nel plot).\n\n');

end


%% =========================================================================
%  HELPER: propagate arcs between CONSECUTIVE nodes of a segment (interior)
%  S_nodes 6×k, T_epochs k×1 (crescenti o decrescenti)
%  Returns cell{k-1} of Np×6 in synodic non-dim
%% =========================================================================
function arcs = propagate_ms_arc(S_nodes, T_epochs, k, c_const)
    arcs = {};
    if k < 2, return; end
    mu = c_const.mu;
    for i = 1:k-1
        seg = arc_syn(S_nodes(:,i), T_epochs(i), T_epochs(i+1), c_const);
        if ~isempty(seg), arcs{end+1} = seg; end  %#ok<AGROW>
    end
end


%% =========================================================================
%  HELPER: propagate from a synodic non-dim state to an arrival epoch
%  (funziona sia forward che backward secondo il segno di dt)
%% =========================================================================
function seg_syn = arc_syn(S_start_syn, T_start, T_end, c_const)
    mu  = c_const.mu;
    opt = odeset('AbsTol',1e-12,'RelTol',1e-12,'MaxStep',3600);

    S0 = synodic2sun_J2000(S_start_syn(:).', T_start, mu);
    S0 = S0(:);
    dt = T_end - T_start;

    [t_prop, S_prop] = ode45(@(t,S) NBODY_J2000(t,S,T_start,c_const), [0 dt], S0, opt);
    seg_syn = sun_J2000_to_synodic(S_prop, T_start + t_prop, mu);
end


%% =========================================================================
%  HELPER: like arc_syn but applies a dv [km/s, column] in J2000 at the start
%  (usato per l'arco post-DSM1)
%% =========================================================================
function seg_syn = arc_syn_with_dv(S_start_syn, T_start, dv_km, T_end, c_const)
    mu  = c_const.mu;
    opt = odeset('AbsTol',1e-12,'RelTol',1e-12,'MaxStep',3600);

    S0 = synodic2sun_J2000(S_start_syn(:).', T_start, mu);
    S0 = S0(:);
    S0(4:6) = S0(4:6) + dv_km(:);
    dt = T_end - T_start;

    [t_prop, S_prop] = ode45(@(t,S) NBODY_J2000(t,S,T_start,c_const), [0 dt], S0, opt);
    seg_syn = sun_J2000_to_synodic(S_prop, T_start + t_prop, mu);
end


%% =========================================================================
%  HELPER: propagate from a J2000 state [km;km/s] to an arrival epoch
%  (usato per l'arco dalla cometa, backward)
%% =========================================================================
function seg_syn = arc_J2000(S_start_J2000, T_start, T_end, c_const)
    mu  = c_const.mu;
    opt = odeset('AbsTol',1e-12,'RelTol',1e-12,'MaxStep',3600);

    dt = T_end - T_start;
    [t_prop, S_prop] = ode45(@(t,S) NBODY_J2000(t,S,T_start,c_const), [0 dt], S_start_J2000(:), opt);
    seg_syn = sun_J2000_to_synodic(S_prop, T_start + t_prop, mu);
end


%% =========================================================================
%  HELPER: inverse conversion J2000 [km,km/s] -> synodic non-dim
%  Inverso esatto di synodic2sun_J2000 (incluso termine Ldot)
%  S_J2000: N×6   epochs: N×1 [ET s]
%% =========================================================================
function S_syn = sun_J2000_to_synodic(S_J2000, epochs, mu)
    epochs   = epochs(:).';    % 1×N
    N        = numel(epochs);

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

    S_syn = zeros(N, 6);
    for i = 1:N
        R   = [e1(:,i), e2(:,i), e3(:,i)];
        r_J = S_J2000(i, 1:3).';
        v_J = S_J2000(i, 4:6).';

        % Posizione: r_syn_adim_mu = R'*r_J / L,  r_syn = r_syn_adim_mu - [mu;0;0]
        r_syn_adim_mu = R' * r_J / L_all(i);
        r_syn         = r_syn_adim_mu - [mu; 0; 0];

        % Velocity: v_J = R*v_rot + omega×r_J + r_J*Ldot/L
        % → v_rot = R' * (v_J - omega×r_J - r_J*Ldot/L)
        omega_vec = e3(:,i) * om_all(i);
        v_rot_km  = R' * (v_J - cross(omega_vec, r_J) - r_J * Ldot_all(i)/L_all(i));
        v_syn     = v_rot_km / V_all(i);

        S_syn(i,:) = [r_syn.' v_syn.'];
    end
end


%% =========================================================================
%  HELPER: plot the NBODY-propagated arcs (one legend entry per segment)
%% =========================================================================
function plot_propagated_arcs(prop_arcs)
    col = {[0.05 0.40 0.90], ...   % T1: blu
           [0.05 0.70 0.25], ...   % T2: verde
           [0.90 0.40 0.00], ...   % T3: arancio
           [0.75 0.00 0.75]};      % T4: viola
    arc_labels = {'T1: NBODY (halo→DSM1)', ...
                  'T2: NBODY (DSM1→flyby)', ...
                  'T3: NBODY (flyby→DSM2, rev)', ...
                  'T4: NBODY (DSM2→comet, rev)'};

    for ti = 1:4
        arcs = prop_arcs{ti};
        if isempty(arcs), continue; end
        first_plotted = true;
        for ai = 1:numel(arcs)
            seg = arcs{ai};
            if isempty(seg), continue; end
            if first_plotted
                plot3(seg(:,1), seg(:,2), seg(:,3), '-', ...
                      'Color', col{ti}, 'LineWidth', 1.8, ...
                      'DisplayName', arc_labels{ti});
                first_plotted = false;
            else
                plot3(seg(:,1), seg(:,2), seg(:,3), '-', ...
                      'Color', col{ti}, 'LineWidth', 1.8, ...
                      'HandleVisibility','off');
            end
        end
    end
end


%% =========================================================================
%  HELPER: plot the celestial bodies and the reference CR3BP trajectory
%% =========================================================================
function plot_bodies_and_orbit(mu, L1_x, L2_x, moon_orbit_x, moon_orbit_y, ...
                               moon_syn, target_syn, has_halo_orb, S_halo_orbit_syn)
    plot3(-mu,  0, 0, 'yo', 'MarkerSize',12, 'MarkerFaceColor','y', 'DisplayName','Sun');
    plot3(1-mu, 0, 0, 'go', 'MarkerSize',10, 'MarkerFaceColor','g', 'DisplayName','Earth');
    plot3(L1_x, 0, 0, 'k+', 'MarkerSize',8,  'LineWidth',1.5, 'DisplayName','L_1');
    plot3(L2_x, 0, 0, 'kx', 'MarkerSize',8,  'LineWidth',1.5, 'DisplayName','L_2');

    plot3(moon_orbit_x, moon_orbit_y, zeros(size(moon_orbit_x)), ...
          '--', 'Color',[0.82 0.82 0.82], 'LineWidth',0.8, 'DisplayName','Moon orbit');

    if has_halo_orb
        plot3(S_halo_orbit_syn(:,1), S_halo_orbit_syn(:,2), S_halo_orbit_syn(:,3), ...
              ':', 'Color',[0.55 0.55 0.55], 'LineWidth',0.9, 'DisplayName','Halo orbit');
    end

    plot3(moon_syn(1),   moon_syn(2),   moon_syn(3), ...
          'ks', 'MarkerSize',9, 'MarkerFaceColor',[0.5 0.5 0.5], 'DisplayName','Moon (flyby)');
    plot3(target_syn(1), target_syn(2), target_syn(3), ...
          'rp', 'MarkerSize',14, 'MarkerFaceColor','r', 'DisplayName','Comet (target)');
end


%% =========================================================================
%  HELPER: plot the reference CR3BP trajectory and MS nodes
%% =========================================================================
function plot_traj_and_nodes(S_traj_syn, S1, S2, S3, S4, k1, k2, k3, k4)

    plot3(S_traj_syn(:,1), S_traj_syn(:,2), S_traj_syn(:,3), ...
          '-', 'Color',[0.80 0.80 0.80], 'LineWidth',1.0, 'DisplayName','CR3BP (ref)');

    col = {[0.05 0.40 0.90], ...   % T1: blu
           [0.05 0.70 0.25], ...   % T2: verde
           [0.90 0.40 0.00], ...   % T3: arancio
           [0.75 0.00 0.75]};      % T4: viola

    node_data  = {S1, S2, S3, S4};
    ks         = [k1 k2 k3 k4];
    seg_labels = {'T1: nodi', 'T2: nodi', ...
                  'T3: nodi (rev)', 'T4: nodi (rev)'};

    for ti = 1:4
        Si = node_data{ti};
        if isempty(Si) || ks(ti) == 0, continue; end

        plot3(Si(1,:), Si(2,:), Si(3,:), 'o', ...
              'Color', col{ti}, 'MarkerSize', 7, ...
              'MarkerFaceColor', col{ti}, 'DisplayName', seg_labels{ti});
    end

    % Pre-DSM1 (last T1 node) and pre-DSM2 (first T3 node) highlighted
    if k1 > 0
        plot3(S1(1,k1), S1(2,k1), S1(3,k1), 'd', ...
              'MarkerSize',11, 'Color',[0.0 0.0 0.7], ...
              'MarkerFaceColor',[0.0 0.0 0.7], 'DisplayName','Pre-DSM1');
    end
    if k3 > 0
        plot3(S3(1,1), S3(2,1), S3(3,1), 'd', ...
              'MarkerSize',11, 'Color',[0.7 0.2 0.0], ...
              'MarkerFaceColor',[0.7 0.2 0.0], 'DisplayName','Pre-DSM2');
    end
end


%% =========================================================================
%  HELPER: convert a J2000 position -> synodic non-dim (position only)
%% =========================================================================
function r_syn = j2000_to_synodic_pos(r_J2000, epoch, mu)
    st_e  = cspice_spkezr('EARTH', epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e   = st_e(1:3);   v_e = st_e(4:6);
    L     = norm(r_e);
    x_hat = r_e / L;
    z_hat = cross(r_e, v_e);   z_hat = z_hat / norm(z_hat);
    y_hat = cross(z_hat, x_hat);
    R     = [x_hat, y_hat, z_hat];
    r_syn = (R' * r_J2000(1:3) / L - [mu; 0; 0])';
end
