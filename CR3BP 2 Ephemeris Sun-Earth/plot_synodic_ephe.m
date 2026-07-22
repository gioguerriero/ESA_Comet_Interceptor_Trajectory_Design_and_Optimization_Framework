function plot_synodic_ephe(x0, c, S_halo, target_position, epoch_comet_flyby, x_opt, S_halo_orbit_syn)
% plot_synodic_ephe  Trajectory in the non-dimensional Sun-Earth synodic frame.
%
%   Propagates with NBODY_J2000 and converts to the CR3BP Sun-Earth synodic
%   frame. If x_opt is provided it is overlaid (solid line). If S_halo_orbit_syn
%   is provided the reference Halo orbit is plotted. Produces a figure.
%
%   Inputs:
%     x0                - 17-element initial guess (see layout below)
%     c                 – struct costanti CR3BP (c.mu, c.Tstar, ...)
%     S_halo            - halo state in synodic non-dim (6x1 or 1x6)
%     target_position   - comet position [km] (3x1 or 1x3)
%     epoch_comet_flyby - comet arrival epoch [ET s]
%     x_opt             – (opzionale) 17-element soluzione ottimizzata;
%                         passare [] se non disponibile
%     S_halo_orbit_syn  - (optional) halo orbit Nx6 in the synodic
%                         non-dimensional barycentre-centred frame (from S_halo in main)
%
%   x layout (17 variabili, con scalings):
%     x(1:3)   – dv_inj   * 1e+3     [m/s scaled]   → / 1e3 = km/s
%     x(4)     – tof_halo2dsm1_d     [days]
%     x(5:7)   – dv_dsm1  * 10       [km/s scaled]  → / 10  = km/s
%     x(8)     – tof_dsm12flyby_d    [days]
%     x(9)     – epoch_flyby / 1e+7  [ET s scaled]  → * 1e7 = ET s
%     x(10:12) – vel_comet_arr       [km/s]
%     x(13)    – tof_dsm22comet_d / 100  [days scaled] → * 100 = days
%     x(14:16) – dv_dsm2  * 10       [km/s scaled]  → / 10  = km/s
%     x(17)    – t_halo              [days]         (shift sulla Halo, CR3BP)

mu_cr3bp = c.mu;   % mu_earth / (mu_sun + mu_earth)  [-]

% S_halo is in synodic non-dim -> 6x1 column
S_halo_syn = S_halo(:);

has_opt      = nargin >= 6 && ~isempty(x_opt);
has_halo_orb = nargin >= 7 && ~isempty(S_halo_orbit_syn);

%% ===== PROPAGA INITIAL GUESS ============================================
arcs0 = propagate_arcs(x0, S_halo_syn, target_position, epoch_comet_flyby, c, mu_cr3bp);

% Compute total TOF and delta-v of the initial guess
% tof_total_0 = (x0(4) + x0(8)) + x0(13)*100;  % tof_inj2dsm1 + tof_dsm12flyby + tof_dsm22comet [days]
tof_total_0 = (x0(4) + x0(8) + (epoch_comet_flyby - x0(9)*1e+8) / 86400);
dv_inj_0    = norm(reshape(x0(1:3),3,1) / 1e+3);
dv_dsm1_0   = norm(reshape(x0(5:7),3,1) / 10);
dv_dsm2_0   = norm(reshape(x0(14:16),3,1) / 10);
dv_total_0  = (dv_inj_0 + dv_dsm1_0 + dv_dsm2_0) * 1e3;  % [m/s]

%% ===== PROPAGA SOLUZIONE OTTIMIZZATA (se fornita) =======================
if has_opt
    arcs_opt = propagate_arcs(x_opt, S_halo_syn, target_position, epoch_comet_flyby, c, mu_cr3bp);

    % Calcola TOF e ΔV totali soluzione
    % tof_total_opt = (x_opt(4) + x_opt(8)) + x_opt(13)*100;
    tof_total_opt = (x_opt(4) + x_opt(8) + (epoch_comet_flyby - x_opt(9)*1e+8) / 86400);
    dv_inj_opt    = norm(reshape(x_opt(1:3),3,1) / 1e+3);
    dv_dsm1_opt   = norm(reshape(x_opt(5:7),3,1) / 10);
    dv_dsm2_opt   = norm(reshape(x_opt(14:16),3,1) / 10);
    dv_total_opt  = (dv_inj_opt + dv_dsm1_opt + dv_dsm2_opt) * 1e3;

    epoch_flyby_opt = x_opt(9) * 1e+8;
    st_moon_opt = cspice_spkezr('MOON', epoch_flyby_opt, 'ECLIPJ2000', 'NONE', 'SUN');
    moon_syn_opt = j2000_to_synodic_pos(st_moon_opt(1:3)', epoch_flyby_opt, mu_cr3bp);
else
    tof_total_opt = [];
    dv_total_opt  = [];
    moon_syn_opt  = [];
end

%% ===== CONVERTE IN FRAME SINODICO =======================================
syn0 = arcs_to_synodic(arcs0, mu_cr3bp);
if has_opt
    syn_opt = arcs_to_synodic(arcs_opt, mu_cr3bp);
end

%% ===== POSIZIONI CORPI CELESTI IN FRAME SINODICO ========================
epoch_flyby = x0(9) * 1e+8;
st_moon = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN');
moon_syn = j2000_to_synodic_pos(st_moon(1:3)', epoch_flyby, mu_cr3bp);

% L1 e L2 (approssimazione primo ordine)
gamma = (mu_cr3bp/3)^(1/3);
L1_x  =  1 - mu_cr3bp - gamma;
L2_x  =  1 - mu_cr3bp + gamma;

% Lunar orbit (about 384400 km)
r_moon_ad = 384400 / c.Lstar;
theta_moon_circ = linspace(0, 2*pi, 100);
moon_orbit_x = r_moon_ad * cos(theta_moon_circ) + (1 - mu_cr3bp);
moon_orbit_y = r_moon_ad * sin(theta_moon_circ);
moon_orbit_z = zeros(size(theta_moon_circ));

% Target cometa
target_syn = j2000_to_synodic_pos(target_position, epoch_comet_flyby, mu_cr3bp);

%% ===== FIGURA ============================================================

% Titolo con metriche
if has_opt
    title_str = sprintf('NBODY Synodic | Guess: TOF=%.1f d, dV=%.0f m/s | Opt: TOF=%.1f d, dV=%.0f m/s', ...
                        tof_total_0, dv_total_0, tof_total_opt, dv_total_opt);
else
    title_str = sprintf('NBODY Synodic | Guess: TOF=%.1f d, dV=%.0f m/s', tof_total_0, dv_total_0);
end

figure('Name','Synodic frame – NBODY propagation','Color','w', ...
       'Position',[100 100 1200 800]);
hold on; grid on; axis equal;
title(title_str, 'FontSize', 11, 'FontWeight','bold');
xlabel('x_{syn} [-]');  ylabel('y_{syn} [-]');  zlabel('z_{syn} [-]');

% --- Sun, Earth, L1, L2 -------------------------------------------------
plot3(-mu_cr3bp,   0, 0, 'yo', 'MarkerSize',12, 'MarkerFaceColor','y', 'DisplayName','Sun');
plot3(1-mu_cr3bp,  0, 0, 'go', 'MarkerSize',10, 'MarkerFaceColor','g', 'DisplayName','Earth');
plot3(L1_x,  0, 0, 'k+', 'MarkerSize',8, 'LineWidth',1.5, 'DisplayName','L_1');
plot3(L2_x,  0, 0, 'kx', 'MarkerSize',8, 'LineWidth',1.5, 'DisplayName','L_2');

% --- Halo orbit (if provided) -------------------------------------------
if has_halo_orb
    plot3(S_halo_orbit_syn(:,1), S_halo_orbit_syn(:,2), S_halo_orbit_syn(:,3), ...
          'Color',[0.6 0.6 0.6], 'LineWidth',0.8, 'LineStyle',':', ...
          'DisplayName','Halo orbit');
end

% --- Lunar orbit -------------------------------------------------------
plot3(moon_orbit_x, moon_orbit_y, moon_orbit_z, ...
      'Color',[0.8 0.8 0.8], 'LineWidth',0.8, 'LineStyle','--', ...
      'DisplayName','Moon orbit');

% --- Initial-guess trajectories (dashed) ----------------------------
draw_arcs_3d(syn0, '--', 0.8, 'Guess');

% --- Refined trajectory (solid) ----------------------------------------
if has_opt
    draw_arcs_3d(syn_opt, '-', 1.5, 'Refined');
end

% --- delta-v arrows (refined, if available) --------------------------------
arrow_len = 0.008;
if has_opt
    dv_inj_opt_vec = reshape(x_opt(1:3), 3, 1) / 1e+3;
    dv_dsm1_opt_vec = reshape(x_opt(5:7), 3, 1) / 10;
    dv_dsm2_opt_vec = reshape(x_opt(14:16), 3, 1) / 10;

    % delta-v injection
    dv_dir_inj = dv_j2000_to_synodic_unit(dv_inj_opt_vec, arcs_opt.ep_inj);
    quiver3(syn_opt.inj(1), syn_opt.inj(2), syn_opt.inj(3), ...
            dv_dir_inj(1)*arrow_len, dv_dir_inj(2)*arrow_len, dv_dir_inj(3)*arrow_len, ...
            0, 'Color',[0 0.8 0], 'LineWidth',3, 'MaxHeadSize',1.5, ...
            'DisplayName', sprintf('ΔV_{inj} (%.1f m/s)', norm(dv_inj_opt_vec)*1e3));

    % ΔV DSM1
    dv_dir_dsm1 = dv_j2000_to_synodic_unit(dv_dsm1_opt_vec, arcs_opt.ep_dsm1);
    quiver3(syn_opt.dsm1(1), syn_opt.dsm1(2), syn_opt.dsm1(3), ...
            dv_dir_dsm1(1)*arrow_len, dv_dir_dsm1(2)*arrow_len, dv_dir_dsm1(3)*arrow_len, ...
            0, 'Color',[0.8 0 0.8], 'LineWidth',3, 'MaxHeadSize',1.5, ...
            'DisplayName', sprintf('ΔV_{DSM1} (%.1f m/s)', norm(dv_dsm1_opt_vec)*1e3));

    % ΔV DSM2
    dv_dir_dsm2 = dv_j2000_to_synodic_unit(dv_dsm2_opt_vec, arcs_opt.ep_dsm2);
    quiver3(syn_opt.dsm2(1), syn_opt.dsm2(2), syn_opt.dsm2(3), ...
            dv_dir_dsm2(1)*arrow_len, dv_dir_dsm2(2)*arrow_len, dv_dir_dsm2(3)*arrow_len, ...
            0, 'Color',[1 0.5 0], 'LineWidth',3, 'MaxHeadSize',1.5, ...
            'DisplayName', sprintf('ΔV_{DSM2} (%.1f m/s)', norm(dv_dsm2_opt_vec)*1e3));

    % Moon at the refined flyby
    plot3(moon_syn_opt(1), moon_syn_opt(2), moon_syn_opt(3), ...
          'ko', 'MarkerSize',8, 'MarkerFaceColor',[0.5 0.5 0.5], ...
          'DisplayName','Moon (flyby refined)');
end

% --- Comet target -------------------------------------------------------
plot3(target_syn(1), target_syn(2), target_syn(3), ...
      'rp', 'MarkerSize',14, 'MarkerFaceColor','r', ...
      'DisplayName','Comet (target)');

% --- Moon at the initial-guess flyby ------------------------------------------
plot3(moon_syn(1), moon_syn(2), moon_syn(3), ...
      'ks', 'MarkerSize',8, 'MarkerFaceColor',[0.7 0.7 0.7], ...
      'DisplayName','Moon (flyby guess)');

%nlegend('Location','best','FontSize',8, 'NumColumns',2);
view(2);

% --- Print flyby gap -------------------------------------------------
gap_guess = norm(arcs0.S_flyby_pre(1:3) - arcs0.S_flyby_post(1:3));
fprintf('  Flyby gap - guess: %.1f km\n', gap_guess);
if has_opt
    gap_opt = norm(arcs_opt.S_flyby_pre(1:3) - arcs_opt.S_flyby_post(1:3));
    fprintf('  Flyby gap - opt:   %.1f km\n', gap_opt);
end

end

%% =========================================================================
%  FUNZIONI LOCALI
%% =========================================================================

function arcs = propagate_arcs(x, S_halo_syn, target_position, epoch_comet_flyby, c, mu_cr3bp)
%PROPAGATE_ARCS  Propaga tutti gli archi con NBODY_J2000.
%   Returns a struct with Nx6 trajectories and Nx1 epoch vectors.
%   S_halo_syn is in synodic non-dim; converted to J2000 at the current initial_epoch.

    opt = odeset('AbsTol',1e-9,'RelTol',1e-9);

    % --- Unpack (x a 17 variabili, con scalings) -------------------------
    dv_inj           = reshape(x(1:3),   3, 1) / 1e+3;    % [km/s]
    tof_inj2dsm1_d   = x(4);                               % [days]
    dv_dsm1          = reshape(x(5:7),   3, 1) / 10;      % [km/s]
    tof_dsm12flyby_d = x(8);                               % [days]
    epoch_flyby      = x(9) * 1e+8;                        % [ET s]
    vel_comet_arr    = reshape(x(10:12), 3, 1);            % [km/s]
    tof_dsm22comet_d = x(13) * 100;                        % [days]
    dv_dsm2          = reshape(x(14:16), 3, 1) / 10;      % [km/s]
    tof_onHalo_d     = x(17);                              % [days]

    tof_inj2dsm1_s   = tof_inj2dsm1_d   * 86400;
    tof_dsm12flyby_s = tof_dsm12flyby_d * 86400;
    tof_dsm22comet_s = tof_dsm22comet_d * 86400;
    tof_flyby2dsm2_s = epoch_comet_flyby - epoch_flyby - tof_dsm22comet_s;

    % Partenza dalla Halo a initial_epoch (dopo eventuale shift)
    initial_epoch = epoch_flyby - tof_dsm12flyby_s - tof_inj2dsm1_s;

    % --- Arc 1: shift along the Halo in CR3BP (synodic non-dim) ---------------
    tof_onHalo_ad = tof_onHalo_d * 86400 / c.Tstar;
    if abs(tof_onHalo_ad) > 3600 / c.Tstar
        opt_cr3bp = odeset('AbsTol',1e-10,'RelTol',1e-10);
        [t1, S1_syn] = ode45(@(t,S) CR3BP(t,S, mu_cr3bp), ...
                              [0 tof_onHalo_ad], S_halo_syn, opt_cr3bp);
        S_halo_syn_end = S1_syn(end,:)';
    else
        t1 = [];  S1_syn = [];
        S_halo_syn_end = S_halo_syn;
    end
    S1 = S1_syn;   % in synodic non-dim (flag handled in arcs_to_synodic)

    % Convert the (post-shift) Halo state from synodic non-dim to J2000
    S_halo_ic = synodic2sun_J2000(S_halo_syn_end', initial_epoch, mu_cr3bp)';   % 6×1
    S_inj_pre = S_halo_ic;

    % --- Injection --------------------------------------------------------
    S_inj_post = S_inj_pre;
    S_inj_post(4:6) = S_inj_post(4:6) + dv_inj;

    % --- Arco 2: injection → DSM1 ----------------------------------------
    ep2 = initial_epoch;
    [t2, S2] = ode45(@(t,S) NBODY_J2000(t,S, ep2, c), ...
                     [0 tof_inj2dsm1_s], S_inj_post, opt);
    S_dsm1_pre = S2(end,:)';

    % --- DSM1 -------------------------------------------------------------
    S_dsm1_post = S_dsm1_pre;
    S_dsm1_post(4:6) = S_dsm1_post(4:6) + dv_dsm1;

    % --- Arco 3: DSM1 → flyby --------------------------------------------
    ep3 = initial_epoch + tof_inj2dsm1_s;
    [t3, S3] = ode45(@(t,S) NBODY_J2000(t,S, ep3, c), ...
                     [0 tof_dsm12flyby_s], S_dsm1_post, opt);
    S_flyby_pre = S3(end,:)';

    % --- Arco 4 backward: cometa → DSM2 ----------------------------------
    S_comet_ic = [target_position(:); vel_comet_arr];
    [t4, S4] = ode45(@(t,S) NBODY_J2000(t,S, epoch_comet_flyby, c), ...
                     [0 -tof_dsm22comet_s], S_comet_ic, opt);
    S_dsm2_post = S4(end,:)';

    % --- DSM2 -------------------------------------------------------------
    S_dsm2_pre = S_dsm2_post;
    S_dsm2_pre(4:6) = S_dsm2_pre(4:6) - dv_dsm2;

    % --- Arco 5 backward: DSM2 → flyby -----------------------------------
    ep5 = epoch_comet_flyby - tof_dsm22comet_s;
    [t5, S5] = ode45(@(t,S) NBODY_J2000(t,S, ep5, c), ...
                     [0 -tof_flyby2dsm2_s], S_dsm2_pre, opt);
    S_flyby_post = S5(end,:)';

    % --- Epochs for each point of each arc ------------------------------
    arcs.initial_epoch = initial_epoch;
    arcs.epoch_comet   = epoch_comet_flyby;

    arcs.t1 = t1;  arcs.S1 = S1;  arcs.ep1 = [];   % arc 1 already in synodic non-dim
    arcs.t2 = t2;  arcs.S2 = S2;  arcs.ep2 = ep2 + t2;
    arcs.t3 = t3;  arcs.S3 = S3;  arcs.ep3 = ep3 + t3;
    arcs.t4 = t4;  arcs.S4 = S4;  arcs.ep4 = epoch_comet_flyby + t4;
    arcs.t5 = t5;  arcs.S5 = S5;  arcs.ep5 = ep5 + t5;

    arcs.S_halo_ic    = S_halo_ic;
    arcs.S_inj_pre    = S_inj_pre;
    arcs.S_dsm1_pre   = S_dsm1_pre;
    arcs.S_flyby_pre  = S_flyby_pre;
    arcs.S_flyby_post = S_flyby_post;
    arcs.S_dsm2_post  = S_dsm2_post;
    arcs.S_comet      = S_comet_ic;

    % Epoche dei punti chiave
    arcs.ep_halo  = initial_epoch;
    arcs.ep_inj   = initial_epoch;          % injection coincide col punto Halo
    arcs.ep_dsm1  = ep3;
    arcs.ep_flyby = epoch_flyby;
    arcs.ep_dsm2  = ep5;
    arcs.ep_comet = epoch_comet_flyby;
end


function syn = arcs_to_synodic(arcs, mu)
% arcs_to_synodic  Convert all trajectories to the synodic non-dim frame.

    fields = {'1','2','3','4','5'};
    for k = 1:numel(fields)
        f  = fields{k};
        Sf = arcs.(['S' f]);
        ef = arcs.(['ep' f]);
        if isempty(Sf)
            syn.(['r' f]) = [];
        elseif k == 1
            % Arc 1 is already in synodic non-dim (barycentre-centred, CR3BP)
            syn.(['r' f]) = Sf(:, 1:3);
        else
            syn.(['r' f]) = traj_to_synodic(Sf, ef, mu);
        end
    end

    % Punti chiave
    key_fields = {'halo','inj','dsm1','flyby_pre','flyby_post','dsm2','comet'};
    key_states = {arcs.S_halo_ic, arcs.S_inj_pre, arcs.S_dsm1_pre, ...
                  arcs.S_flyby_pre, arcs.S_flyby_post, arcs.S_dsm2_post, ...
                  arcs.S_comet};
    key_epochs = {arcs.ep_halo, arcs.ep_inj, arcs.ep_dsm1, ...
                  arcs.ep_flyby, arcs.ep_flyby, arcs.ep_dsm2, arcs.ep_comet};

    for k = 1:numel(key_fields)
        S6 = key_states{k}(:)';   % forza 1×6
        r3 = traj_to_synodic(S6, key_epochs{k}, mu);   % 1×3
        syn.(key_fields{k}) = r3;
    end
end


function r_syn = traj_to_synodic(S_J2000, epochs, mu)
% traj_to_synodic  Convert N J2000 states to synodic non-dim positions.
%
%   S_J2000 – N×6  [km; km/s]  Sun-centred ECLIPJ2000
%   epochs  – N×1  [ET s]
%   mu      – CR3BP mass ratio  [-]
%   r_syn   - N×3  barycentre-centred synodic positions [-]

    n = size(S_J2000, 1);

    % Batch SPICE: 6×N
    st_e = cspice_spkezr('EARTH', epochs(:)', 'ECLIPJ2000', 'NONE', 'SUN');
    r_e  = st_e(1:3, :);   % 3×N
    v_e  = st_e(4:6, :);   % 3×N

    L = sqrt(sum(r_e.^2, 1));   % 1×N

    x_hat = r_e ./ L;
    z_raw = cross(r_e, v_e);
    z_hat = z_raw ./ sqrt(sum(z_raw.^2, 1));
    y_hat = cross(z_hat, x_hat);

    r_syn = zeros(n, 3);
    for i = 1:n
        R = [x_hat(:,i), y_hat(:,i), z_hat(:,i)];
        r_sc_syn = R' * S_J2000(i, 1:3)' / L(i);
        r_syn(i,:) = (r_sc_syn - [mu; 0; 0])';
    end
end


function r_syn = j2000_to_synodic_pos(r_J2000, epoch, mu)
% j2000_to_synodic_pos  Convert a single position vector.
    st_e  = cspice_spkezr('EARTH', epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e   = st_e(1:3);  v_e = st_e(4:6);
    L     = norm(r_e);
    x_hat = r_e / L;
    z_hat = cross(r_e, v_e);  z_hat = z_hat / norm(z_hat);
    y_hat = cross(z_hat, x_hat);
    R     = [x_hat, y_hat, z_hat];
    r_sc  = R' * r_J2000(:) / L - [mu; 0; 0];
    r_syn = r_sc(:)';
end


function dv_syn_unit = dv_j2000_to_synodic_unit(dv_J2000, epoch)
% dv_j2000_to_synodic_unit  Rotate a delta-v from J2000 to the synodic frame and normalize it.
%   Useful to visualize the maneuver direction in the synodic plane.
%   dv_J2000 – 3×1 [km/s]
%   epoch    - ET s (impulse epoch)
    st_e  = cspice_spkezr('EARTH', epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e   = st_e(1:3);  v_e = st_e(4:6);
    L     = norm(r_e);
    x_hat = r_e / L;
    z_hat = cross(r_e, v_e);  z_hat = z_hat / norm(z_hat);
    y_hat = cross(z_hat, x_hat);
    R     = [x_hat, y_hat, z_hat];
    dv_syn       = R' * dv_J2000(:);   % J2000 -> dimensional synodic
    dv_syn_unit  = dv_syn / norm(dv_syn);
end


function draw_arcs_3d(syn, ls, lw, label_prefix)
% draw_arcs_3d  Plot the arcs and key points in 3D in the synodic frame.

    % Archi forward
    colors_fwd     = {[0.5 0.7 1.0], [0.0 0.4 0.9], [0.0 0.0 0.7]};
    arc_labels_fwd = {'Halo shift', 'Inj \rightarrow DSM1', 'DSM1 \rightarrow Flyby'};
    for k = 1:3
        r = syn.(['r' num2str(k)]);
        if isempty(r), continue; end
        plot3(r(:,1), r(:,2), r(:,3), ls, 'Color', colors_fwd{k}, ...
              'LineWidth', lw, 'DisplayName', [label_prefix ' ' arc_labels_fwd{k}]);
    end

    % Archi backward (invertiti per mostrare flyby→cometa)
    colors_bwd     = {[1.0 0.4 0.4], [0.7 0.0 0.0]};
    arc_labels_bwd = {'Flyby \rightarrow DSM2', 'DSM2 \rightarrow Comet'};
    for k = [5 4]
        r = flipud(syn.(['r' num2str(k)]));
        if isempty(r), continue; end
        idx = k - 3;   % k=5 -> idx=2 ('Flyby->DSM2'), k=4 -> idx=1 ('DSM2->Comet')
        plot3(r(:,1), r(:,2), r(:,3), ls, 'Color', colors_bwd{idx}, ...
              'LineWidth', lw, 'DisplayName', [label_prefix ' ' arc_labels_bwd{idx}]);
    end

    % Punti chiave
    mk      = {'cs',  'b^',   'bd',      'mo',          'mx',          'rd',        'rp'    };
    knames  = {'halo','inj',  'dsm1',    'flyby_pre',   'flyby_post',  'dsm2',      'comet' };
    klabels = {'Halo','Inj.', 'DSM1',    'Flyby (fwd)', 'Flyby (bwd)', 'DSM2',      'Comet'};
    mfc     = {'c',   'b',    [0 0.6 1], 'm',           'none',        [1 0.4 0],   'r'     };
    for k = 1:numel(knames)
        r = syn.(knames{k});
        plot3(r(1), r(2), r(3), mk{k}, 'MarkerSize', 7, ...
              'MarkerFaceColor', mfc{k}, ...
              'DisplayName', [label_prefix ' ' klabels{k}]);
    end
end


function draw_arcs_xz(syn, ls, lw, label_prefix)
% draw_arcs_xz  Plot the arcs in the x-z plane of the synodic frame.
    colors_fwd = {[0.5 0.7 1.0], [0.0 0.4 0.9], [0.0 0.0 0.7]};
    colors_bwd = {[1.0 0.4 0.4], [0.7 0.0 0.0]};

    arc_labels = {'Attesa','Inj→DSM1','DSM1→Flyby'};
    for k = [1 2 3]
        r = syn.(['r' num2str(k)]);
        if isempty(r), continue; end
        plot(r(:,1), r(:,3), ls, 'Color', colors_fwd{k}, ...
             'LineWidth', lw, 'DisplayName', [label_prefix ' ' arc_labels{k}]);
    end
    arc_labels_bwd = {'Flyby->DSM2','DSM2->Comet'};
    for k = [5 4]
        r = flipud(syn.(['r' num2str(k)]));
        if isempty(r), continue; end
        idx = k - 3;
        plot(r(:,1), r(:,3), ls, 'Color', colors_bwd{idx}, ...
             'LineWidth', lw, 'DisplayName', [label_prefix ' ' arc_labels_bwd{idx}]);
    end
    mk = {'cs','b^','bd','mo','mx','rd','rp'};
    knames = {'halo','inj','dsm1','flyby_pre','flyby_post','dsm2','comet'};
    for k = 1:numel(knames)
        r = syn.(knames{k});
        plot(r(1), r(3), mk{k}, 'MarkerSize', 7, 'HandleVisibility','off');
    end
end


function draw_bodies_3d(mu, moon_syn, L1_x, L2_x)
% draw_bodies_3d  Plot Sun, Earth, Moon, L1, L2 in 3D.
%   Sun, Earth, L1, L2 lie in the z=0 plane.
%   The Moon may have a non-zero z component (real, from ephemerides).
    plot3(-mu,   0, 0, 'yo', 'MarkerSize',12, 'MarkerFaceColor','y',         'DisplayName','Sun');
    plot3(1-mu,  0, 0, 'go', 'MarkerSize',10, 'MarkerFaceColor','g',         'DisplayName','Earth');
    plot3(moon_syn(1), moon_syn(2), moon_syn(3), 'ko', 'MarkerSize',6, ...
          'MarkerFaceColor',[0.5 0.5 0.5],                                   'DisplayName','Moon (flyby)');
    plot3(L1_x,  0, 0, 'k+', 'MarkerSize',8, 'LineWidth',1.5,               'DisplayName','L_1');
    plot3(L2_x,  0, 0, 'kx', 'MarkerSize',8, 'LineWidth',1.5,               'DisplayName','L_2');
end
