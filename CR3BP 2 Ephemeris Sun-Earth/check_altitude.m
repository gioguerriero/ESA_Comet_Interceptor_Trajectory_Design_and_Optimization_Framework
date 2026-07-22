function [h_min_earth, h_flyby, info] = check_altitude(out, c, varargin)
% check_altitude  Minimum Earth altitude and lunar flyby altitude of an
%   ephemeris-refined solution (output of variables_organizer).
%
%   [h_min_earth, h_flyby] = check_altitude(out, c)
%   [h_min_earth, h_flyby, info] = check_altitude(out, c, 'Name',Value,...)
%
%   Inputs:
%     out  - variables_organizer struct (post optimization_ephe). Fields used:
%              .departure.state.J2000, .departure.epoch
%              .injection.J2000_kms
%              .dsm1.state_pre.J2000, .dsm1.dv.J2000_kms, .dsm1.epoch
%              .flyby.state_post.J2000, .flyby.epoch, .flyby.vinf_in, .flyby.vinf_out
%              .dsm2.state_pre.J2000, .dsm2.dv.J2000_kms, .dsm2.epoch
%              .comet.epoch
%     c    - constants struct (rEarth, rMoon, muMoon, mEarth, mSun, G).
%
%   Options (Name-Value):
%     'dt_coarse_days' (0.1)  coarse grid step [days]
%     'win_days'       (2)    refinement half-window around the CA [days]
%     'dt_fine_s'      (60)   fine grid step [s]
%     'verbose'        (true) print the summary
%     'plot'           (true) draw the trajectory in the synodic non-dim frame
%     'S_halo'         ([])   parking halo orbit (Nx6 synodic non-dim):
%                             if provided, drawn dashed in the plot
%
%   Outputs:
%     h_min_earth - minimum altitude above the Earth surface [km]
%     h_flyby     - lunar flyby altitude (zeroSOI) [km]
%     info        - struct with details (CA epoch, distance, segment, etc.)
%                   also includes info.max_sun_dist_km / .max_sun_dist_AU:
%                   maximum heliocentric distance reached along the trajectory
%                   (coarse grid over the 4 segments; SS(:,1:3) is already
%                   Sun-centered ECLIPJ2000, so the distance from the Sun is
%                   simply the norm of the propagated state).
%
%   Method (Earth): the minimum Earth distance is found in two passes.
%     1) COARSE: each segment (dep->DSM1, DSM1->flyby, flyby->DSM2, DSM2->comet)
%        is propagated with NBODY_J2000 on a regular grid and the approximate
%        closest approach is located.
%     2) FINE: from the EXACT segment start state it re-propagates at high
%        accuracy to (CA_epoch - win) and then over a dense grid from
%        (CA_epoch - win) to (CA_epoch + win), finding the true distance.
%     This keeps the minimum free of coarse-integration errors.
%
%   Method (flyby): the lunar flyby altitude is computed from the zeroSOI
%     bending angle between vinf_in and vinf_out (patched flyby, model-consistent).
%
%   Dynamics: NBODY_J2000 (Sun + Earth), the same as the refinement.

% -------------------- parameters --------------------
p = inputParser;
p.addParameter('dt_coarse_days', 0.1);
p.addParameter('win_days',       2);
p.addParameter('dt_fine_s',      60);
p.addParameter('verbose',        true);
p.addParameter('plot',           true);
p.addParameter('S_halo',         []);     % parking halo orbit (N×6 synodic non-dim)
p.addParameter('SaveBlender',    true);   % export blender_input.json
p.addParameter('BlenderFile',    'blender_input.json');
p.addParameter('N_blender',      5000);   % trajectory points in the json
p.addParameter('T_halo_adim',    3.083023392967281);  % halo period [non-dim]
p.addParameter('N_halo_pts',     2000);   % points of the 2 halo orbits
p.parse(varargin{:});
dt_coarse_days = p.Results.dt_coarse_days;
win_days       = p.Results.win_days;
dt_fine_s      = p.Results.dt_fine_s;
verbose        = p.Results.verbose;
do_plot        = p.Results.plot;
S_halo         = p.Results.S_halo;
save_blender   = p.Results.SaveBlender;
blender_file   = p.Results.BlenderFile;
N_blender      = p.Results.N_blender;
T_halo_adim    = p.Results.T_halo_adim;
N_halo_pts     = p.Results.N_halo_pts;

opt_coarse = odeset('AbsTol',1e-11,'RelTol',1e-11);
opt_fine   = odeset('AbsTol',1e-13,'RelTol',1e-13);
AU_km      = 149597870.7;    % 1 AU exact (IAU)

% ==================================================================
%  1) LUNAR FLYBY ALTITUDE (zeroSOI, from bending angle)
% ==================================================================
vin  = out.flyby.vinf_in(:);
vout = out.flyby.vinf_out(:);
vinf_mag = 0.5 * (norm(vin) + norm(vout));          % mean |vinf| (in≈out)
cos_d    = dot(vin, vout) / (norm(vin) * norm(vout));
cos_d    = max(-1, min(1, cos_d));                  % numerical clamp
delta    = acos(cos_d);                             % bending angle [rad]
rp_flyby = c.muMoon / vinf_mag^2 * (1/sin(delta/2) - 1);   % periapsis radius [km]
h_flyby  = rp_flyby - c.rMoon;

% ==================================================================
%  2) BUILD THE 4 SEGMENTS (POST-maneuver state + epochs)
% ==================================================================
% Seg 1: departure (+ injection) -> DSM1
S0_1 = out.departure.state.J2000(:);
S0_1(4:6) = S0_1(4:6) + out.injection.J2000_kms(:);
% Seg 2: DSM1 (+ dv_dsm1) -> flyby
S0_2 = out.dsm1.state_pre.J2000(:);
S0_2(4:6) = S0_2(4:6) + out.dsm1.dv.J2000_kms(:);
% Seg 3: post-flyby -> DSM2  (vinf rotation already folded into state_post)
S0_3 = out.flyby.state_post.J2000(:);
% Seg 4: DSM2 (+ dv_dsm2) -> comet
S0_4 = out.dsm2.state_pre.J2000(:);
S0_4(4:6) = S0_4(4:6) + out.dsm2.dv.J2000_kms(:);

seg(1) = local_seg(S0_1, out.departure.epoch, out.dsm1.epoch,  'dep->DSM1');
seg(2) = local_seg(S0_2, out.dsm1.epoch,      out.flyby.epoch, 'DSM1->flyby');
seg(3) = local_seg(S0_3, out.flyby.epoch,     out.dsm2.epoch,  'flyby->DSM2');
seg(4) = local_seg(S0_4, out.dsm2.epoch,      out.comet.epoch, 'DSM2->comet');

% ==================================================================
%  3) COARSE PASS: locate the closest approach to Earth
% ==================================================================
best.d    = inf;
best.seg  = 0;
best.epca = NaN;
max_sun_dist_km = 0;         % maximum heliocentric distance along the trajectory [km]

for s = 1:numel(seg)
    dt = seg(s).ep1 - seg(s).ep0;
    if abs(dt) < 1            % negligible-duration segment -> skip
        continue;
    end
    npt   = max(2, ceil(abs(dt)/86400 / dt_coarse_days) + 1);
    tgrid = linspace(0, dt, npt).';
    [~, SS] = ode45(@(t,S) NBODY_J2000(t,S, seg(s).ep0, c), tgrid, seg(s).S0, opt_coarse);
    ep  = seg(s).ep0 + tgrid;
    dE  = local_earth_dist(SS(:,1:3), ep);
    [dmin_s, imin] = min(dE);
    if dmin_s < best.d
        best.d    = dmin_s;
        best.seg  = s;
        best.epca = ep(imin);
    end

    % SS(:,1:3) e' gia' Sun-centered ECLIPJ2000 (r_sc in NBODY_J2000.m),
    % so the heliocentric distance is simply its norm.
    d_sun_seg = vecnorm(SS(:,1:3), 2, 2);
    max_sun_dist_km = max(max_sun_dist_km, max(d_sun_seg));
end
max_sun_dist_AU = max_sun_dist_km / AU_km;

if best.seg == 0
    error('check_altitude: no valid segment to propagate.');
end

% ==================================================================
%  4) FINE PASS: refinement around the CA (high accuracy)
% ==================================================================
s   = best.seg;
ep_a = min(seg(s).ep0, seg(s).ep1);    % segment boundaries (sorted)
ep_b = max(seg(s).ep0, seg(s).ep1);

% window +/-win, clamped to the segment boundaries
ep_lo = max(ep_a, best.epca - win_days*86400);
ep_hi = min(ep_b, best.epca + win_days*86400);

% precise state at ep_lo: re-propagate from the EXACT segment start
if abs(ep_lo - seg(s).ep0) < 1
    S_lo = seg(s).S0;
else
    [~, Sp] = ode45(@(t,S) NBODY_J2000(t,S, seg(s).ep0, c), [0, ep_lo - seg(s).ep0], seg(s).S0, opt_fine);
    S_lo = Sp(end,:).';
end

% dense grid over [ep_lo, ep_hi] at high accuracy
npt_f = max(2, ceil((ep_hi - ep_lo) / dt_fine_s) + 1);
tg    = linspace(0, ep_hi - ep_lo, npt_f).';
[~, SSf] = ode45(@(t,S) NBODY_J2000(t,S, ep_lo, c), tg, S_lo, opt_fine);
epf = ep_lo + tg;
dEf = local_earth_dist(SSf(:,1:3), epf);
[d_ca, ica] = min(dEf);

h_min_earth = d_ca - c.rEarth;

% ==================================================================
%  INFO + STAMPA
% ==================================================================
info.dist_ca_km       = d_ca;
info.alt_ca_km        = h_min_earth;
info.epoch_ca         = epf(ica);
info.segment          = seg(best.seg).name;
info.dist_coarse_km   = best.d;
info.flyby_rp_km      = rp_flyby;
info.flyby_delta_deg  = rad2deg(delta);
info.flyby_vinf_kms   = vinf_mag;
info.max_sun_dist_km  = max_sun_dist_km;
info.max_sun_dist_AU  = max_sun_dist_AU;

% ==================================================================
%  BLENDER EXPORT (.json) - trajectory + halo + Moon + maneuvers
%  (non-dimensional Sun-Earth synodic frame)
% ==================================================================
if save_blender
    local_export_blender(out, c, seg, blender_file, ...
                         N_blender, T_halo_adim, N_halo_pts, verbose);
    info.blender_file = blender_file;
end

% ==================================================================
%  5) PLOT - non-dimensional Sun-Earth synodic frame
% ==================================================================
if do_plot
    mu        = c.mu;
    earth_syn = [1-mu, 0, 0];          % Earth fixed in the synodic frame
    opt_plot  = odeset('AbsTol',1e-11,'RelTol',1e-11);

    figure('Color','w','Name','check_altitude — synodic trajectory', ...
           'Position',[100 100 1200 640]);
    hold on; grid on; box on;

    % --- dense trajectory, segment by segment (smooth) ---
    h_traj = [];
    for s = 1:numel(seg)
        dt = seg(s).ep1 - seg(s).ep0;
        if abs(dt) < 1, continue; end
        npt = max(200, round(abs(dt)/86400 * 30));     % ~30 points/day
        tg  = linspace(0, dt, npt).';
        [~, SS] = ode45(@(t,S) NBODY_J2000(t,S, seg(s).ep0, c), tg, seg(s).S0, opt_plot);
        r_syn = local_to_synodic_pos(SS(:,1:3), seg(s).ep0 + tg, mu);
        hh = plot3(r_syn(:,1), r_syn(:,2), r_syn(:,3), '-', ...
                   'Color',[0 0.45 0.74], 'LineWidth',1.3);
        if isempty(h_traj)
            h_traj = hh;  set(hh,'DisplayName','Trajectory');
        else
            set(hh,'HandleVisibility','off');
        end
    end

    % --- parking Halo orbit (dashed) ---
    % S_halo is already in synodic non-dim CR3BP (same frame as the plot).
    if ~isempty(S_halo)
        rH = [S_halo(:,1:3); S_halo(1,1:3)];        % close the periodic loop
        plot3(rH(:,1), rH(:,2), rH(:,3), '--', ...
              'Color',[0.20 0.60 0.20], 'LineWidth',1.1, ...
              'DisplayName','Halo parking orbit');
    end

    % --- 1 lunar orbit before the flyby (dashed) ---
    % In the rotating Sun-Earth frame the Moon returns to itself after one
    % SYNODIC month (~29.53 d), not sidereal: that closes the loop.
    T_moon_s = 29.530589 * 86400;                   % mese sinodico [s]
    ep_moon  = linspace(out.flyby.epoch - T_moon_s, out.flyby.epoch, 500).';
    MS       = cspice_spkezr('MOON', ep_moon.', 'ECLIPJ2000', 'NONE', 'SUN');  % 6×N
    rMoon_syn = local_to_synodic_pos(MS(1:3,:).', ep_moon, mu);
    plot3(rMoon_syn(:,1), rMoon_syn(:,2), rMoon_syn(:,3), '--', ...
          'Color',[0.50 0.50 0.50], 'LineWidth',1.1, ...
          'DisplayName','Moon orbit (1 rev pre-flyby)');

    % --- Earth at real scale (sphere of radius rEarth/Lstar) ---
    [xs,ys,zs] = sphere(48);
    rE_DU = c.rEarth / c.Lstar;
    surf(xs*rE_DU + earth_syn(1), ys*rE_DU + earth_syn(2), zs*rE_DU + earth_syn(3), ...
         'EdgeColor','none', 'FaceColor',[0.10 0.45 0.90], 'FaceAlpha',1.0, ...
         'DisplayName','Earth (true scale)');

    % --- Sole ---
    plot3(-mu, 0, 0, 'o', 'MarkerSize',12, 'MarkerFaceColor',[1 0.8 0], ...
          'MarkerEdgeColor','k', 'DisplayName','Sun');

    % --- eventi (departure, DSM1, flyby, DSM2, comet) ---
    ev_r = [out.departure.state.J2000(1:3);
            out.dsm1.state_pre.J2000(1:3);
            out.flyby.state_pre.J2000(1:3);
            out.dsm2.state_pre.J2000(1:3);
            out.comet.state.J2000(1:3)];
    ev_ep   = [out.departure.epoch; out.dsm1.epoch; out.flyby.epoch; out.dsm2.epoch; out.comet.epoch];
    ev_syn  = local_to_synodic_pos(ev_r, ev_ep, mu);
    ev_name = {'Departure','DSM1','Flyby','DSM2','Comet'};
    plot3(ev_syn(:,1), ev_syn(:,2), ev_syn(:,3), 'ks', ...
          'MarkerFaceColor',[0.95 0.85 0.1], 'MarkerSize',6, 'DisplayName','Events');
    text(ev_syn(:,1), ev_syn(:,2), ev_syn(:,3), ev_name, ...
         'FontSize',8, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');

    % --- closest approach: point + Earth-spacecraft segment ---
    r_ca = local_to_synodic_pos(SSf(ica,1:3), epf(ica), mu);
    plot3([earth_syn(1) r_ca(1)], [earth_syn(2) r_ca(2)], [earth_syn(3) r_ca(3)], ...
          'r-', 'LineWidth',1.6, 'DisplayName','Earth-SC segment (CA)');
    plot3(r_ca(1), r_ca(2), r_ca(3), 'o', 'MarkerSize',7, ...
          'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'DisplayName','Closest approach');

    info.r_ca_syn = r_ca;

    axis equal;
    set(gca,'Position',[0.055 0.09 0.58 0.84]);   % leave room on the right for the panel
    xlabel('x [DU]'); ylabel('y [DU]'); zlabel('z [DU]');
    title('Trajectory — Sun-Earth synodic frame (nondim)');
    legend('Location','best');
    view(2);
    hold off;

    % ==============================================================
    %  SIDE STATS PANEL — card pulita (font proporzionale, due colonne)
    % ==============================================================
    d2s    = 86400;
    tof_dd = (out.dsm1.epoch  - out.departure.epoch) / d2s;
    tof_df = (out.flyby.epoch - out.dsm1.epoch)      / d2s;
    tof_fd = (out.dsm2.epoch  - out.flyby.epoch)     / d2s;
    tof_dc = (out.comet.epoch - out.dsm2.epoch)      / d2s;
    tof_tot= (out.comet.epoch - out.departure.epoch) / d2s;

    dv_i = out.injection.norm_ms;
    dv_1 = out.dsm1.dv.norm_ms;
    dv_2 = out.dsm2.dv.norm_ms;
    dv_t = dv_i + dv_1 + dv_2;

    s_dep = cspice_et2utc(out.departure.epoch, 'C', 0);   dep_str = s_dep(1:end-3);  % togli i secondi
    s_arr = cspice_et2utc(out.comet.epoch,     'C', 0);   arr_str = s_arr(1:end-3);

    % dedicated panel axis (0..1 coordinates), invisible
    axp = axes('Position',[0.665 0.05 0.315 0.90]);
    hold(axp,'on');  set(axp,'XLim',[0 1],'YLim',[0 1]);  axis(axp,'off');

    % card di sfondo: bianca, bordo morbido, angoli arrotondati
    rectangle(axp,'Position',[0.005 0.005 0.99 0.99],'Curvature',0.04, ...
        'FaceColor','w','EdgeColor',[0.86 0.86 0.86],'LineWidth',1.2);

    % titolo
    y = 0.955;
    text(axp,0.5,y,'Mission Summary','FontSize',13,'FontWeight','bold', ...
        'Color',[0.15 0.15 0.15],'HorizontalAlignment','center','VerticalAlignment','middle');
    y = y - 0.032;
    line(axp,[0.06 0.94],[y y],'Color',[0.80 0.80 0.80],'LineWidth',1.0);
    y = y - 0.040;

    % --- Time of flight ---
    y = local_panel_header(axp,y,'Time of Flight  [days]');
    y = local_panel_row(axp,y,'Departure \rightarrow DSM1', sprintf('%.1f',tof_dd));
    y = local_panel_row(axp,y,'DSM1 \rightarrow Flyby',     sprintf('%.1f',tof_df));
    y = local_panel_row(axp,y,'Flyby \rightarrow DSM2',     sprintf('%.1f',tof_fd));
    y = local_panel_row(axp,y,'DSM2 \rightarrow Comet',     sprintf('%.1f',tof_dc));
    y = local_panel_row(axp,y,'Total',                      sprintf('%.1f',tof_tot), true);

    % --- Delta-v budget ---
    y = y - 0.014;
    y = local_panel_header(axp,y,'\Deltav Budget  [m/s]');
    y = local_panel_row(axp,y,'Injection', sprintf('%.1f',dv_i));
    y = local_panel_row(axp,y,'DSM1',      sprintf('%.1f',dv_1));
    y = local_panel_row(axp,y,'DSM2',      sprintf('%.1f',dv_2));
    y = local_panel_row(axp,y,'Total',     sprintf('%.1f',dv_t), true);

    % --- Epochs ---
    y = y - 0.014;
    y = local_panel_header(axp,y,'Epochs  (UTC)');
    y = local_panel_row(axp,y,'Departure', dep_str);
    y = local_panel_row(axp,y,'Arrival',   arr_str);

    % --- Flyby & closest approach ---
    y = y - 0.014;
    y = local_panel_header(axp,y,'Flyby & Closest Approach');
    y = local_panel_row(axp,y,'Lunar flyby alt.', sprintf('%.0f km',  h_flyby));
    y = local_panel_row(axp,y,'v_\infty (flyby)', sprintf('%.3f km/s', vinf_mag));
    y = local_panel_row(axp,y,'Earth CA alt.',    sprintf('%.0f km',  h_min_earth));
end

if verbose
    fprintf('\n--- CHECK ALTITUDE ---\n');
    fprintf('  Min Earth distance : %.2f km   (altitudine %.2f km)\n', d_ca, h_min_earth);
    fprintf('    @ %s  | segment: %s\n', cspice_et2utc(epf(ica),'C',3), seg(best.seg).name);
    fprintf('    coarse -> fine    : %.2f km -> %.2f km  (correzione %.2f km)\n', ...
        best.d, d_ca, best.d - d_ca);
    fprintf('  Lunar flyby alt    : %.2f km   (rp %.2f km, delta %.2f deg, vinf %.4f km/s)\n', ...
        h_flyby, rp_flyby, rad2deg(delta), vinf_mag);
    fprintf('  Max heliocentric distance : %.4f AU  (%.0f km)\n', ...
        max_sun_dist_AU, max_sun_dist_km);
    fprintf('-----------------------\n\n');
end

end


% ======================================================================
%  HELPER: build a segment struct
% ======================================================================
function s = local_seg(S0, ep0, ep1, name)
    s.S0   = S0(:);
    s.ep0  = ep0;
    s.ep1  = ep1;
    s.name = name;
end


% ======================================================================
%  HELPER: Earth distance of Sun-centered ECLIPJ2000 states
%    r_sc   : N×3 [km]   (heliocentric positions)
%    epochs : N×1 [ET s]
%    d      : N×1 [km]
% ======================================================================
function d = local_earth_dist(r_sc, epochs)
    ES = cspice_spkezr('EARTH', epochs(:).', 'ECLIPJ2000', 'NONE', 'SUN');  % 6×N
    rE = ES(1:3,:).';                       % N×3
    d  = vecnorm(r_sc - rE, 2, 2);          % N×1
end


% ======================================================================
%  HELPER: Sun-centered ECLIPJ2000 position -> non-dimensional synodic
%    (position only; "pulsating" frame normalized by L(epoch), consistent
%     with synodic2sun_J2000 used elsewhere in the code)
%    r_J    : N×3 [km]   epochs : N×1 [ET s]
%    r_syn  : N×3 [DU]   (barycentre at the origin, Earth at (1-mu,0,0))
% ======================================================================
function r_syn = local_to_synodic_pos(r_J, epochs, mu)
    epochs = epochs(:).';
    N  = numel(epochs);
    ES = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');   % 6×N
    rE = ES(1:3,:);   vE = ES(4:6,:);
    L  = sqrt(sum(rE.^2, 1));
    h  = cross(rE, vE);   hn = sqrt(sum(h.^2, 1));
    e1 = rE ./ L;
    e3 = h  ./ hn;
    e2 = cross(e3, e1);

    r_J   = reshape(r_J, N, 3);
    r_syn = zeros(N, 3);
    for i = 1:N
        R = [e1(:,i), e2(:,i), e3(:,i)];
        r_syn(i,:) = (R' * r_J(i,:).' / L(i) - [mu;0;0]).';
    end
end


% ======================================================================
%  PANEL HELPER: section heading (with a separator line)
% ======================================================================
function y2 = local_panel_header(ax, y, txt)
    accent = [0.00 0.42 0.72];
    text(ax, 0.06, y, txt, 'FontSize',10.5, 'FontWeight','bold', 'Color',accent, ...
         'Interpreter','tex', 'HorizontalAlignment','left', 'VerticalAlignment','middle');
    line(ax, [0.06 0.94], [y-0.018 y-0.018], 'Color',[0.88 0.88 0.88], 'LineWidth',0.8);
    y2 = y - 0.046;
end


% ======================================================================
%  PANEL HELPER: "label .......... value" row (two columns)
% ======================================================================
function y2 = local_panel_row(ax, y, label, value, bold)
    if nargin < 5 || isempty(bold), bold = false; end
    if bold, fw = 'bold';   lc = [0.10 0.10 0.10];
    else,    fw = 'normal'; lc = [0.32 0.32 0.32];
    end
    text(ax, 0.09, y, label, 'FontSize',10, 'FontWeight',fw, 'Color',lc, ...
         'Interpreter','tex', 'HorizontalAlignment','left',  'VerticalAlignment','middle');
    text(ax, 0.91, y, value, 'FontSize',10, 'FontWeight',fw, 'Color',[0.10 0.10 0.10], ...
         'Interpreter','tex', 'HorizontalAlignment','right', 'VerticalAlignment','middle');
    y2 = y - 0.040;
end


% ======================================================================
%  HELPER: JSON export for the dynamic Blender plot
%
%  Salva (tutto in frame SINODICO Sun-Earth ADIMENSIONALE):
%    - trajectory : N time-uniform points (x,y,z) + t [s from the 1st point]
%    - halo       : 2 orbite halo prima dell'injection (back-prop CR3BP della
%                   departure.state.syn -> end at the injection point)
%    - moon       : Moon position at the epochs (halo + mission)
%    - maneuvers  : 3 components (synodic non-dim) + magnitude in m/s + epoch + pos
%
%  Time base: t = 0 at the first halo point (= epoch_injection - 2*halo_period).
% ======================================================================
function local_export_blender(out, c, seg, fname, N, T_halo_adim, N_halo, verbose)

    mu    = c.mu;
    Tstar = c.Tstar;
    optp  = odeset('AbsTol',1e-11,'RelTol',1e-11);
    optc  = odeset('AbsTol',1e-12,'RelTol',1e-12);

    ep_dep   = seg(1).ep0;            % departure = injection epoch
    ep_comet = seg(end).ep1;

    % --- time base: t0 = start of the 2 halo orbits ---
    period2_s = 2 * T_halo_adim * Tstar;     % durata 2 orbite [s]
    t0_ET     = ep_dep - period2_s;

    %% ---- MISSION TRAJECTORY (dense -> resampled to N uniform) ----
    EP = [];  R = zeros(0,3);  first = true;
    for s = 1:numel(seg)
        dt = seg(s).ep1 - seg(s).ep0;
        if abs(dt) < 1, continue; end
        npt = max(200, round(abs(dt)/86400 * 30));     % ~30 points/day
        tg  = linspace(0, dt, npt).';
        [~, SS] = ode45(@(t,S) NBODY_J2000(t,S, seg(s).ep0, c), tg, seg(s).S0, optp);
        ep  = seg(s).ep0 + tg;
        rs  = local_to_synodic_pos(SS(:,1:3), ep, mu);
        if ~first
            ep(1) = [];  rs(1,:) = [];   % avoid a duplicated epoch at the boundary
        end
        EP = [EP; ep];  R = [R; rs];  first = false;   %#ok<AGROW>
    end

    [EPu, ia] = unique(EP);  Ru = R(ia,:);             % strettamente crescente
    ep_u = linspace(ep_dep, ep_comet, N).';            % uniform grid
    xq = interp1(EPu, Ru(:,1), ep_u, 'pchip');
    yq = interp1(EPu, Ru(:,2), ep_u, 'pchip');
    zq = interp1(EPu, Ru(:,3), ep_u, 'pchip');
    t_traj = ep_u - t0_ET;

    %% ---- 2 HALO ORBITS (CR3BP back-prop of departure.state.syn) ----
    S_dep_syn = out.departure.state.syn(:);
    tau = linspace(0, -2*T_halo_adim, N_halo).';       % tempo adim (backward)
    [~, SH] = ode45(@(t,S) CR3BP(t,S,mu), tau, S_dep_syn, optc);
    ep_h = ep_dep + tau * Tstar;                        % epoche reali
    ep_h = flipud(ep_h);  SH = flipud(SH(:,1:3));       % ordine cronologico
    t_halo = ep_h - t0_ET;

    %% ---- MOON over the whole timeline (halo + mission) ----
    ep_all = [ep_h; ep_u];
    MS     = cspice_spkezr('MOON', ep_all.', 'ECLIPJ2000', 'NONE', 'SUN');  % 6×M
    rMoon  = local_to_synodic_pos(MS(1:3,:).', ep_all, mu);
    t_moon = ep_all - t0_ET;

    %% ---- MANOVRE ----
    names = {'injection', 'dsm1', 'dsm2'};
    eps_m = [out.departure.epoch; out.dsm1.epoch; out.dsm2.epoch];
    pos_m = {out.departure.state.syn(1:3), ...
             out.dsm1.state_pre.syn(1:3), ...
             out.dsm2.state_pre.syn(1:3)};
    dv_m  = {out.injection.syn_adim, ...
             out.dsm1.dv.syn_adim, ...
             out.dsm2.dv.syn_adim};
    dvms  = [out.injection.norm_ms; out.dsm1.dv.norm_ms; out.dsm2.dv.norm_ms];

    man = struct('name',{},'epoch_et',{},'t',{},'epoch_utc',{}, ...
                 'pos_syn',{},'dv_syn_adim',{},'dv_ms',{});
    for k = 1:3
        man(k).name        = names{k};
        man(k).epoch_et    = eps_m(k);
        man(k).t           = eps_m(k) - t0_ET;
        man(k).epoch_utc   = cspice_et2utc(eps_m(k), 'ISOC', 3);
        man(k).pos_syn     = pos_m{k}(:).';
        man(k).dv_syn_adim = dv_m{k}(:).';
        man(k).dv_ms       = dvms(k);
    end

    %% ---- ASSEMBLA STRUCT JSON ----
    J.frame   = 'Sun-Earth synodic, nondimensional (DU). Barycenter at origin, Earth at (1-mu,0,0), Sun at (-mu,0,0).';
    J.mu      = mu;
    J.Lstar_km = c.Lstar;
    J.Tstar_s  = Tstar;
    J.t0_et    = t0_ET;
    J.t0_utc   = cspice_et2utc(t0_ET, 'ISOC', 3);
    J.units    = struct('position','DU (nondim)', 'time','s from first point', ...
                        'dv_components','synodic nondim', 'dv_amplitude','m/s');

    J.trajectory = struct('t',t_traj(:).', 'x',xq(:).', 'y',yq(:).', 'z',zq(:).');
    J.halo       = struct('t',t_halo(:).', 'x',SH(:,1).', 'y',SH(:,2).', 'z',SH(:,3).');
    J.moon       = struct('t',t_moon(:).', 'x',rMoon(:,1).', 'y',rMoon(:,2).', 'z',rMoon(:,3).');
    J.maneuvers  = man;

    %% ---- SCRITTURA FILE ----
    try
        txt = jsonencode(J, 'PrettyPrint', true);   % R2021a+
    catch
        txt = jsonencode(J);
    end
    fid = fopen(fname, 'w');
    if fid < 0
        warning('local_export_blender: impossibile aprire "%s"', fname);
        return;
    end
    fwrite(fid, txt, 'char');
    fclose(fid);

    if verbose
        fprintf('  [blender] scritto %s : traj %d pts, halo %d pts, moon %d pts\n', ...
                fname, numel(t_traj), numel(t_halo), numel(t_moon));
    end
end
