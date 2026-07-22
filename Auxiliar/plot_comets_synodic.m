function pos_syn_ad = plot_comets_synodic(comets, mu, varargin)
% PLOT_COMETS_SYNODIC  Plot comet encounter positions in adimensional synodic frame.
%
%   plot_comets_synodic(comets, mu)
%   plot_comets_synodic(comets, mu, 'PlotManifolds', true, 'S_halo', S_halo, ...
%                        'UnstableDir', unstable_dir, 'C', c)
%
%   Inputs:
%       comets  - cell array of structs, each with fields:
%                   .comet_pos  [1x3] km, heliocentric ECLIPJ2000
%                   .epoch      scalar, seconds past J2000 (ET)
%                   .name       string, used in legend
%       mu      - Sun-Earth mass parameter (optional, default: mEarth/(mSun+mEarth))
%
%   The synodic frame is the Sun-Earth CR3BP rotating frame at each
%   comet's encounter epoch:
%       - Origin: barycenter  (≈ Sun for Sun-Earth)
%       - x-axis: Sun → Earth direction at encounter epoch
%       - Lengths normalised by Earth-Sun distance at that epoch
%       - Sun at (-mu, 0), Earth at (1-mu, 0)
%
%   Optional manifold overlay (Name,Value), all in the SAME synodic adim
%   frame as the comets/Sun/Earth (no extra rotation needed):
%       'PlotManifolds' (false)  - true to overlay unstable manifold arcs
%                                  departing the parking halo orbit
%       'S_halo'        ([])     - Nx6 halo orbit states (sinodico adim),
%                                  REQUIRED if PlotManifolds = true
%       'UnstableDir'   ([])     - Nx3 unstable eigenvector directions
%                                  (one per S_halo row), REQUIRED if
%                                  PlotManifolds = true
%       'C'             ([])     - constants struct (needs .mu, .Tstar,
%                                  .Vstar), REQUIRED if PlotManifolds = true
%       'EpsVelMs'      (15)     - injection velocity perturbation [m/s]
%                                  applied along the unstable eigenvector
%       'TmaxYears'     (6)      - forward propagation time [years]
%       'NManifolds'    (200)    - number of departure points sampled
%                                  (equally spaced by index) along S_halo
%       'NptsManifold'  (400)    - samples per propagated arc (plot res.)
%
%   Example (comets only):
%       comets = {C2023X1, C2013US10, C2001Q4, C2008A1};
%       plot_comets_synodic(comets, c.mu);
%
%   Example (comets + unstable manifolds):
%       load('S_halo.mat'); load('unstable_dir.mat');
%       plot_comets_synodic(comets, c.mu, 'PlotManifolds', true, ...
%           'S_halo', S_halo, 'UnstableDir', unstable_dir, 'C', c);

%% ---- optional manifold-overlay parameters -----------------------------
p = inputParser;
p.addParameter('PlotManifolds', false);
p.addParameter('S_halo',        []);
p.addParameter('UnstableDir',   []);
p.addParameter('C',             []);
p.addParameter('EpsVelMs',      15);
p.addParameter('TmaxYears',     6);
p.addParameter('NManifolds',    200);
p.addParameter('NptsManifold',  400);
p.parse(varargin{:});
o = p.Results;

if o.PlotManifolds
    if isempty(o.S_halo) || isempty(o.UnstableDir) || isempty(o.C)
        error(['plot_comets_synodic: ''PlotManifolds'',true richiede anche ' ...
               '''S_halo'', ''UnstableDir'' e ''C''.']);
    end
end

if nargin < 2 || isempty(mu)
    mSun   = 1.9885e30;    % kg
    mEarth = 5.9722e24;    % kg
    mu = mEarth / (mSun + mEarth);
end

markers = {'o','s','d','^','p','h','v','>','<','*'};
colors  = [0.00 0.45 0.74;
           0.85 0.33 0.10;
           0.47 0.67 0.19;
           0.49 0.18 0.56;
           0.93 0.69 0.13;
           0.30 0.75 0.93;
           0.64 0.08 0.18;
           0.07 0.62 0.25;
           0.80 0.40 0.00;
           0.40 0.40 0.40];

n = numel(comets);

pos_syn_ad = zeros(n, 3);   % output: [n x 3], adimensional synodic positions

figure; hold on; grid on; axis equal; box on;
xlabel('X [DU]','FontSize',12);
ylabel('Y [DU]','FontSize',12);
zlabel('Z [DU]','FontSize',12);
set(gca,'FontSize',11);

%% ---- unstable manifold overlay (drawn first, underneath everything) ---
if o.PlotManifolds
    Ch = o.C;
    eps_vel_ad = o.EpsVelMs / (1e3 * Ch.Vstar);
    t_end      = o.TmaxYears * 365 * 24 * 3600 / Ch.Tstar;
    tspan      = linspace(0, t_end, o.NptsManifold);
    opt_man    = odeset('AbsTol', 1e-9, 'RelTol', 1e-9);

    idx_dep = round(linspace(1, size(o.S_halo,1), o.NManifolds));

    h_manifold = gobjects(0);
    for kk = 1:numel(idx_dep)
        ih = idx_dep(kk);
        x0 = o.S_halo(ih,:)';
        vu = o.UnstableDir(ih,:)';
        vu = vu / norm(vu);
        x0_pert = [x0(1:3); x0(4:6) + eps_vel_ad*vu];

        try
            [~, S_man] = ode45(@(t,S) CR3BP(t,S,Ch.mu), tspan, x0_pert, opt_man);
        catch
            continue   % skip arcs that fail to integrate (e.g. Earth close approach)
        end

        hh = plot3(S_man(:,1), S_man(:,2), S_man(:,3), '-', ...
            'Color', [0.20 0.55 0.85], 'LineWidth', 0.6, ...
            'HandleVisibility', 'off');
        hh.Color(4) = 0.20;   % alpha: 80% transparency
        h_manifold(end+1) = hh; %#ok<AGROW>
    end

    if ~isempty(h_manifold)
        set(h_manifold(1), 'HandleVisibility', 'on', 'DisplayName', 'Unstable manifold');
    end
end

h_comets = gobjects(n,1);

for k = 1:n
    comet = comets{k};

    earth_state = cspice_spkezr('EARTH', comet.epoch, 'ECLIPJ2000', 'NONE', 'SUN');
    earth_pos   = earth_state(1:3)';

    theta = atan2(earth_pos(2), earth_pos(1));
    Rz    = [ cos(theta)  sin(theta)  0;
             -sin(theta)  cos(theta)  0;
              0           0           1];

    Lstar_k        = norm(earth_pos(1:2));           % actual AU at this epoch
    pos_syn        = Rz * comet.comet_pos(1:3)';     % rotate to synodic x-axis
    pos_syn_ad_k   = pos_syn / Lstar_k;              % normalise
    pos_syn_ad(k,:) = pos_syn_ad_k';                 % store in output matrix

    mk  = markers{mod(k-1, numel(markers)) + 1};
    col = colors(mod(k-1, size(colors,1)) + 1, :);

    h_comets(k) = plot3(pos_syn_ad_k(1), pos_syn_ad_k(2), pos_syn_ad_k(3), mk, ...
        'Color', col, ...
        'MarkerFaceColor', col, ...
        'MarkerSize', 9, ...
        'LineWidth', 1.5, ...
        'DisplayName', comet.name);
end

% --- Celestial bodies ---
h_sun   = plot(-mu, 0, 'o', 'MarkerSize', 14, ...
    'MarkerFaceColor', [1 0.8 0], 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Sun');
h_earth = plot(1-mu, 0, 'o', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
    'DisplayName', 'Earth');

% --- Earth orbit (unit circle in synodic frame) ---
th_vec      = linspace(0, 2*pi, 300);
h_earth_orb = plot(cos(th_vec), sin(th_vec), ...
    '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, ...
    'DisplayName', 'Earth orbit');

if o.PlotManifolds && ~isempty(h_manifold)
    legend([h_manifold(1); h_comets; h_sun; h_earth; h_earth_orb], ...
        'Location', 'best', 'FontSize', 10);
else
    legend([h_comets; h_sun; h_earth; h_earth_orb], ...
        'Location', 'best', 'FontSize', 10);
end

view(2);
end
