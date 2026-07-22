function res = earth_altitude_check(global_results, c, S_halo, unstable_dir, eps_vel_ms, varargin)
% earth_altitude_check  Check the minimum Earth flyby altitude.
%
%   For each global-search solution it reconstructs the full CR3BP trajectory
%   (Halo->DSM1->lunar flyby->DSM2->comet, same arcs as vinf_escape_study.m /
%   plot_paper_trajectory.m), computes the minimum distance from the Earth
%   centre along the whole path and converts it to a minimum ALTITUDE above
%   the Earth surface [km]:
%
%       d_geo_min = min_trajectory |r_syn - [1-mu,0,0]|      (non-dim)
%       h_min     = d_geo_min * Lstar - rEarth               (km)
%
%   Used to verify that no trajectory passes too close to (or below) the
%   Earth. Plots the DISTRIBUTION of minimum altitudes [km] with the
%   reference lines at 10 km, 50 km, 100 km and GEO (~35786 km).
%
%   NB: the lunar-flyby / Moon->comet arcs always stay >~380000 km from Earth
%   (the Moon is at 384400 km), so the minimum altitude almost always falls on
%   the Halo->flyby leg; the other arcs are still included in the minimum for
%   robustness.
%
% -------------------------------------------------------------------------
%  INPUT
% -------------------------------------------------------------------------
%   global_results - cell array (or struct array) of global-search solutions.
%                    Each solution sol must have:
%                       sol.x_opt              [fpa oop halo_th tof_moon tof_halo]
%                       sol.result_fmincon     [Vinf th fpa oop tof beta dv(1:3)]
%                       sol.fmincon_vinf_results.theta_moon / .vinf_fmincon
%   c              - constants struct (c.mu, c.Vstar, c.Lstar, c.rEarth, ...)
%   S_halo         - Nx6 halo states (synodic non-dim)
%   unstable_dir   - Nx3 unstable directions (velocity) on the halo
%   eps_vel_ms     - injection magnitude from the manifold [m/s] (search_params.eps_vel_ms)
%
%   The closest approach is then REFINED: from the coarse-grid minimum point it
%   re-propagates forward and backward up to +/-RefineDays days with Nrefine
%   points (high resolution) to pinpoint the true minimum. The window is CLAMPED
%   to the arc containing the minimum, so the propagation never crosses a DSM /
%   flyby node (discontinuous velocity): the refinement stays within the same
%   dynamics.
%
%   Options ('Name',value):
%     'NArc'       (600)   - ode45 points for each of the 4 arcs
%     'RefineDays' (5)     - +/- window for the refinement [days]
%     'Nrefine'    (2000)  - points per branch in the refinement
%     'AbsTol'     (1e-10) - ode45 absolute tolerance
%     'RelTol'     (1e-10) - ode45 relative tolerance
%     'Plot'       (true)  - draw the altitude distribution
%     'Verbose'    (true)  - print the text summary
%
% -------------------------------------------------------------------------
%  OUTPUT (struct res)
% -------------------------------------------------------------------------
%   res.min_alt_km    - minimum altitude above Earth per solution [km]
%   res.idx_valid     - indices (in global_results) of the computed solutions
%   res.n_below_1000  - number of solutions with h_min < 1000 km
%   res.n_below_10k   - number of solutions with h_min < 10 000 km
%   res.n_below_geo   - number of solutions with h_min < GEO (~35786 km)
%   res.n_below_100k  - number of solutions with h_min < 100 000 km
%   res.n_collision   - number of solutions with h_min < 0 km (Earth impact)
%   res.geo_alt_km    - GEO altitude used as reference [km]
%   res.max_sun_dist_AU         - maximum heliocentric distance per solution [AU]
%   res.max_sun_dist_overall_AU - maximum heliocentric distance over the whole dataset [AU]
%   res.max_sun_dist_idx        - index (in global_results) of the corresponding solution
%
%   EXAMPLE (from the command window):
%     addpath('Auxiliar/'); addpath('Moon2Comet/');
%     S = load('Synthetic1_global_search.mat');   % -> global_results, c, search_params
%     load('S_halo.mat'); load('unstable_dir.mat');
%     res = earth_altitude_check(S.global_results, S.c, S_halo, unstable_dir, ...
%                                S.search_params.eps_vel_ms);

%% ---- opzioni ----------------------------------------------------------
p = inputParser;
p.addParameter('NArc',    600);
p.addParameter('AbsTol',  1e-10);
p.addParameter('RelTol',  1e-10);
p.addParameter('RefineDays', 5);      % +/- window for the closest-approach refinement [days]
p.addParameter('Nrefine',    2000);   % points per branch in the refinement (high resolution)
p.addParameter('Plot',    true);
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

%% ---- costanti ---------------------------------------------------------
mu     = c.mu;
Vstar  = c.Vstar;
Lstar  = c.Lstar;
rEarth = c.rEarth;
earth  = [1-mu; 0; 0];                  % posizione Terra (sinodico)
sun    = [-mu;  0; 0];                  % posizione Sole  (sinodico)
AU_km  = 149597870.7;                   % 1 AU esatta (IAU), NON c.Lstar (~1.496e8 km, approssimato)
eps_vel_ad = eps_vel_ms / (1e3 * Vstar);
opt    = odeset('AbsTol', o.AbsTol, 'RelTol', o.RelTol);
Narc   = o.NArc;
Nrefine      = o.Nrefine;
delta_t_days = o.RefineDays;
delta_t_nd   = delta_t_days * 24 * 3600 / c.Tstar;   % refinement window [CR3BP non-dim]
geo_alt_km = 42164 - rEarth;            % GEO altitude (~35786 km)

if isstruct(global_results), global_results = num2cell(global_results); end
Nsol = numel(global_results);

min_alt_km      = nan(Nsol,1);
max_sun_dist_AU = nan(Nsol,1);

%% ---- loop over the solutions ---------------------------------------------
for k = 1:Nsol
    try
        sol = global_results{k};
        xo  = sol.x_opt;                 % [fpa oop halo_th tof_moon tof_halo]
        rf  = sol.result_fmincon;        % [Vinf th fpa oop tof beta dv(1:3)]
        fvr = sol.fmincon_vinf_results;

        ref_fpa = xo(1);  ref_oop = xo(2);  ref_halo_th = xo(3);
        ref_tof_moon = xo(4);  ref_tof_halo = xo(5);
        Vinf_m2c = rf(1); theta_m2c = rf(2); fpa_m2c = rf(3);
        oop_m2c  = rf(4); tof_m2c   = rf(5); beta_m2c = rf(6);
        dv_dsm2  = rf(7:9); dv_dsm2 = dv_dsm2(:);

        % --- Arc 1: Halo departure (with injection) -> DSM1 ---
        [init_st, idx_h] = state_finder(ref_halo_th, S_halo);
        init_st = init_st(:);
        vu = unstable_dir(idx_h,:)';  vu = vu / norm(vu);
        S0_arc1 = [init_st(1:3); init_st(4:6) + eps_vel_ad*vu];
        [~, S1] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, ref_tof_halo, Narc), S0_arc1, opt);

        % --- Arc 2: DSM1 -> lunar flyby (backward from the Moon, then flipped) ---
        moon_h2m = moon_state(fvr.theta_moon, c)';
        vinf_h2m = vinf_rotation(moon_h2m, fvr.vinf_fmincon, ref_fpa, ref_oop);
        S0_flyby = moon_h2m;  S0_flyby(4:6) = S0_flyby(4:6) + vinf_h2m ./ Vstar;
        [~, S2b] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, -ref_tof_moon, Narc), S0_flyby, opt);
        S2 = flipud(S2b);

        % --- Arc 3: lunar flyby -> DSM2 ---
        moon_m2c = moon_state(theta_m2c, c);
        vinf_m2c = vinf_rotation(moon_m2c', Vinf_m2c, fpa_m2c, oop_m2c);
        S0_arc3  = [moon_m2c(1:3)'; moon_m2c(4:6)' + vinf_m2c ./ Vstar];
        [~, S3]  = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, tof_m2c*beta_m2c, Narc), S0_arc3, opt);

        % --- Arc 4: DSM2 -> comet ---
        S0_arc4 = [S3(end,1:3)'; S3(end,4:6)' + dv_dsm2];
        [~, S4] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, tof_m2c*(1-beta_m2c), Narc), S0_arc4, opt);

        % --- minimum Earth distance (coarse grid), ARC BY ARC ---
        %   Keep the arcs separate to NOT lose the identity of the arc the minimum
        %   belongs to: this is needed for the clamped refinement (below).
        arcs  = {S1, S2, S3, S4};
        T_arc = [ref_tof_halo, ref_tof_moon, tof_m2c*beta_m2c, tof_m2c*(1-beta_m2c)];  % durate [adim]

        d_min_ad = Inf;  a_best = 0;  j_best = 0;     % coarse minimum, arc and local index
        d_sun_max_ad = 0;                             % maximum heliocentric distance [non-dim]
        for a = 1:numel(arcs)
            d_a = vecnorm(arcs{a}(:,1:3) - earth.', 2, 2);   % geocentric distance on the arc [non-dim]
            [d_a_min, j_a] = min(d_a);
            if d_a_min < d_min_ad
                d_min_ad = d_a_min;  a_best = a;  j_best = j_a;
            end

            d_helio = vecnorm(arcs{a}(:,1:3) - sun.', 2, 2); % heliocentric distance on the arc [non-dim]
            d_sun_max_ad = max(d_sun_max_ad, max(d_helio));
        end
        max_sun_dist_AU(k) = d_sun_max_ad * Lstar / AU_km;   % maximum distance from the Sun [AU]

        % --- local refinement CLAMPED TO THE ARC of the minimum ---
        %   From the coarse minimum point, re-propagate at high resolution forward
        %   and backward, limiting the window to the time that stays INSIDE the arc:
        %   cosi' la propagazione (CR3BP puro) non attraversa mai un nodo DSM /
        %   flyby, dove la velocita' e' discontinua. Tutti gli archi sono
        %   campionati uniformemente (S2 resta uniforme dopo flipud), quindi
        %   the local index j directly gives the fraction of the arc traversed.
        S_min = arcs{a_best}(j_best, 1:6);
        Ta    = T_arc(a_best);
        dt_bwd = min(delta_t_nd, (j_best-1)      / (Narc-1) * Ta);  % time to the arc start
        dt_fwd = min(delta_t_nd, (Narc-j_best)   / (Narc-1) * Ta);  % time to the arc end

        if dt_fwd > 0
            [~, S_ref_post] = ode45(@(t,S) CR3BP(t,S,mu), ...
                                    linspace(0,  dt_fwd, Nrefine), S_min(:), opt);
            d_min_ad = min(d_min_ad, min(vecnorm(S_ref_post(:,1:3) - earth.', 2, 2)));
        end
        if dt_bwd > 0
            [~, S_ref_pre]  = ode45(@(t,S) CR3BP(t,S,mu), ...
                                    linspace(0, -dt_bwd, Nrefine), S_min(:), opt);
            d_min_ad = min(d_min_ad, min(vecnorm(S_ref_pre(:,1:3) - earth.', 2, 2)));
        end

        min_alt_km(k) = d_min_ad * Lstar - rEarth;            % altitude above Earth [km]

    catch ME
        if o.Verbose
            fprintf('[skip] solution %d: %s\n', k, ME.message);
        end
    end
end

%% ---- aggregazione -----------------------------------------------------
idx_valid = find(~isnan(min_alt_km));
alt       = min_alt_km(idx_valid);

res.min_alt_km    = min_alt_km;
res.idx_valid     = idx_valid;
res.n_valid       = numel(idx_valid);
res.n_total       = Nsol;
res.n_below_1000  = sum(alt < 1000);
res.n_below_10k   = sum(alt < 10000);
res.n_below_geo   = sum(alt < geo_alt_km);
res.n_below_100k  = sum(alt < 100000);
res.n_collision   = sum(alt < 0);
res.geo_alt_km    = geo_alt_km;

res.max_sun_dist_AU = max_sun_dist_AU;
if ~isempty(idx_valid)
    [res.max_sun_dist_overall_AU, i_sun_max] = max(max_sun_dist_AU(idx_valid));
    res.max_sun_dist_idx = idx_valid(i_sun_max);   % index in global_results
else
    res.max_sun_dist_overall_AU = NaN;
    res.max_sun_dist_idx = [];
end

%% ---- log --------------------------------------------------------------
if o.Verbose
    fprintf('\n============== EARTH ALTITUDE CHECK ==============\n');
    fprintf(' Solutions: %d total | %d computed\n', Nsol, numel(idx_valid));
    if ~isempty(alt)
        fprintf(' h_min [km]: min=%.1f | median=%.1f | max=%.1f\n', ...
            min(alt), median(alt), max(alt));
        fprintf(' Below 1000 km : %d | below 10k km : %d | below GEO : %d | below 100k km : %d\n', ...
            res.n_below_1000, res.n_below_10k, res.n_below_geo, res.n_below_100k);
        if res.n_collision > 0
            fprintf(2,' -> WARNING: %d trajectories IMPACT the Earth (h_min < 0).\n', ...
                res.n_collision);
        else
            fprintf(' -> No trajectory goes below the Earth surface.\n');
        end
        fprintf(' Maximum heliocentric distance reached: %.4f AU (solution index %d)\n', ...
            res.max_sun_dist_overall_AU, res.max_sun_dist_idx);
    end
    fprintf('=================================================\n\n');
end

%% ---- plot -------------------------------------------------------------
if o.Plot && ~isempty(alt)
    figure('Color','w','Name','Minimum Earth altitude distribution');
    hold on; box on; grid on;

    % Only positive altitudes for the log scale (collisions, h<0, cannot be
    % plotted in log: they are reported in the title/log).
    alt_pos = alt(alt > 0);
    lo = min([alt_pos; 1000]);
    hi = max([alt_pos; 100000]);
    edges = logspace(log10(lo) - 0.1, log10(hi) + 0.1, 40);

    histogram(alt_pos, edges, 'FaceColor',[0.30 0.55 0.85], ...
        'FaceAlpha',0.65, 'EdgeColor',[0.2 0.2 0.2], 'DisplayName','min altitude');

    % Reference lines
    xline(1000,       'Color',[0.85 0.15 0.15], 'LineStyle','-',  'LineWidth',1.8, ...
        'Label','1000 km',  'LabelVerticalAlignment','top', 'DisplayName','1000 km');
    xline(10000,      'Color',[0.95 0.55 0.10], 'LineStyle','-',  'LineWidth',1.6, ...
        'Label','10k km',   'LabelVerticalAlignment','top', 'DisplayName','10 000 km');
    xline(geo_alt_km, 'Color',[0.35 0.25 0.75], 'LineStyle','--', 'LineWidth',1.8, ...
        'Label','GEO',      'LabelVerticalAlignment','top', 'DisplayName',sprintf('GEO (%.0f km)',geo_alt_km));
    xline(100000,     'Color',[0.10 0.55 0.30], 'LineStyle','-',  'LineWidth',1.6, ...
        'Label','100k km',  'LabelVerticalAlignment','top', 'DisplayName','100 000 km');

    set(gca, 'XScale', 'log');
    xlabel('Minimum Earth altitude  [km]');
    ylabel('Number of solutions');
    if res.n_collision > 0
        title(sprintf('Minimum Earth altitude — %d solutions (%d impacts, h<0 not shown)', ...
            numel(idx_valid), res.n_collision));
    else
        title(sprintf('Minimum Earth altitude — %d solutions', numel(idx_valid)));
    end
    legend('Location','best');
    hold off;
end

end
