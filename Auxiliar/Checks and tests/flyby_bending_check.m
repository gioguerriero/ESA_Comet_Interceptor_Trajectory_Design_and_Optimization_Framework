function res = flyby_bending_check(global_results, c, varargin)
% flyby_bending_check  Feasibility check of the lunar flyby bending.
%
% For each global-search solution it computes:
%   - the ACTUAL lunar flyby bending (angle between the incoming and outgoing
%     v_inf at the Moon);
%   - the MAXIMUM feasible bending for that v_inf, at the minimum allowed
%     periapsis (r_p_min = rMoon + h_min);
%   - the DIFFERENCE (actual - maximum), in degrees.
%
% If the maximum bending is respected in all solutions, ALL differences must
% be <= 0. A difference > 0 flags a physically unrealizable flyby (it would
% require a periapsis below r_p_min).
%
% Physics (same conventions as bplane_from_vinf.m / verify_flyby_geometry.m):
%   v_inf_in  = vinf_rotation(moon_h2m, vinf_fmincon, fpa_in,  oop_in)   [km/s]
%   v_inf_out = vinf_rotation(moon_m2c, Vinf_m2c,     fpa_out, oop_out)  [km/s]
%   delta_actual = acos( <v_in,v_out> / (|v_in||v_out|) )          (actual bending)
%     v_inf = 0.5*(|v_in| + |v_out|)                                 (equal in patched conics)
%     e_min = 1 + r_p_min * v_inf^2 / muMoon                         (eccentricity at r_p_min)
%     delta_max = 2*asin(1/e_min)                                    (max bending at r_p_min)
%     diff = delta_actual - delta_max        (<= 0  => feasible flyby)
%
%   The angle between two vectors is rotation-invariant, so it can be computed
%   directly in the synodic frame (no need to switch to ECLIPJ2000).
%
% -------------------------------------------------------------------------
%  INPUT
% -------------------------------------------------------------------------
%   global_results - cell array (or struct array) of global-search solutions.
%                    Each solution sol must have:
%                       sol.x_opt              [fpa_in oop_in halo_th tof_moon tof_halo]
%                       sol.result_fmincon     [Vinf th fpa_out oop_out tof beta dv(1:3)]
%                       sol.fmincon_vinf_results.theta_moon / .vinf_fmincon
%   c              - constants struct (uses c.muMoon, c.rMoon; moon_state uses c)
%
%   Options ('Name',value):
%     'HminKm'  (750)   - minimum lunar flyby altitude [km] -> r_p_min = rMoon + HminKm
%     'Plot'    (true)  - draw the plot of the differences
%     'Verbose' (true)  - print the text summary
%
% -------------------------------------------------------------------------
%  OUTPUT (struct res)
% -------------------------------------------------------------------------
%   res.diff_deg          - difference (actual - maximum) per solution [deg]
%   res.delta_actual_deg  - actual bending per solution [deg]
%   res.delta_max_deg     - maximum feasible bending per solution [deg]
%   res.vinf_kms          - flyby |v_inf| per solution [km/s]
%   res.rp_actual_km      - periapsis implied by the actual bending [km]
%   res.idx_valid         - indices (in global_results) of the computed solutions
%   res.n_infeasible      - number of solutions with diff > 0 (bending NOT respected)
%   res.rp_min_km         - minimum periapsis used [km]
%
%   EXAMPLE (from the command window):
%     addpath('Auxiliar/'); addpath('Moon2Comet/');
%     S = load('Synthetic1_global_search.mat');   % -> global_results, c
%     res = flyby_bending_check(S.global_results, S.c);

%% ---- options ----------------------------------------------------------
p = inputParser;
p.addParameter('HminKm',  750);    % minimum lunar flyby altitude [km]
p.addParameter('Plot',    true);
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

muMoon = c.muMoon;                 % Moon GM [km^3/s^2]
rp_min = c.rMoon + o.HminKm;       % minimum allowed periapsis [km]

if isstruct(global_results), global_results = num2cell(global_results); end
Nsol = numel(global_results);

% Preallocation (NaN = solution not computable)
diff_deg         = nan(Nsol,1);
delta_actual_deg = nan(Nsol,1);
delta_max_deg    = nan(Nsol,1);
vinf_kms         = nan(Nsol,1);
rp_actual_km     = nan(Nsol,1);

for k = 1:Nsol
    try
        sol = global_results{k};
        xo  = sol.x_opt;                 % [fpa_in oop_in halo_th tof_moon tof_halo]
        rf  = sol.result_fmincon;        % [Vinf th fpa_out oop_out tof beta dv(1:3)]
        fvr = sol.fmincon_vinf_results;

        % --- incoming v_inf at the flyby (Halo->Moon branch) ---
        moon_h2m = moon_state(fvr.theta_moon, c)';                       % Moon state [6x1]
        vinf_in  = vinf_rotation(moon_h2m, fvr.vinf_fmincon, xo(1), xo(2));  % [3x1, km/s]

        % --- outgoing v_inf from the flyby (Moon->Comet branch) ---
        moon_m2c = moon_state(rf(2), c);                                 % Moon state [1x6]
        vinf_out = vinf_rotation(moon_m2c', rf(1), rf(3), rf(4));        % [3x1, km/s]

        vin_n  = norm(vinf_in);
        vout_n = norm(vinf_out);

        % --- actual bending = angle between the two v_inf ---
        cos_d = dot(vinf_in, vinf_out) / (vin_n * vout_n);
        delta_actual = acos(max(-1, min(1, cos_d)));      % [rad]

        % --- maximum feasible bending at the minimum periapsis ---
        vinf   = 0.5*(vin_n + vout_n);                    % flyby |v_inf| [km/s]
        e_min  = 1 + rp_min * vinf^2 / muMoon;            % eccentricity at r_p_min
        delta_max = 2*asin(1/e_min);                      % [rad]

        % --- periapsis implied by the actual bending (cross-check) ---
        %   r_p = (muMoon/vinf^2) * (1/sin(delta_actual/2) - 1)
        if delta_actual > 0
            rp_act = (muMoon / vinf^2) * (1/sin(delta_actual/2) - 1);
        else
            rp_act = Inf;   % no deflection -> infinite periapsis
        end

        % --- store ---
        delta_actual_deg(k) = rad2deg(delta_actual);
        delta_max_deg(k)    = rad2deg(delta_max);
        diff_deg(k)         = rad2deg(delta_actual - delta_max);
        vinf_kms(k)         = vinf;
        rp_actual_km(k)     = rp_act;

    catch ME
        if o.Verbose
            fprintf('[skip] solution %d: %s\n', k, ME.message);
        end
    end
end

%% ---- aggregation -----------------------------------------------------
idx_valid = find(~isnan(diff_deg));

res.diff_deg         = diff_deg;
res.delta_actual_deg = delta_actual_deg;
res.delta_max_deg    = delta_max_deg;
res.vinf_kms         = vinf_kms;
res.rp_actual_km     = rp_actual_km;
res.idx_valid        = idx_valid;
res.rp_min_km        = rp_min;
res.n_infeasible     = sum(diff_deg(idx_valid) > 0);
res.n_valid          = numel(idx_valid);
res.n_total          = Nsol;

% Maximum violation (diff>0) and index of the corresponding solution
viol_mask = diff_deg(idx_valid) > 0;
if any(viol_mask)
    viol_idx_in_valid          = idx_valid(viol_mask);
    [res.max_violation_deg, i_max] = max(diff_deg(viol_idx_in_valid));
    res.max_violation_idx      = viol_idx_in_valid(i_max);   % index in global_results
else
    res.max_violation_deg = NaN;
    res.max_violation_idx = [];
end

%% ---- log --------------------------------------------------------------
if o.Verbose
    dv = diff_deg(idx_valid);
    fprintf('\n============== FLYBY BENDING CHECK ==============\n');
    fprintf(' r_p_min = %.0f km (rMoon %.0f + h_min %.0f)\n', rp_min, c.rMoon, o.HminKm);
    fprintf(' Solutions: %d total | %d computed\n', Nsol, numel(idx_valid));
    if ~isempty(dv)
        fprintf(' diff (actual - maximum) [deg]: min=%.2f | max=%.2f | mean=%.2f\n', ...
            min(dv), max(dv), mean(dv));
        fprintf(' Solutions with bending NOT respected (diff>0): %d / %d\n', ...
            res.n_infeasible, numel(idx_valid));
        if res.n_infeasible == 0
            fprintf(' -> OK: the maximum bending is respected in all solutions.\n');
        else
            fprintf(2,' -> ATTENZIONE: %d flyby non fattibili (diff>0).\n', res.n_infeasible);
            fprintf(2,' -> Maximum violation: +%.2f deg (solution index %d in global_results)\n', ...
                res.max_violation_deg, res.max_violation_idx);
        end
    end
    fprintf('================================================\n\n');
end

%% ---- plot -------------------------------------------------------------
if o.Plot && ~isempty(idx_valid)
    dv = diff_deg(idx_valid);
    figure('Color','w','Name','Flyby bending: actual - maximum');
    hold on; box on; grid on;

    % Feasible solutions (diff<=0) in green, infeasible (diff>0) in red
    feas = dv <= 0;
    scatter(find(feas),  dv(feas),  28, [0.10 0.60 0.30], 'filled', ...
        'DisplayName','feasible (diff \leq 0)');
    if any(~feas)
        scatter(find(~feas), dv(~feas), 34, [0.85 0.15 0.15], 'filled', ...
            'DisplayName','NOT feasible (diff > 0)');
    end

    % Zero threshold line (feasibility boundary)
    yline(0, 'k--', 'LineWidth', 1.4, 'DisplayName','feasibility threshold (0)');

    xlabel('Solution (index among valid ones)');
    ylabel('\Delta bending = \delta_{actual} - \delta_{max}  [deg]');
    if res.n_infeasible > 0
        title(sprintf('Flyby bending check - %d/%d infeasible flybys (h_{min}=%.0f km, max violation = +%.2f deg)', ...
            res.n_infeasible, numel(idx_valid), o.HminKm, res.max_violation_deg));
    else
        title(sprintf('Flyby bending check — %d/%d infeasible flybys (h_{min}=%.0f km)', ...
            res.n_infeasible, numel(idx_valid), o.HminKm));
    end
    legend('Location','best');
    hold off;
end

end
