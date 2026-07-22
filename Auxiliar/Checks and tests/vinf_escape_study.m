function res = vinf_escape_study(global_results, c, S_halo, unstable_dir, ...
                                 eps_vel_ms, control_radius, selected_comet, varargin)
% vinf_escape_study  Ballistic escape v-infinity at the control radius.
%
%   For each trajectory in the dataset (Monte Carlo / global search) it
%   reconstructs the CR3BP arcs, removes from the final HELIOCENTRIC energy (at
%   the control radius rc*) the ENERGY injected by the propulsive maneuvers
%   (injection, DSM1, and DSM2 if performed inside rc*), and returns the
%   "ballistic" v_inf (v_inf_ball), to compare against reference values:
%
%        v_inf,fb (Rc < 1 AU) = 1044.2 m/s      (v_inf,dir = 816.9 -> +227.3)
%        v_inf,fb (Rc > 1 AU) =  983.7 m/s      (v_inf,dir = 825.2 -> +158.5)
%
%   It does NOT isolate the lunar flyby alone: it isolates the NON-propulsive
%   part of the trajectory (ballistic injection + coast + lunar flyby), exactly
%   like the end-to-end ballistic v_inf,fb of Ranuschio. The (ballistic) lunar
%   flyby contribution is included: eps_f is measured AFTER the flyby, at rc*.
%
% =========================================================================
%  DEFINITION OF v_inf (HELIOCENTRIC, as Ranuschio Eq. 4.3-4.4, p.53-54)
% =========================================================================
%   Ranuschio computes a HELIOCENTRIC v_inf (Sun-centred inertial frame), NOT
%   geocentric. It propagates out to rc* = 0.2 (non-dim = 0.2 AU) from Earth,
%   where the motion is by then Sun-dominated and the osculating heliocentric
%   semi-major axis "a" has converged (Fig. 4.3). Then:
%
%       eps_helio = 0.5*|V_helio|^2 - mu_sun/|r_helio|        (specific energy)
%       a         = -mu_sun / (2*eps_helio)                   (osculating semi-major axis)
%       v_inf     = sqrt( 2*mu_sun/rE - mu_sun/a ) - sqrt(mu_sun/rE)   (Eq. 4.4)
%
%   with rE = 1 AU (fixed reference: heliocentric velocity evaluated at Earth's
%   orbit). Sign: v_inf>0 => a>1AU (outer target), v_inf<0 => a<1AU (inner
%   target). Ranuschio reports the MAGNITUDE in Tab. 5.2.
%
%   In CR3BP non-dimensional units (Lstar=1AU, G(mSun+mEarth)=1):
%       mu_sun = 1 - mu   (verified in CR3BP.m),   rE = 1,
%       sqrt(mu_sun/rE) = Earth circular velocity ~ 1 (~Vstar).
%
% =========================================================================
%  FRAME (heliocentric velocities in the synodic axes)
% =========================================================================
%   Everything is in the non-dimensional SYNODIC axes, but the velocities used
%   for the energy are INERTIAL relative to the Sun:
%       r_helio = r_syn - [-mu,0,0] = r_syn + [mu,0,0]     (Sun at [-mu,0,0])
%       V_helio = v_syn + omega x r_helio,   omega = [0,0,1] adim
%             = [ v_syn(1)-r_helio(2);  v_syn(2)+r_helio(1);  v_syn(3) ]
%   Norms and dot products are invariant under the rigid rotation synodic <->
%   inertial axes: no need to rotate to a fixed frame. The impulsive Delta-v is
%   the same vector in both frames (position is continuous).
%
% =========================================================================
%  MANEUVER CLEANING (heliocentric energy, Oberth included)
% =========================================================================
%   Increment of heliocentric specific energy from an impulsive burn:
%       d_eps_i = dot(V_helio_pre_i, dv_i) + 0.5*norm(dv_i)^2
%   (V_helio_pre ~ Earth orbital velocity => strong leverage on the heliocentric
%    energy: it is exactly the boost we want to remove.)
%
%       eps_ball  = eps_f - sum_i(d_eps_i)      (only maneuvers "before" the escape)
%       v_inf_ball = Eq.4.4(eps_ball)
%
%   APPROXIMATION: removing the burns would change the coasts; here the standard
%   energy balance (sum of impulsive jumps) is used, consistent with Ranuschio's
%   ballistic v_inf,fb. Natural coast drift and flyby jump: included.
%
% =========================================================================
%  IN/OUT OF THE CONTROL RADIUS (temporal rule, not only spatial)
% =========================================================================
%   - Injection and DSM1 are on the Halo->flyby leg (inbound): they ALWAYS occur
%     before the final escape => always subtracted.
%   - DSM2 is on the flyby->comet leg (outbound): subtracted ONLY if the last
%     outward crossing of rc* falls AFTER DSM2 (crossing on the DSM2->comet arc).
%     This keeps the measured eps_f and the subtracted DSM2 always consistent.
%   The crossing is detected on the GEOCENTRIC distance (|r_syn - [1-mu,0,0]| = rc*).
%
% -------------------------------------------------------------------------
%  INPUT
% -------------------------------------------------------------------------
%   global_results - cell array (or struct array) of global-search solutions.
%                    Each solution sol must have:
%                       sol.x_opt              [fpa oop halo_th tof_moon tof_halo]
%                       sol.result_fmincon     [Vinf th fpa oop tof beta dv(1:3)]
%                       sol.fmincon_vinf_results.theta_moon / .vinf_fmincon
%   c              - constants struct (c.mu, c.Vstar, c.Lstar, c.rMoon_ad, ...)
%   S_halo         - Nx6 halo states (synodic non-dim)
%   unstable_dir   - Nx3 unstable directions (velocity) on the halo
%   eps_vel_ms     - injection magnitude from the manifold [m/s]  (search_params.eps_vel_ms)
%   control_radius - NON-DIMENSIONAL control radius rc* (=/Lstar). Ranuschio: 0.2
%   selected_comet - comet struct (region tag). .comet_pos [km, heliocentric]
%
%   Options ('Name',value):
%     'Plot'(true) 'NArc'(600) 'AbsTol'(1e-10) 'RelTol'(1e-10) 'Verbose'(true)
%     'RefFbLt1AU'(1044.2) 'RefFbGt1AU'(983.7)      % v_inf,fb Ranuschio [m/s]
%     'RefDirLt1AU'(816.9) 'RefDirGt1AU'(825.2)     % v_inf,dir Ranuschio [m/s]
%
% -------------------------------------------------------------------------
%  OUTPUT (struct res)
% -------------------------------------------------------------------------
%   res.vinf_ball_ms   - ballistic v_inf (SIGNED) per valid solution [m/s]
%   res.vinf_ball_abs_ms - magnitude (for direct comparison with Ranuschio)
%   res.idx_valid, res.mean_ms, res.std_ms, res.mean_abs_ms, res.std_abs_ms
%   res.region ('Rc<1AU'|'Rc>1AU'), res.ref_fb_ms, res.ref_dir_ms, res.R_comet_AU
%   res.per_sol        - per-trajectory diagnostics (all)
%   res.anomalies      - non-real v_inf_ball cases (radicand<0): count, %, ...
%
%   Example:
%     S = load('Synthetic1_global_search.mat');   % global_results, c, selected_comet, search_params
%     load('S_halo.mat'); load('unstable_dir.mat');
%     res = vinf_escape_study(S.global_results, S.c, S_halo, unstable_dir, ...
%                 S.search_params.eps_vel_ms, 0.2, S.selected_comet);

%% ---- options ------------------------------------------------------------
% inputParser handles the optional ('Name', value) pairs passed at the end of
% the function (varargin). Each addParameter defines a default: if the user does
% not specify that name, the value indicated here is used.
p = inputParser;
p.addParameter('Plot',    true);      % true -> draw the final PDF of v_inf_ball
p.addParameter('NArc',    600);       % number of ode45 samples for each of the 4 arcs
p.addParameter('AbsTol',  1e-10);     % ode45 absolute integration tolerance
p.addParameter('RelTol',  1e-10);     % ode45 relative integration tolerance
p.addParameter('RefFbLt1AU', 1044.2); % v_inf,fb Ranuschio for comets with Rc < 1 AU [m/s]
p.addParameter('RefFbGt1AU',  983.7); % v_inf,fb Ranuschio for comets with Rc > 1 AU [m/s]
p.addParameter('RefDirLt1AU', 816.9); % v_inf,dir Ranuschio (direct escape) Rc < 1 AU [m/s]
p.addParameter('RefDirGt1AU', 825.2); % v_inf,dir Ranuschio (direct escape) Rc > 1 AU [m/s]
p.addParameter('RemoveManeuvers', true);      % true  = v_inf,ball (maneuver energy
%   removed, comparable with Ranuschio). false = "raw" v_inf at the control
%   radius, maneuvers INCLUDED (no cleaning): computes only the v_inf as-is.
p.addParameter('SubtractInjection', true);   % false = keep the injection in the
%   ballistic baseline (the 15 m/s manifold injection is analogous to the
%   "quasi-ballistic departure" that Ranuschio does NOT remove). Effect <~15 m/s.
%   Ignored if RemoveManeuvers = false.
p.addParameter('Verbose', true);      % true -> print the text summary at the end of the function
p.parse(varargin{:});
o = p.Results;    % struct with all resolved options (default or user-passed)

%% ---- constants -------------------------------------------------------
mu     = c.mu;                   % Sun-Earth mass ratio (non-dimensional, CR3BP)
mu_sun = 1 - mu;                 % heliocentric gravitational parameter [non-dim];
                                  % in the CR3BP the Sun "weighs" (1-mu), the Earth "mu"
Vstar  = c.Vstar;                % characteristic CR3BP velocity [km/s], used to
                                  % re-dimensionalize v_inf from non-dim to m/s
rE     = 1;                      % Earth orbit [non-dim] = 1 AU (by CR3BP definition)
earth  = [1-mu; 0; 0];           % Earth position in the synodic frame -> used to
                                  % detect the control-radius crossing
sun    = [-mu;  0; 0];           % Sun position in the synodic frame -> used as
                                  % origin to compute the heliocentric energy
eps_vel_ad = eps_vel_ms / (1e3 * Vstar);   % injection magnitude: from m/s to non-dim
opt    = odeset('AbsTol', o.AbsTol, 'RelTol', o.RelTol);   % common ode45 options
Narc   = o.NArc;                 % resolution (number of points) of each propagated arc
rc     = control_radius;         % control radius rc* [non-dim] (Ranuschio: 0.2)

% The rest of the function works on a cell array: if the user passes a struct
% array it is converted here once and for all.
if isstruct(global_results), global_results = num2cell(global_results); end
Nsol = numel(global_results);    % total number of trajectories/solutions to process

%% ---- region tag (comet heliocentric distance) ------------------------
% Ranuschio distinguishes two target categories: comets with an encounter at a
% heliocentric distance below or above 1 AU, with different reference values.
% Here the Sun-comet distance is computed and the category (and thus the
% reference values) to compare against is selected automatically.
R_comet_AU = norm(selected_comet.comet_pos(:)) / c.Lstar;
if R_comet_AU < 1
    region = 'Rc<1AU';  ref_fb = o.RefFbLt1AU;  ref_dir = o.RefDirLt1AU;
else
    region = 'Rc>1AU';  ref_fb = o.RefFbGt1AU;  ref_dir = o.RefDirGt1AU;
end

%% ---- loop over the solutions ---------------------------------------------
% Preallocate the diagnostics struct array, one entry per dataset trajectory
% (valid or not): assigning the last element (Nsol) with these fields makes
% MATLAB automatically extend the array to that length.
per_sol(Nsol) = struct('valid',false,'reason','','vinf_ball_ms',NaN, ...
    'vinf_f_ms',NaN,'a_ball_AU',NaN,'a_f_AU',NaN,'sum_d_eps_ad',NaN, ...
    'eps_f_ad',NaN,'eps_ball_ad',NaN,'maneuvers',[],'dsm2_subtracted',false);

for k = 1:Nsol
    % try/catch: if a single trajectory fails (ode45 integration, missing
    % fields, etc.) it is skipped and the loop continues with the others instead
    % of crashing the whole function. The failure reason is stored in .reason.
    try
        sol = global_results{k};
        xo  = sol.x_opt;                 % design variables of the Halo->flyby leg:
                                          % [fpa oop halo_th tof_moon tof_halo]
        rf  = sol.result_fmincon;        % design variables of the flyby->comet leg:
                                          % [Vinf th fpa oop tof beta dv(1:3)]
        fvr = sol.fmincon_vinf_results;  % lunar flyby geometry (theta_moon, vinf_fmincon, ...)

        % --- Halo -> DSM1 -> flyby leg variables ---
        ref_fpa = xo(1);  ref_oop = xo(2);  ref_halo_th = xo(3);
        ref_tof_moon = xo(4);  ref_tof_halo = xo(5);
        % --- flyby -> DSM2 -> comet leg variables ---
        Vinf_m2c = rf(1); theta_m2c = rf(2); fpa_m2c = rf(3);
        oop_m2c  = rf(4); tof_m2c   = rf(5); beta_m2c = rf(6);
        dv_dsm2  = rf(7:9); dv_dsm2 = dv_dsm2(:);   % DSM2 Delta-v vector (synodic non-dim)

        % ================= CR3BP arc reconstruction ======================
        % The solution stores only the design variables (angles, times of
        % flight, etc.), not the states along the trajectory: here the 4 arcs
        % (Halo->DSM1->flyby->DSM2->comet) are re-propagated exactly as in
        % plot_paper_trajectory.m, to recover position/velocity at each maneuver
        % and along the whole path.

        % --- Arc 1: Halo departure (with injection) -> DSM1 ---
        % Find the Halo-orbit point closest to angle ref_halo_th and its
        % unstable direction (manifold eigenvector).
        [init_st, idx_h] = state_finder(ref_halo_th, S_halo);
        init_st = init_st(:);
        vu = unstable_dir(idx_h,:)';  vu = vu / norm(vu);   % normalized unstable direction
        dv_inj    = eps_vel_ad * vu;         % injection Delta-v (small push along the manifold)
        r_inj     = init_st(1:3);            % position at injection
        v_inj_pre = init_st(4:6);            % velocity BEFORE injection (on the Halo orbit)
        S0_arc1   = [r_inj; v_inj_pre + dv_inj];   % arc-1 initial state (after injection)
        [~, S1] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, ref_tof_halo, Narc), S0_arc1, opt);

        % --- Arc 2: DSM1 -> lunar flyby ---
        % This arc is reconstructed by propagating BACKWARD in time from the
        % lunar-flyby state (directly parameterized by the design variables),
        % then flipped (flipud) to obtain the arc in the correct direction
        % (DSM1 -> flyby).
        moon_h2m = moon_state(fvr.theta_moon, c)';           % Moon state at the flyby instant
        vinf_h2m = vinf_rotation(moon_h2m, fvr.vinf_fmincon, ref_fpa, ref_oop);  % incoming v_inf at the flyby
        S0_flyby = moon_h2m;  S0_flyby(4:6) = S0_flyby(4:6) + vinf_h2m ./ Vstar;
        [~, S2b] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, -ref_tof_moon, Narc), S0_flyby, opt);
        S2 = flipud(S2b);    % now S2 goes from DSM1 (start) to flyby (end), correct time order

        % --- DSM1: maneuver straddling arc 1 and arc 2 ---
        % The position at the end of arc1 and start of arc2 coincides (same
        % physical point); the DSM1 Delta-v is the velocity "jump" needed to
        % join the two separately propagated arcs.
        r_dsm1     = S1(end,1:3).';
        v_dsm1_pre = S1(end,4:6).';                       % velocity BEFORE DSM1 (end of arc1)
        dv_dsm1    = (S2(1,4:6) - S1(end,4:6)).';          % Delta-v = v_after (start of arc2) - v_before

        % --- Arc 3: flyby lunare -> DSM2 ---
        moon_m2c = moon_state(theta_m2c, c);                % Moon state at the flyby instant (outbound)
        vinf_m2c = vinf_rotation(moon_m2c', Vinf_m2c, fpa_m2c, oop_m2c);   % outgoing v_inf from the flyby
        S0_arc3  = [moon_m2c(1:3)'; moon_m2c(4:6)' + vinf_m2c ./ Vstar];
        [~, S3]  = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, tof_m2c*beta_m2c, Narc), S0_arc3, opt);

        % --- DSM2: maneuver at the end of arc 3 ---
        % Unlike DSM1 (obtained as the difference between two arcs), DSM2 is a
        % direct output of the solution (dv_dsm2 = rf(7:9)), the result of a
        % Lambert arc computed during the optimization.
        r_dsm2     = S3(end,1:3).';
        v_dsm2_pre = S3(end,4:6).';                        % velocity BEFORE DSM2 (end of arc3)

        % --- Arc 4: DSM2 -> comet ---
        S0_arc4 = [S3(end,1:3)'; S3(end,4:6)' + dv_dsm2];   % velocity AFTER DSM2
        [~, S4] = ode45(@(t,S) CR3BP(t,S,mu), linspace(0, tof_m2c*(1-beta_m2c), Narc), S0_arc4, opt);

        % ================= control-radius crossing (geocentric) ======
        % Look for the LAST "outward" crossing (Earth distance rising past rc*)
        % along arc 4; if none exists (the spacecraft was already outside rc*
        % before DSM2), look on arc 3. This also determines whether DSM2 occurred
        % INSIDE the control radius (dsm2_before_exit = true, arc4) or outside
        % (false, the exit happened earlier, on arc3).
        [state_cross_4, ok4] = last_outward_crossing(S4, earth, rc);
        if ok4
            state_cross = state_cross_4;  dsm2_before_exit = true;
        else
            [state_cross_3, ok3] = last_outward_crossing(S3, earth, rc);
            if ok3
                state_cross = state_cross_3;  dsm2_before_exit = false;
            else
                % No crossing found: the trajectory never leaves the control
                % radius (an anomalous/non-physical case for this study).
                per_sol(k).reason = 'no outward crossing of rc*';
                continue
            end
        end

        % eps_f = "osculating" (instantaneous) HELIOCENTRIC specific energy at
        % the exact point where the trajectory crosses rc*. It is the "raw"
        % energy (maneuvers still included), the equivalent of Ranuschio's
        % Eq. 4.3 but applied to OUR trajectory.
        rh_f  = state_cross(1:3) - sun;                 % position relative to the Sun
        Vh_f  = v_helio_inertial(state_cross(4:6), rh_f);   % inertial heliocentric velocity
        eps_f = 0.5*norm(Vh_f)^2 - mu_sun / norm(rh_f);     % specific energy = kinetic + potential

        % ================= maneuver energy contributions ============
        % Decide which maneuvers to subtract from the energy:
        %  - if RemoveManeuvers=false, subtract NOTHING ("raw" v_inf with all
        %    maneuvers included: the "as-is" comparison mode).
        %  - injection: subtracted only if RemoveManeuvers=true AND SubtractInjection=true.
        %  - DSM1: always on the inbound leg, subtracted whenever
        %    RemoveManeuvers=true (it always occurs before the final exit).
        %  - DSM2: subtracted only if RemoveManeuvers=true AND DSM2 occurred
        %    before leaving the control radius (dsm2_before_exit).
        sub_inj  = o.RemoveManeuvers && o.SubtractInjection;
        sub_dsm1 = o.RemoveManeuvers;
        sub_dsm2 = o.RemoveManeuvers && dsm2_before_exit;
        % mk_man computes, for each maneuver, the heliocentric energy jump d_eps
        % (Oberth formula) and records whether it should be subtracted.
        man(1) = mk_man('inj',  r_inj,  v_inj_pre,  dv_inj,  sub_inj,  sun, earth, rc, mu_sun);
        man(2) = mk_man('dsm1', r_dsm1, v_dsm1_pre, dv_dsm1, sub_dsm1, sun, earth, rc, mu_sun);
        man(3) = mk_man('dsm2', r_dsm2, v_dsm2_pre, dv_dsm2, sub_dsm2, sun, earth, rc, mu_sun);

        % Sum only the d_eps of the maneuvers marked "subtract" (the others, if
        % RemoveManeuvers=false or if they do not meet the conditions above,
        % contribute 0 to the sum).
        sum_d_eps = 0;
        for m = 1:numel(man)
            if man(m).subtract, sum_d_eps = sum_d_eps + man(m).d_eps; end
        end
        eps_ball = eps_f - sum_d_eps;    % "clean" energy (equivalent ballistic)

        % ================= Eq. 4.4: energia -> v_inf ======================
        % Convert both the raw energy (eps_f) and the ballistic one (eps_ball)
        % to v_inf with Ranuschio's formula (Eq. 4.4), so both quantities are
        % computed identically and are comparable.
        [vinf_f_ad,    a_f_ad]    = vinf_from_eps(eps_f,    mu_sun, rE);
        [vinf_ball_ad, a_ball_ad] = vinf_from_eps(eps_ball, mu_sun, rE);

        % Store all the intermediate diagnostics for this trajectory, useful to
        % inspect single or anomalous cases afterwards.
        per_sol(k).eps_f_ad    = eps_f;
        per_sol(k).eps_ball_ad = eps_ball;
        per_sol(k).sum_d_eps_ad= sum_d_eps;
        per_sol(k).maneuvers   = man;
        per_sol(k).dsm2_subtracted = dsm2_before_exit;
        per_sol(k).a_f_AU      = a_f_ad;      % "raw" osculating semi-major axis (non-dim = AU)
        per_sol(k).a_ball_AU   = a_ball_ad;   % "ballistic" osculating semi-major axis
        if isreal(vinf_f_ad) && isfinite(vinf_f_ad)
            per_sol(k).vinf_f_ms = vinf_f_ad * Vstar * 1e3;   % raw v_inf, from non-dim to m/s
        end

        % The trajectory is "valid" only if vinf_ball_ad is a real and finite
        % number: vinf_from_eps returns NaN when the Eq.4.4 radicand is negative
        % (the orbit with energy eps_ball would never reach 1 AU, a non-physical
        % situation for this comparison).
        if isreal(vinf_ball_ad) && isfinite(vinf_ball_ad)
            per_sol(k).valid = true;
            per_sol(k).vinf_ball_ms = vinf_ball_ad * Vstar * 1e3;   % ballistic v_inf, from non-dim to m/s
        else
            per_sol(k).reason = 'v_inf_ball not real (Eq.4.4 radicand < 0)';
        end

    catch ME
        % Any unexpected error (e.g. ode45 not converging, missing fields in the
        % sol struct, etc.) is caught here: trajectory k stays marked
        % invalid=false with the error message in .reason.
        per_sol(k).reason = sprintf('reconstruction error: %s', ME.message);
    end
end

%% ---- aggregation -----------------------------------------------------
% Extract only the valid trajectories (valid_mask) and build the final results
% vector on which the summary statistics are computed.
valid_mask   = [per_sol.valid];
vinf_ball_ms = [per_sol(valid_mask).vinf_ball_ms];
idx_valid    = find(valid_mask);

res.vinf_ball_ms     = vinf_ball_ms(:);          % SIGNED ballistic v_inf [m/s] (column)
res.vinf_ball_abs_ms = abs(vinf_ball_ms(:));     % magnitude, for direct comparison with Ranuschio
res.idx_valid        = idx_valid(:);             % indices (in global_results) of the valid solutions
res.mean_ms          = mean(vinf_ball_ms);       % signed mean
res.std_ms           = std(vinf_ball_ms);        % signed standard deviation
res.mean_abs_ms      = mean(abs(vinf_ball_ms));  % mean of the magnitude -> compare with v_inf,fb
res.std_abs_ms       = std(abs(vinf_ball_ms));   % standard deviation of the magnitude
res.region           = region;                   % 'Rc<1AU' or 'Rc>1AU'
res.ref_fb_ms         = ref_fb;                  % Ranuschio reference value (flyby) used
res.ref_dir_ms        = ref_dir;                 % Ranuschio reference value (direct) used
res.R_comet_AU        = R_comet_AU;              % heliocentric distance of the target comet [AU]
res.control_radius_ad = rc;                      % control radius used [non-dim]
res.per_sol           = per_sol;                 % full diagnostics, ALL trajectories

% ---- anomaly report (Eq.4.4 radicand < 0) ---------------------------
% Count separately how many trajectories are "invalid" for the specific reason
% "negative radicand" (maneuvers dominating the energy balance), distinguishing
% them from other failure types (numerical errors, no crossing found, etc.), so
% they can be inspected separately.
is_anom  = arrayfun(@(s) ~s.valid && contains(s.reason,'radicand'), per_sol);
anom_idx = find(is_anom);
res.anomalies.count   = numel(anom_idx);
res.anomalies.percent = 100 * numel(anom_idx) / max(Nsol,1);
res.anomalies.idx     = anom_idx(:);
res.anomalies.sum_d_eps_ad = arrayfun(@(s) s.sum_d_eps_ad, per_sol(anom_idx));
res.anomalies.eps_f_ad     = arrayfun(@(s) s.eps_f_ad,     per_sol(anom_idx));

res.n_total   = Nsol;                     % total number of trajectories processed
res.n_valid   = numel(idx_valid);         % number of trajectories with a valid v_inf_ball
res.n_invalid = Nsol - numel(idx_valid);  % all the others (anomalies + errors + no crossing)

%% ---- log --------------------------------------------------------------
% Print a readable text summary to screen (only if Verbose=true): mode used,
% comet/region, Ranuschio reference values, control radius, validity counts,
% final statistics, and the deviation vs Ranuschio.
if o.Verbose
    if o.RemoveManeuvers, mode_str = 'v_inf,ball (maneuvers REMOVED)';
    else,                 mode_str = 'raw v_inf (maneuvers INCLUDED)'; end
    fprintf('\n================ VINF ESCAPE STUDY (eliocentrico, Eq.4.4) ================\n');
    fprintf(' Mode: %s\n', mode_str);
    fprintf(' Comet: %s | R_comet = %.3f AU -> region %s\n', ...
        get_name(selected_comet), R_comet_AU, region);
    fprintf(' Ranuschio: v_inf,fb = %.1f m/s | v_inf,dir = %.1f m/s | benefit flyby = %+.1f m/s\n', ...
        ref_fb, ref_dir, ref_fb - ref_dir);
    fprintf(' Control radius rc* = %.4f adim = %.0f km = %.3f AU\n', rc, rc*c.Lstar, rc);
    fprintf(' Solutions: %d total | %d valid | %d invalid (%d anomalies)\n', ...
        Nsol, res.n_valid, res.n_invalid, res.anomalies.count);
    if res.n_valid > 0
        fprintf(' |v_inf_ball| : mean = %.1f m/s | std = %.1f m/s | [%.1f, %.1f]\n', ...
            res.mean_abs_ms, res.std_abs_ms, min(abs(vinf_ball_ms)), max(abs(vinf_ball_ms)));
        fprintf(' Delta vs Ranuschio v_inf,fb = %+.1f m/s\n', res.mean_abs_ms - ref_fb);
    end
    fprintf('=========================================================================\n\n');
end

%% ---- plot PDF ---------------------------------------------------------
% Draw the histogram/PDF of |v_inf_ball| with the Ranuschio reference lines
% overlaid (only if Plot=true and there are valid solutions).
if o.Plot && res.n_valid > 0
    plot_vinf_pdf(abs(vinf_ball_ms), ref_fb, ref_dir, region, ...
                  get_name(selected_comet), res.mean_abs_ms);
end

end % ===================== end of main function ======================


%% ======================================================================
%  SUBFUNCTIONS
%  ======================================================================

function Vh = v_helio_inertial(v_syn, r_helio)
% v_helio_inertial  Convert a synodic velocity to an inertial velocity relative
% to the Sun, WITHOUT changing axes (the synodic axes remain, only the "point of
% view" becomes inertial). Adding the frame-rotation drag term of the synodic
% frame (omega = [0;0;1]):
%   V_helio = v_syn + omega x r_helio
% Expanding the cross product with omega=[0;0;1] gives the component-by-component
% form used below (faster than a generic cross()).
    v_syn = v_syn(:);  r_helio = r_helio(:);
    Vh = [ v_syn(1) - r_helio(2);
           v_syn(2) + r_helio(1);
           v_syn(3) ];
end

function [vinf, a] = vinf_from_eps(eps, mu_sun, rE)
% vinf_from_eps  Ranuschio Eq. 4.4: from the heliocentric specific energy eps,
% first obtain the osculating semi-major axis "a" (orbital-energy definition:
% eps = -mu_sun/(2a)), then v_inf as "the heliocentric velocity that the orbit
% of semi-major axis a would have at distance rE, minus Earth's circular
% velocity at rE" (i.e. the velocity excess over Earth's natural motion).
% Sign: vinf>0 if a>rE (outbound orbit, outer target),
%       vinf<0 if a<rE (inner orbit, inner target).
    a = -mu_sun / (2*eps);                 % osculating semi-major axis [non-dim=AU]
    disc = 2*mu_sun/rE - mu_sun/a;         % = V_helio(rE)^2 (vis-viva valutata a rE)
    if disc < 0
        % If the radicand is negative, the orbit with energy "eps" never passes
        % through rE=1AU (physically it would not reach Earth's orbit): return
        % NaN and the trajectory will be marked as an anomaly.
        vinf = NaN;  return
    end
    vinf = sqrt(disc) - sqrt(mu_sun/rE);
end

function m = mk_man(name, r_pre, v_pre, dv, subtract, sun, earth, rc, mu_sun) %#ok<INUSD>
% mk_man  Build the diagnostics struct of ONE impulsive maneuver and compute
% its heliocentric energy contribution d_eps (Oberth-effect formula):
% d_eps = dot(V_helio_pre, dv) + 0.5*|dv|^2. This is NOT just the Delta-v
% magnitude: it also depends on how fast (and in which direction) the spacecraft
% was travelling BEFORE the maneuver.
    r_pre = r_pre(:);  v_pre = v_pre(:);  dv = dv(:);
    r_helio = r_pre - sun;                        % position relative to the Sun at the burn
    Vh      = v_helio_inertial(v_pre, r_helio);    % heliocentric velocity BEFORE the maneuver
    d_eps   = dot(Vh, dv) + 0.5*norm(dv)^2;        % heliocentric specific-energy jump [non-dim]
    r_geo   = norm(r_pre - earth);                 % Earth distance (diagnostic only, to tell
                                                    % whether the maneuver occurred inside or outside rc*)
    m.name     = name;          % label ('inj' / 'dsm1' / 'dsm2')
    m.r_geo_ad = r_geo;         % geocentric distance at the burn [non-dim]
    m.inside   = r_geo < rc;    % spatial flag: inside the control radius? (diagnostic only)
    m.dv_ad    = norm(dv);      % Delta-v magnitude [non-dim]
    m.d_eps    = d_eps;         % energy contribution computed above [non-dim]
    m.subtract = subtract;      % flag: should this maneuver be subtracted from the final balance?
end

function [state_cross, ok] = last_outward_crossing(S, earth, rc)
% last_outward_crossing  Find, along a propagated arc S (Nx6, synodic), the
% LAST instant where the Earth distance goes from BELOW to ABOVE rc (i.e. an
% "outward" crossing of the control radius). If there are several along the arc,
% take the latest one (cross(end)): it is the physically relevant event defining
% the "final exit". The exact state at radius rc is estimated by linear
% interpolation between the sample just before and just after the crossing (the
% velocity is continuous within the arc, so interpolation is a good
% approximation when the samples are dense enough).
    rgeo  = vecnorm(S(:,1:3) - earth.', 2, 2);    % Earth distance at each arc sample
    cross = find(rgeo(1:end-1) < rc & rgeo(2:end) >= rc);   % indices i with a crossing between samples i and i+1
    if isempty(cross), state_cross = [];  ok = false;  return; end
    i = cross(end);                                % last (latest) crossing found
    f = (rc - rgeo(i)) / (rgeo(i+1) - rgeo(i));     % interpolation fraction between the two samples
    state_cross = (S(i,:) + f*(S(i+1,:) - S(i,:))).';   % state (6x1) interpolated exactly at rc
    ok = true;
end

function name = get_name(selected_comet)
% get_name  Extract the comet name from the struct (for titles/legends), with a
% default value if the .name field is absent.
    if isfield(selected_comet,'name'), name = char(string(selected_comet.name));
    else, name = 'comet'; end
end

function plot_vinf_pdf(vinf_abs_ms, ref_fb, ref_dir, region, cname, mean_ms)
% plot_vinf_pdf  Draw the normalized histogram (PDF) of |v_inf_ball| over the
% valid samples, with an overlaid kernel-density estimate (if ksdensity is
% available) and three vertical reference lines: Ranuschio's v_inf,fb,
% Ranuschio's v_inf,dir, and the mean of our data.
    figure('Color','w','Name',sprintf('v_{inf,ball} PDF — %s',cname));
    hold on; box on; grid on;

    % x-axis range fitted to both the data and the reference values, with a
    % 40 m/s margin per side.
    lo = min([vinf_abs_ms(:); ref_fb; ref_dir]) - 40;
    hi = max([vinf_abs_ms(:); ref_fb; ref_dir]) + 40;
    edges = linspace(lo, hi, 40);

    histogram(vinf_abs_ms, edges, 'Normalization','pdf', ...
        'FaceColor',[0.30 0.55 0.85], 'FaceAlpha',0.55, 'EdgeColor',[0.2 0.2 0.2]);

    % Kernel-density estimate (continuous curve), only if the function exists
    % (requires the Statistics Toolbox); wrapped in try/catch for safety.
    if exist('ksdensity','file')
        xg = linspace(lo, hi, 400);
        try, plot(xg, ksdensity(vinf_abs_ms, xg), 'Color',[0.10 0.30 0.60], 'LineWidth',1.8); catch, end
    end

    % Vertical reference lines (v_inf,fb, v_inf,dir, our mean), drawn across the
    % current plot height (yl).
    yl = ylim;
    h_fb  = plot([ref_fb  ref_fb ], yl, 'r-',  'LineWidth',2.2);
    h_dir = plot([ref_dir ref_dir], yl, 'Color',[0.9 0.5 0.1], 'LineStyle','-.', 'LineWidth',1.8);
    h_mu  = plot([mean_ms mean_ms], yl, 'k--','LineWidth',1.8);
    ylim(yl);

    xlabel('|v_{\infty, ball}|  [m/s]');  ylabel('PDF');
    title(sprintf('Ballistic escape v_{\\infty} — %s (%s)', cname, region));
    legend([h_fb h_dir h_mu], ...
        {sprintf('Ranuschio v_{inf,fb} = %.1f', ref_fb), ...
         sprintf('Ranuschio v_{inf,dir} = %.1f', ref_dir), ...
         sprintf('mean = %.1f m/s', mean_ms)}, 'Location','best');
    hold off;
end
