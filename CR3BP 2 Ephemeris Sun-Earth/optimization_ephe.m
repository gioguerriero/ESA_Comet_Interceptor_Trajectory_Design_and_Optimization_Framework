function [result_fmincon, ultima_iter] = optimization_ephe(x0, S_halo, max_dv, target_position, epoch_comet_flyby, hp, c, max_dv_inj, multiple_shooting, k_vec, f_nodes, d_max, dt_max, optimality_tolerance, constraint_tolerance, min_days_between)
% optimization_ephe  Ephemeris-refinement fmincon driver (multiple-shooting).
% Sets up event-ordering linear constraints and control-point position bounds,
% then minimizes the total maneuver delta-v under the closure constraints.
%
% Inputs:
%   x0                   - initial guess (base + node states)
%   S_halo               - halo orbit states
%   max_dv               - total delta-v budget [m/s]
%   target_position      - target comet position [km]
%   epoch_comet_flyby    - comet arrival epoch [ET s]
%   hp                   - minimum lunar flyby altitude [km]
%   c                    - constants struct
%   max_dv_inj           - maximum injection delta-v [m/s]
%   multiple_shooting    - 1 = multiple shooting, 0 = single shooting
%   k_vec                - [k1 k2 k3 k4] nodes per segment
%   f_nodes              - fixed node time fractions (from build_initial_guess_MS)
%   d_max                - position bound on each control point [synodic non-dim]
%   dt_max               - (unused here) epoch bound helper
%   optimality_tolerance - fmincon OptimalityTolerance
%   constraint_tolerance - fmincon ConstraintTolerance
%   min_days_between     - min event spacing [inj->DSM1, DSM1->flyby, flyby->DSM2, DSM2->comet] [days]
%
% Outputs:
%   result_fmincon - fmincon output struct
%   ultima_iter    - optimized design vector

% Multiple shooting: epoch constraints (chronological order for segments 1-2,
% reversed for segments 3-4) and +/-d_max bounds on each control point's
% position components.

if multiple_shooting

    k1 = k_vec(1);  k2 = k_vec(2);
    k3 = k_vec(3);  k4 = k_vec(4);

    n_base = 17;                              % DV base (dv + 4 event epochs)
    off1 = n_base;
    off2 = off1 + 6*k1;
    off3 = off2 + 6*k2;
    off4 = off3 + 6*k3;
    n_tot = off4 + 6*k4;                      % total length of x (node states only)

    % Indices of the node STATES (node epochs are NOT variables: they are
    % reconstructed from f_nodes and the 4 event epochs).
    idx_S1 = (off1+1):(off1+6*k1);
    idx_S2 = (off2+1):(off2+6*k2);
    idx_S3 = (off3+1):(off3+6*k3);
    idx_S4 = (off4+1):(off4+6*k4);

    epoch_idx = [7 15 16 17];   % event-epoch indices: flyby, dep, dsm1, dsm2

    %% ===== LINEAR INEQUALITY CONSTRAINTS (event ordering) =======
    % The ONLY variable epochs are the 4 physical events:
    %   dep = x(15), dsm1 = x(16), flyby = x(7), dsm2 = x(17).
    % In "days before comet arrival" (x = (epoch_comet - epoch)/86400) the
    % chronological chain is monotonically DECREASING:
    %   x_dep > x_dsm1 > x_flyby > x_dsm2 > 0 (comet).
    % Node epochs are reconstructed by fixed fraction (T = Ts+f*(Te-Ts)) and are
    % therefore auto-ordered: no constraints on the nodes are needed.
    % For each pair (E=earlier, L=later): x_E - x_L >= dt_min
    %   <=>  x_L - x_E <= -dt_min   (row(L)=+1, row(E)=-1).
    % Minimum spacing [days] per consecutive event pair (user-configurable):
    %   b(1): inj->dsm1, b(2): dsm1->flyby, b(3): flyby->dsm2
    A = zeros(3, n_tot);
    A(1,16) =  1;  A(1,15) = -1;         % dsm1  - dep    (dep  > dsm1)
    A(2, 7) =  1;  A(2,16) = -1;         % flyby - dsm1   (dsm1 > flyby)
    A(3,17) =  1;  A(3, 7) = -1;         % dsm2  - flyby  (flyby> dsm2)
    b = -min_days_between(1:3).';        % [-d_inj_dsm1; -d_dsm1_flyby; -d_flyby_dsm2]
    Aeq = [];
    beq = [];

    %% ===== UPPER / LOWER BOUNDS =====================================
    lb = -inf(n_tot, 1);
    ub =  inf(n_tot, 1);

    % Event epochs (days before comet): +/- 4 days from the guess
    %   x(7)=flyby, x(15)=dep, x(16)=dsm1, x(17)=dsm2
    for ie = epoch_idx
        lb(ie) = x0(ie) - 4;
        ub(ie) = x0(ie) + 4;
    end

    % t_halo (x(14), days, NOT shifted): +/- 5 days from the guess
    lb(14) = x0(14) - 5;
    ub(14) = x0(14) + 5;

    % dsm2 must occur BEFORE comet arrival: in "days before comet" the comet is
    % x=0, so we require x_dsm2 >= d_dsm2_comet (DSM2->comet spacing).
    lb(17) = max(lb(17), min_days_between(4));

    % Bounds on the 3 position components of each control point
    % (synodic non-dim, +/-d_max). Velocities left free.
    pos_comp = [1 2 3];
    blocks   = {idx_S1, idx_S2, idx_S3, idx_S4};
    ks       = [k1 k2 k3 k4];
    for bi = 1:4
        if ks(bi) == 0, continue; end
        S_idx_mat = reshape(blocks{bi}, 6, ks(bi));   % 6×k indices of the block
        for j = 1:ks(bi)
            for p = pos_comp
                idx = S_idx_mat(p, j);
                lb(idx) = x0(idx) - d_max;
                ub(idx) = x0(idx) + d_max;
            end
        end
    end

    nonlcon = @(x) nonlinear_constraints_ephe_multiple(x, S_halo, max_dv, target_position, epoch_comet_flyby, hp, c, max_dv_inj, k_vec, f_nodes);

end   % closes the "if multiple_shooting" block

% else
%     % Nonlinear constraints (if needed)
%     % Bounds
%     %   x layout (17 variabili):
%     %     x(1)     – tof_wait_d        [days]   guess ± 10 giorni
%     %     x(2:4)   – dv_inj            [km/s]   no bound (constraint sul modulo)
%     %     x(5)     – tof_inj2dsm1_d   [days]   guess ± 20 giorni
%     %     x(6:8)   – dv_dsm1           [km/s]   no bound
%     %     x(9)     – tof_dsm12flyby_d [days]   guess ± 20 giorni
%     %     x(10)    – epoch_flyby       [ET s / 1e+7]   guess ± 5 giorni
%     %     x(11:13) – vel_comet_arr     [km/s]   no bound
%     %     x(14)    – tof_dsm22comet_d [days]   guess ± 90 giorni
%     %     x(15:17) – dv_dsm2           [km/s]   no bound
% 
%     % Linear constraints (none)
%     A = [];
%     b = [];
%     Aeq = [];
%     beq = [];
% 
%     lb = -inf(17,1);
%     ub =  inf(17,1);
% 
%     % tof_wait_d  [days]
%     % lb(1)  = x0(1) - 8;
%     % ub(1)  = x0(1) + 8;
% 
%     % tof_inj2dsm1_d  [days]
%     lb(4)  = max(x0(4) - 10, 5);
%     ub(4)  = x0(4) + 10;
% 
%     % tof_dsm12flyby_d  [days]
%     lb(8)  = max(x0(8) - 30,10);
%     ub(8)  = x0(8) + 30;
% 
%     % epoch_flyby  [ET s]  ± 1 giorni
%     lb(9) = x0(9) - 2*86400/1e+8;
%     ub(9) = x0(9) + 2*86400/1e+8;
% 
%     % tof_dsm22comet_d  [days]
%     lb(13) = max(x0(13) - 100/100,10/100);
%     ub(13) = x0(13) + 100/100;
% 
%     % Spostamento sulla Halo [days]
%     lb(17) = x0(17) - 5;
%     ub(17) = x0(17) + 5;
% 
%     nonlcon = @(x) nonlinear_constraints_ephe(x, S_halo, max_dv, target_position, epoch_comet_flyby, hp, c, max_dv_inj);
% end

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...              % SQP is robust
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',100000, ...
    'OptimalityTolerance',optimality_tolerance, ...
    'StepTolerance',1e-8, ...
    'ConstraintTolerance',constraint_tolerance,...
    'FiniteDifferenceType','central', ...
    'OutputFcn', @ms_output_fcn);
% 'FiniteDifferenceStepSize',0.5e-4
% 'FunctionTolerance', 1e-4,...

% Call fmincon
[x_opt, fval, exitflag, output, lambda] = fmincon( ...
    @(x) objective_function_ephe(x, epoch_comet_flyby, multiple_shooting), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);

% Output
% result_fmincon.x_opt = x_opt;
% result_fmincon.fval  = fval;
% result_fmincon.exitflag = exitflag;
result_fmincon = output;
ultima_iter = x_opt;


end


% =========================================================================
% Custom OutputFcn: stop when the solution is feasible AND stagnant.
% Avoids spinning to reduce first-order optimality (which, with finite
% differences over adaptive ode45, does not converge organically).
% =========================================================================
function stop = ms_output_fcn(~, optimValues, state)
% ms_output_fcn  fmincon output function that halts once the iterate is
% feasible and fval has stopped improving.
    stop = false;

    % --- parameters (edit here to change the behavior) ---
    feas_tol   = 1e-3;   % constrviolation threshold to consider a point "feasible"
    stag_tol   = 0.1/1e+3;   % minimum fval change counted as "progress" -> in [m/s]
    stag_iters = 10;     % consecutive iterations without progress -> stop

    % TO BE REMOVED
    stag_tol   = 20/1e+3;   % minimum fval change counted as "progress" -> in [m/s]
    stag_iters = 1;     % consecutive iterations without progress -> stop

    persistent fval_hist

    switch state
        case 'init'
            fval_hist = [];

        case 'iter'
            fval_hist(end+1) = optimValues.fval; %#ok<AGROW>

            is_feasible = optimValues.constrviolation < feas_tol;

            n = numel(fval_hist);
            if n >= stag_iters
                window = fval_hist(end-stag_iters+1:end);
                is_stagnant = (max(window) - min(window)) < stag_tol;
            else
                is_stagnant = false;
            end

            if is_feasible && is_stagnant
                fprintf('\n[OutputFcn] Stop anticipato: feasibile (%.2e) e stagnante da %d iter.\n', ...
                        optimValues.constrviolation, stag_iters);
                stop = true;
            end

        case 'done'
            fval_hist = [];
    end
end
