function [result_fmincon, ultima_iter, info] = refinement_pre_flyby(x0, c_const, state_tcm, epoch_tcm, state_inj, epoch_inj, dsm1_epoch_ref)
% refinement_pre_flyby  Run the pre-flyby (injection + DSM1) fmincon refinement
% in the full ephemeris model.
%
% Inputs:
%   x0             - initial guess [inj(1:3) dsm1(4:6) dsm1_epoch_var(7)]
%   c_const        - constants struct
%   state_tcm      - state at TCM1, heliocentric J2000 [km, km/s]
%   epoch_tcm      - TCM1 epoch [ET s]
%   state_inj      - state at injection, heliocentric J2000 [km, km/s]
%   epoch_inj      - injection epoch [ET s]
%   dsm1_epoch_ref - reference DSM1 epoch [ET s]
%
% Outputs:
%   result_fmincon - fmincon output struct
%   ultima_iter    - optimized design vector
%   info           - convergence/feasibility info (exitflag, constrviolation, fval, success)

lb = [];
ub = [];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Nonlinear constraints (if needed)
nonlcon = @(x) NC_refinement_pre_flyby(x, c_const, state_tcm, epoch_tcm, state_inj, epoch_inj, dsm1_epoch_ref);

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...              % SQP is robust
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',100000, ...
    'OptimalityTolerance',1e-4, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1, ...
    'OutputFcn', @pre_flyby_output_fcn);

% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) OF_refinement_pre_flyby(x), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);


result_fmincon = output;
ultima_iter = x_opt;

% Convergence / feasibility info (used by run_refinement to flag the stage)
info.exitflag       = exitflag;
info.constrviolation = output.constrviolation;
info.fval           = fval;
% exitflag == -1 = stopped by our OutputFcn (feasible + stagnant): treat it as
% success if the solution is actually feasible.
info.success        = (exitflag > 0 || exitflag == -1) && (output.constrviolation <= options.ConstraintTolerance);


end


% =========================================================================
% Custom OutputFcn: stop when the solution is feasible AND stagnant.
% (Same mechanism used in optimization_ephe.)
% =========================================================================
function stop = pre_flyby_output_fcn(~, optimValues, state)
% pre_flyby_output_fcn  fmincon output function that halts the solver once the
% iterate is feasible and fval has stopped improving.
    stop = false;

    % --- tunable parameters ---
    feas_tol   = 1;      % threshold on constrviolation (aligned to ConstraintTolerance)
    stag_tol   = 0.1/1e+3;   % minimum fval change counted as "progress"
    stag_iters = 10;     % consecutive iterations without progress -> stop

    persistent fval_hist
    switch state
        case 'init'
            fval_hist = [];
        case 'iter'
            fval_hist(end+1) = optimValues.fval; %#ok<AGROW>
            is_feasible = optimValues.constrviolation < feas_tol;
            if numel(fval_hist) >= stag_iters
                window = fval_hist(end-stag_iters+1:end);
                is_stagnant = (max(window) - min(window)) < stag_tol;
            else
                is_stagnant = false;
            end
            if is_feasible && is_stagnant
                fprintf('\n[pre_flyby OutputFcn] Stop: feasibile (%.2e) e stagnante da %d iter.\n', ...
                        optimValues.constrviolation, stag_iters);
                stop = true;
            end
        case 'done'
            fval_hist = [];
    end
end

