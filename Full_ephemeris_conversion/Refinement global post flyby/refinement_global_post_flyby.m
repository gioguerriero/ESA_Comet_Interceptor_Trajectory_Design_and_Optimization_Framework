function [result_fmincon, ultima_iter, info] = refinement_global_post_flyby(x0, c_const, epoch_tcm1, comet_pos, comet_epoch, t_tcm2, epoch_flyby, state_pre_tcm, h_min_moon, dsm2_epoch_ref)
% refinement_global_post_flyby  Run the global (B-plane targeted) post-flyby
% fmincon refinement in the full ephemeris model.
%
% Inputs:
%   x0             - initial guess (B-plane targets, DSM2, epoch, comet velocity, TCM2)
%   c_const        - constants struct
%   epoch_tcm1     - TCM1 epoch [ET s]
%   comet_pos      - target comet position [km]
%   comet_epoch    - comet rendezvous epoch [ET s]
%   t_tcm2         - time after flyby when TCM2 is applied [days]
%   epoch_flyby    - lunar flyby epoch [ET s]
%   state_pre_tcm  - state before TCM1, heliocentric J2000 [km, km/s]
%   h_min_moon     - minimum lunar flyby altitude [km]
%   dsm2_epoch_ref - reference DSM2 epoch [ET s]
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
nonlcon = @(x) NC_refinement_global_post_flyby(x, c_const, epoch_tcm1, comet_pos, comet_epoch, t_tcm2, epoch_flyby, state_pre_tcm, h_min_moon, dsm2_epoch_ref);

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ... 
    'Algorithm','sqp', ...              % SQP is robust sqp
    'MaxIterations',300, ...
    'MaxFunctionEvaluations',100000, ...
    'OptimalityTolerance',1e-3, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1,...
    'FiniteDifferenceType','central', ...
    'OutputFcn', @post_flyby_output_fcn);

% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) OF_refinement_global_post_flyby(x, c_const, epoch_tcm1, comet_pos, comet_epoch, t_tcm2, epoch_flyby, state_pre_tcm, h_min_moon), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);


result_fmincon = output;
ultima_iter = x_opt;

% Convergence / feasibility info (used by run_refinement to flag the stage)
info.exitflag        = exitflag;
info.constrviolation = output.constrviolation;
info.fval            = fval;
% exitflag == -1 = stopped by our OutputFcn (feasible + stagnant): treat it as
% success if the solution is actually feasible.
info.success         = (exitflag > 0 || exitflag == -1) && (output.constrviolation <= options.ConstraintTolerance);


end


% =========================================================================
% Custom OutputFcn: stop when the solution is feasible AND stagnant.
% (Same mechanism used in optimization_ephe.)
% =========================================================================
function stop = post_flyby_output_fcn(~, optimValues, state)
% post_flyby_output_fcn  fmincon output function that halts the solver once the
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
                fprintf('\n[post_flyby OutputFcn] Stop: feasibile (%.2e) e stagnante da %d iter.\n', ...
                        optimValues.constrviolation, stag_iters);
                stop = true;
            end
        case 'done'
            fval_hist = [];
    end
end

