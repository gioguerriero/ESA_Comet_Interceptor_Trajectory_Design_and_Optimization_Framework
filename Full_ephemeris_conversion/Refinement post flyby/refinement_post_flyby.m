function [result_fmincon, ultima_iter, info] = refinement_post_flyby(x0, c_const, state_tcm, dv_tcm, epoch_tcm, comet_pos, comet_epoch, t_tcm2, epoch_flyby)
% refinement_post_flyby  Run the post-flyby (TCM2 + DSM2) fmincon refinement
% in the full ephemeris model, targeting the comet position.
%
% Inputs:
%   x0          - initial guess [tcm2(1:3) dsm2(4:6) dsm2_epoch(7)]
%   c_const     - constants struct
%   state_tcm   - state at TCM1, heliocentric J2000 [km, km/s]
%   dv_tcm      - TCM1 delta-v already applied [km/s]
%   epoch_tcm   - TCM1 epoch [ET s]
%   comet_pos   - target comet position [km]
%   comet_epoch - comet rendezvous epoch [ET s]
%   t_tcm2      - time after flyby when TCM2 is applied [days]
%   epoch_flyby - lunar flyby epoch [ET s]
%
% Outputs:
%   result_fmincon - optimized design vector
%   ultima_iter    - optimized design vector (same as result_fmincon)
%   info           - convergence/feasibility info (exitflag, constrviolation, fval, success)

lb = [];
ub = [];

% Linear constraints (none)
A = [];
b = [];
Aeq = [];
beq = [];

% Nonlinear constraints (if needed)
nonlcon = @(x) NC_refinement_post_flyby(x, c_const, state_tcm, dv_tcm, epoch_tcm, comet_pos, comet_epoch, t_tcm2, epoch_flyby);

% fmincon Options
options = optimoptions('fmincon', ...
    'Display','iter', ...
    'Algorithm','sqp', ...              % SQP is robust
    'MaxIterations',200, ...
    'MaxFunctionEvaluations',20000, ...
    'OptimalityTolerance',1e-4, ...
    'StepTolerance',1e-10, ...
    'ConstraintTolerance',1e-6);

% Call fmincon
[x_opt, fval, exitflag, output] = fmincon( ...
    @(x) OF_refinement_post_flyby(x), ...
    x0, A, b, Aeq, beq, lb, ub, nonlcon, options);


result_fmincon = x_opt;
ultima_iter = x_opt;

% Convergence / feasibility info (used by run_refinement to flag the stage)
info.exitflag        = exitflag;
info.constrviolation = output.constrviolation;
info.fval            = fval;
info.success         = (exitflag > 0) && (output.constrviolation <= options.ConstraintTolerance);


end

