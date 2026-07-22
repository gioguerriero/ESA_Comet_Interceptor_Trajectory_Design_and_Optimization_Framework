function dv = OF_refinement_pre_flyby(x)
% OF_refinement_pre_flyby  Objective for the pre-flyby refinement: total
% magnitude of the two maneuvers (injection + DSM1). TOF stays constant.
%
% Inputs:
%   x  - design variables (uses x(1:3) and x(4:6) as the two delta-v vectors)
%
% Outputs:
%   dv - combined delta-v magnitude to minimize

    % Minimize the total delta-v of the two maneuvers (TOF stays constant)
    dv = norm(x(1:3)) + norm(x(4:6));

    % dv = 1e-2 * dv;
    % dv = 0;

end