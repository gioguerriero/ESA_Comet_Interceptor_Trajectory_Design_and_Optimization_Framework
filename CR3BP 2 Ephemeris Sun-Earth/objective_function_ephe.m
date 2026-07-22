function dv = objective_function_ephe(x, epoch_comet, multiple_shooting)
% objective_function_ephe  Ephemeris-refinement cost: weighted sum of the
% injection, DSM1 and DSM2 delta-v magnitudes.
%
% Inputs:
%   x                 - design-variable vector (layout depends on multiple_shooting)
%   epoch_comet       - comet arrival epoch [ET s] (fixed)
%   multiple_shooting - 1 = multiple-shooting layout, 0 = single-shooting layout
%
% Outputs:
%   dv - scalar cost (weighted total maneuver delta-v)

    % tof = (x(4) + x(8) + (epoch_comet - x(9)*1e+8) / 86400)/1000; % tof in days
    tof = x(7) - x(15); % tof in secondi/1e+8
    tof = 0;

    if multiple_shooting
        dv = norm(x(1:3))  / 1e3 ...   % dv_inj
           + norm(x(4:6))  / 10  ...   % dv_dsm1
           + norm(x(11:13))/ 10;       % dv_dsm2
    else
        dv = norm(x(1:3))  / 1e3 ...   % dv_inj
           + norm(x(5:7))  / 10  ...   % dv_dsm1
           + norm(x(14:16))/ 10;       % dv_dsm2
    end

    dv = 1e-0 * dv;
    % dv = 0;
    
end