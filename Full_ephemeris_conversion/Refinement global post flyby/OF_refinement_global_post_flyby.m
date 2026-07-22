function dv = OF_refinement_global_post_flyby(x, c_const, epoch_tcm1, comet_pos, epoch_comet, t_tcm2, epoch_flyby, state_pre_tcm, h_min_moon)
% OF_refinement_global_post_flyby  Objective for the global (MBH) post-flyby
% refinement: total delta-v of the TCM1 (via B-plane), TCM2 and DSM2 maneuvers.
%
% Inputs:
%   x             - design variables (B-plane targets, DSM2, epoch, comet velocity, TCM2)
%   c_const       - constants struct
%   epoch_tcm1    - TCM1 epoch [ET s]
%   comet_pos     - target comet position [km]
%   epoch_comet   - comet rendezvous epoch [ET s]
%   t_tcm2        - time after flyby when TCM2 is applied [days]
%   epoch_flyby   - lunar flyby epoch [ET s]
%   state_pre_tcm - state before TCM1, heliocentric J2000 [km, km/s]
%   h_min_moon    - minimum lunar flyby altitude [km] (used by the constraints)
%
% Outputs:
%   dv - total delta-v magnitude to minimize [km/s]

%% Optimization with TCM2 and DSM2
    % Minimize the total delta-v of the two maneuvers (TOF stays constant)
    % dv = norm(x(3:5)) + norm(x(6:8));

%% Optimization with TCM2 only
    % dv = norm(x(3:5)); % cost function quando c'è solo DSM2
    % dv = 0;

%% Optimization with TCM2, DSM2 and a control point
BT = x(1)*5000;
BR = x(2)*5000;
dsm2 = x(3:5).*0.1;
dsm2_epoch_variation = x(6);
comet_arrival_velocity = x(7:9).*10;
dv_tcm2 = x(10:12).*0.1;

options.verbose = 0;
[dv_tcm1, dv_mag, info] = bplane_tcm(state_pre_tcm, epoch_tcm1, epoch_flyby, ...
                                     'MOON', BT, BR, ...
                                     c_const.muMoon, c_const,options);

% S_post_tcm1 = state_pre_tcm;
% S_post_tcm1(4:6) = S_post_tcm1(4:6) + dv_tcm1;
% 
% opt_traj = odeset('AbsTol',1e-12,'RelTol',1e-12);
% time_tcm1_tcm2 = epoch_flyby - epoch_tcm1 + t_tcm2*86400;
% tspan = [0, time_tcm1_tcm2];
% [~, S_tcm1_tcm2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm1, c_const), tspan, S_post_tcm1, opt_traj);
% 
% S_pre_tcm2 = S_tcm1_tcm2(end,:);
% 
% dv_tcm2 = state_post_tcm2(4:6) - S_pre_tcm2(4:6)';

dv = norm(dv_tcm1) + norm(dv_tcm2) + norm(dsm2);

% dv = 1e-2 * dv;


end