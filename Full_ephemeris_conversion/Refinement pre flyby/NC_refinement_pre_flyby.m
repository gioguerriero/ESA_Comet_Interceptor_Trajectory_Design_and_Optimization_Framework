function [c, ceq] = NC_refinement_pre_flyby(x, c_const, state_tcm, epoch_tcm, state_inj, epoch_inj, dsm1_epoch_ref)
% NC_refinement_pre_flyby  Nonlinear constraints for the pre-flyby refinement:
% enforce a smooth (position/velocity) patch between the injection and TCM1 arcs.
%
% Inputs:
%   x              - design variables [inj(1:3) dsm1(4:6) dsm1_epoch_var(7)]
%   c_const        - constants struct
%   state_tcm      - state at TCM1, heliocentric J2000 [km, km/s]
%   epoch_tcm      - TCM1 epoch [ET s]
%   state_inj      - state at injection, heliocentric J2000 [km, km/s]
%   epoch_inj      - injection epoch [ET s]
%   dsm1_epoch_ref - reference DSM1 epoch [ET s]
%
% Outputs:
%   c   - inequality constraints (<= 0), maneuver spacing
%   ceq - equality constraints (= 0), position/velocity match at the patch point

c = [];
ceq = [];

inj = x(1:3);
dsm1 = x(4:6);
epoch_dsm1_variation = x(7);
epoch_dsm1 = epoch_dsm1_variation*86400 + dsm1_epoch_ref;

c = [c; (epoch_dsm1 - epoch_tcm + 10*86400)/3600; (epoch_inj + 10*86400 - epoch_dsm1)/3600]; % keep DSM1 well separated from the other maneuvers -> 1h tolerance
% c = [c; norm(tcm2) - 0.1]; % limite al modulo della TCM [km/s]

% forward propagation
opt_traj = odeset('AbsTol',1e-12,'RelTol',1e-12);

state_inj(4:6) = state_inj(4:6) + inj';
tspan1 = [0 (epoch_dsm1 - epoch_inj)/2];

[~, S_inj2pre_patch] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_inj, c_const), tspan1, state_inj, opt_traj);
S_pre_patch = S_inj2pre_patch(end,:);


% now back-propagate from TCM1
tspan2 = [0  -(epoch_tcm - epoch_dsm1)];
[~, S_tcm2dsm] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_tcm, c_const), tspan2, state_tcm, opt_traj);
S_post_dsm = S_tcm2dsm(end,:);

S_pre_dsm = S_post_dsm;
S_pre_dsm(4:6) = S_pre_dsm(4:6) - dsm1';

tspan3 = [0  -(epoch_dsm1 - epoch_inj)/2];
[~, S_dsm2post_patch] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_dsm1, c_const), tspan3, S_pre_dsm, opt_traj);
S_post_patch = S_dsm2post_patch(end,:);



% Match position and velocity at the patch point for a smooth trajectory
ceq = [ceq; (S_pre_patch(1:3) - S_post_patch(1:3))'./1; (S_pre_patch(4:6) - S_post_patch(4:6))'./1e-4]; % ~100 km and 10 m/s tolerance
% ceq = [ceq; (S_pre_tcm(1:3)' - state_tcm(1:3))./100];

end
