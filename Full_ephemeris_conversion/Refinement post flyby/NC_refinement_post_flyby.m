function [c, ceq] = NC_refinement_post_flyby(x, c_const, state_tcm, dv_tcm, epoch_tcm, comet_pos, epoch_comet, t_tcm2, epoch_flyby)
% NC_refinement_post_flyby  Nonlinear constraints for the post-flyby refinement:
% propagate TCM1 -> TCM2 -> DSM2 -> comet and match the comet position.
%
% Inputs:
%   x           - design variables [tcm2(1:3) dsm2(4:6) dsm2_epoch/1e8 (7)]
%   c_const     - constants struct
%   state_tcm   - state at TCM1, heliocentric J2000 [km, km/s]
%   dv_tcm      - TCM1 delta-v applied first [km/s]
%   epoch_tcm   - TCM1 epoch [ET s]
%   comet_pos   - target comet position [km]
%   epoch_comet - comet rendezvous epoch [ET s]
%   t_tcm2      - time after flyby when TCM2 is applied [days]
%   epoch_flyby - lunar flyby epoch [ET s]
%
% Outputs:
%   c   - inequality constraints (<= 0), maneuver spacing
%   ceq - equality constraints (= 0), comet position match

c = [];
ceq = [];

tcm2 = x(1:3);
dsm2 = x(4:6);
epoch_dsm = x(7)*1e+8;

c = [c; (epoch_dsm - epoch_comet + 20*86400)/1e+7; (epoch_tcm + 10*86400 - epoch_dsm)/1e+7]; % keep DSM2 well separated from the other maneuvers
% c = [c; norm(tcm2) - 0.1]; % limite al modulo della TCM [km/s]

% forward propagation
opt_traj = odeset('AbsTol',1e-13,'RelTol',1e-13);

state_tcm(4:6) = state_tcm(4:6) + dv_tcm;
epoch_tcm2 = epoch_flyby+t_tcm2*86400;
tspan1 = [0 epoch_tcm2 - epoch_tcm];


[~, S_tcm2tcm] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_tcm,c_const), tspan1, state_tcm, opt_traj);
S_pre_tcm2 = S_tcm2tcm(end,:);
S_post_tcm2 = S_pre_tcm2';
S_post_tcm2(4:6) = S_post_tcm2(4:6) + tcm2;


tspan2 = [0 epoch_dsm - epoch_tcm2];
[~, S_tcm2dsm] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_tcm2,c_const), tspan2, S_post_tcm2, opt_traj);
S_pre_dsm2 = S_tcm2dsm(end,:);
S_post_dsm2 = S_pre_dsm2';
S_post_dsm2(4:6) = S_post_dsm2(4:6) + dsm2;

tspan3 = [0 epoch_comet - epoch_dsm];
[~, S_dsm2comet] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_dsm,c_const), tspan3, S_post_dsm2, opt_traj);
S_comet = S_dsm2comet(end,:)';

ceq = [ceq; (S_comet(1:3) - comet_pos')./1e+6];


end
