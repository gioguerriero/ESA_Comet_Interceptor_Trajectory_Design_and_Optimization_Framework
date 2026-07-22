function [c, ceq] = NC_refinement_global_post_flyby(x, c_const, epoch_tcm1, comet_pos, epoch_comet, t_tcm2, epoch_flyby, state_pre_tcm, h_min_moon, dsm2_epoch_ref)
% NC_refinement_global_post_flyby  Nonlinear constraints for the global post-flyby
% refinement: B-plane flyby altitude limit plus a smooth patch between the TCM2 arc
% and a back-propagated arc from the comet.
%
% Inputs:
%   x              - design variables (B-plane targets, DSM2, epoch, comet velocity, TCM2)
%   c_const        - constants struct
%   epoch_tcm1     - TCM1 epoch [ET s]
%   comet_pos      - target comet position [km]
%   epoch_comet    - comet rendezvous epoch [ET s]
%   t_tcm2         - time after flyby when TCM2 is applied [days]
%   epoch_flyby    - lunar flyby epoch [ET s]
%   state_pre_tcm  - state before TCM1, heliocentric J2000 [km, km/s]
%   h_min_moon     - minimum lunar flyby altitude [km]
%   dsm2_epoch_ref - reference DSM2 epoch [ET s]
%
% Outputs:
%   c   - inequality constraints (<= 0): maneuver spacing and flyby altitude
%   ceq - equality constraints (= 0): position/velocity patch match

c = [];
ceq = [];

% for the optimization with DSM2 and TCM2 but without a multiple-shooting state
%{
BT = x(1)*1e+3;
BR = x(2)*1e+3;
tcm2 = x(3:5);
dsm2 = x(6:8);
epoch_dsm2 = x(9)*1e+8;
%}

% for the optimization without TCM2
%{
BT = x(1)*1e+3;
BR = x(2)*1e+3;
dsm2 = x(3:5);
epoch_dsm2 = x(6)*1e+8;
%}

% for the optimization without TCM2 but with a post-flyby state (a kind of multiple shooting)
BT = x(1)*5000;
BR = x(2)*5000;
dsm2 = x(3:5)'.*0.1;
dsm2_epoch_variation = x(6);
epoch_dsm2 = dsm2_epoch_variation*86400 + dsm2_epoch_ref;
comet_arrival_velocity = x(7:9).*10;
tcm2 = x(10:12)'.*0.1;

c = [c; (epoch_dsm2 - epoch_comet + 20*86400)/3600; (epoch_flyby + t_tcm2*86400 + 20*86400 - epoch_dsm2)/3600]; % keep DSM2 well separated from the other maneuvers
% c = [c; (epoch_dsm2 - epoch_comet + 20*86400)/1e+7; (epoch_flyby + 20*86400 - epoch_dsm2)/1e+7];

%% Common code amongst methods

options.verbose = 0;
[dv_tcm1, dv_mag, info] = bplane_tcm(state_pre_tcm, epoch_tcm1, epoch_flyby, ...
                                     'MOON', BT, BR, ...
                                     c_const.muMoon, c_const,options);
h_flyby = info.h_flyby;

c = [c; (-h_flyby + h_min_moon)/0.5];

S_post_tcm1 = state_pre_tcm;
S_post_tcm1(4:6) = S_post_tcm1(4:6) + dv_tcm1;

%% Constraints when both TCM2 and DSM2 are present
%{
opt_traj = odeset('AbsTol',1e-12,'RelTol',1e-12);
time_tcm1_tcm2 = epoch_flyby - epoch_tcm1 + t_tcm2*86400;
tspan = [0, time_tcm1_tcm2];
[~, S_tcm1_tcm2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm1, c_const), tspan, S_post_tcm1, opt_traj);

S_pre_tcm2 = S_tcm1_tcm2(end,:);
S_post_tcm2 = S_pre_tcm2';
S_post_tcm2(4:6) = S_post_tcm2(4:6) + tcm2;


epoch_tcm2 = epoch_flyby + t_tcm2*86400;
tspan2 = [0 epoch_dsm2 - epoch_tcm2];
[~, S_tcm2dsm] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_tcm2,c_const), tspan2, S_post_tcm2, opt_traj);
S_pre_dsm2 = S_tcm2dsm(end,:);
S_post_dsm2 = S_pre_dsm2';
S_post_dsm2(4:6) = S_post_dsm2(4:6) + dsm2;

tspan3 = [0 epoch_comet - epoch_dsm2];
[~, S_dsm2comet] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_dsm2,c_const), tspan3, S_post_dsm2, opt_traj);
S_comet = S_dsm2comet(end,:)';

ceq = [ceq; (S_comet(1:3) - comet_pos')./1e+6];
%}


%% Constraints when only DSM2 is present
%{

opt_traj = odeset('AbsTol',1e-12,'RelTol',1e-12);
time_tcm1_dsm2 = epoch_dsm2 - epoch_tcm1;
tspan = [0, time_tcm1_dsm2];
[~, S_tcm1_dsm2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm1, c_const), tspan, S_post_tcm1, opt_traj);

S_pre_dsm2 = S_tcm1_dsm2(end,:);
S_post_dsm2 = S_pre_dsm2';
S_post_dsm2(4:6) = S_post_dsm2(4:6) + dsm2;

tspan2 = [0 epoch_comet - epoch_dsm2];
[~, S_dsm2comet] = ode45(@(t,S) NBODY_J2000_full_ephe(t,S,epoch_dsm2,c_const), tspan2, S_post_dsm2, opt_traj);
S_comet = S_dsm2comet(end,:)';

ceq = [ceq; (S_comet(1:3) - comet_pos')./1e+6];
%}

%% Constraints when a control state is also used

opt_traj = odeset('AbsTol',1e-12,'RelTol',1e-12);
time_tcm1_tcm2 = epoch_flyby - epoch_tcm1 + t_tcm2*86400;
tspan = [0, time_tcm1_tcm2];
[~, S_tcm1_tcm2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm1, c_const), tspan, S_post_tcm1, opt_traj);

% Apply the TCM2 maneuver
S_pre_tcm2 = S_tcm1_tcm2(end,:);
S_post_tcm2 = S_pre_tcm2;
S_post_tcm2(4:6) = S_post_tcm2(4:6) + tcm2;

% Propagate 5 days forward to the patching point
days_patching = 5;
epoch_tcm2 = epoch_flyby + t_tcm2*86400;
[~, S_patch_pre] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm2, c_const), [0 days_patching*86400], S_post_tcm2, opt_traj);
S_patch_pre = S_patch_pre(end,:);

% Now back-propagate from the comet
state_at_comet = [comet_pos(:); comet_arrival_velocity];
[~, S_post_dsm2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_comet, c_const), [0 -(epoch_comet-epoch_dsm2)], state_at_comet, opt_traj);
S_post_dsm2 = S_post_dsm2(end,:);

S_pre_dsm = S_post_dsm2;
S_pre_dsm(4:6) = S_pre_dsm(4:6) - dsm2;

dt_back = (epoch_tcm2 + days_patching*86400) - epoch_dsm2;
[~, S_patch_post] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_dsm2, c_const), [0 dt_back], S_pre_dsm, opt_traj);
S_patch_post = S_patch_post(end,:);

ceq = [ceq; (S_patch_post(1:3)' - S_patch_pre(1:3)')./1]; % ~1 km tolerance
ceq = [ceq; (S_patch_post(4:6)' - S_patch_pre(4:6)')./1e-4]; % ~0.1 m/s tolerance


end
