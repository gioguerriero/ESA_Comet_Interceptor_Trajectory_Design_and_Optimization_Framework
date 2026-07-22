function [c, ceq] = NC_refinement_post_flyby(x, c_const, state_tcm, dv_tcm, epoch_tcm, comet_pos, epoch_comet, t_tcm2, epoch_flyby)
% Nonlinear constraint function for post-flyby trajectory refinement.
%
% Decision variables x (7 elements):
%   x(1:3)  – TCM2 delta-v [km/s]   (post-flyby cleanup maneuver)
%   x(4:6)  – DSM2 delta-v [km/s]   (deep-space maneuver toward comet)
%   x(7)    – epoch of DSM2 [s/1e8] (absolute SPICE ET, scaled)
%
% Fixed inputs:
%   state_tcm   – heliocentric J2000 state at TCM1 epoch [km, km/s]
%   dv_tcm      – TCM1 delta-v already computed [km/s]  (applied first)
%   epoch_tcm   – SPICE ET of TCM1 [s]
%   epoch_flyby – SPICE ET of lunar flyby [s]
%   t_tcm2      – time after flyby when TCM2 is applied [days]
%   comet_pos   – target comet position at epoch_comet [km, 1x3 or 3x1]
%   epoch_comet – SPICE ET of comet rendezvous [s]
%
% Constraints imposed:
%   c   = []   (no inequality constraints)
%   ceq = 0    position match at comet epoch (3 components, scaled by 1e6 km)

c   = [];
ceq = [];

tcm2       = x(1:3);
dsm2       = x(4:6);
epoch_dsm  = x(7) * 1e+8;   % unscale to SPICE ET [s]

opt_traj = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);

% Apply TCM1 to the state at epoch_tcm
state_tcm(4:6) = state_tcm(4:6) + dv_tcm;

% Absolute epoch of TCM2 (fixed offset after flyby)
epoch_tcm2 = epoch_flyby + t_tcm2 * 86400;

% Segment 1: TCM1 epoch -> TCM2 epoch
tspan1 = [0, epoch_tcm2 - epoch_tcm];
[~, S_tcm2tcm] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_tcm, c_const), ...
                        tspan1, state_tcm, opt_traj);
S_pre_tcm2  = S_tcm2tcm(end, :);
S_post_tcm2 = S_pre_tcm2';
S_post_tcm2(4:6) = S_post_tcm2(4:6) + tcm2;   % apply TCM2

% Segment 2: TCM2 epoch -> DSM2 epoch
tspan2 = [0, epoch_dsm - epoch_tcm2];
[~, S_tcm2dsm] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_tcm2, c_const), ...
                        tspan2, S_post_tcm2, opt_traj);
S_pre_dsm2  = S_tcm2dsm(end, :);
S_post_dsm2 = S_pre_dsm2';
S_post_dsm2(4:6) = S_post_dsm2(4:6) + dsm2;   % apply DSM2

% Segment 3: DSM2 epoch -> comet rendezvous epoch
tspan3 = [0, epoch_comet - epoch_dsm];
[~, S_dsm2comet] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_dsm, c_const), ...
                          tspan3, S_post_dsm2, opt_traj);
S_comet = S_dsm2comet(end, :)';

% Position match constraint (scaled by 1e6 km for numerical conditioning)
ceq = [ceq; (S_comet(1:3) - comet_pos(:)) / 1e6];

end
