function [c, ceq] = nonlinear_constraints_ephe(x, S_halo, max_dv, target_position, epoch_comet_flyby, hp, c_const, max_dv_inj)
% nonlinear_constraints_ephe  Constraints for the direct ephemeris optimization
% (single-shooting variant).
%
% Inputs:
%   x                 - 20-element design-variable vector
%   S_halo            - halo state [km; km/s] (fixed parameter)
%   max_dv            - maximum delta-v budget [m/s]
%   target_position   - comet position at arrival [km]
%   epoch_comet_flyby - comet arrival epoch [ET s]
%   hp                - lunar flyby periapsis altitude [km]
%   c_const           - constants struct
%   max_dv_inj        - maximum injection delta-v [m/s]
%
% Outputs:
%   c   - inequality constraints (<= 0)
%   ceq - equality constraints (= 0)
%
% x layout (20 variables):
%   x(1)     - tof_wait_d        [days]  halo waiting time
%   x(2:4)   - dv_inj            [km/s]  injection delta-v
%   x(5)     - tof_inj2dsm1_d   [days]  injection -> DSM1
%   x(6:8)   - dv_dsm1           [km/s]  DSM1 delta-v
%   x(9)     - tof_dsm12flyby_d [days]  DSM1 -> flyby
%   x(10)    - epoch_flyby       [ET s]  lunar flyby epoch
%   x(11:16) - S_comet_arr       [km; km/s]  comet arrival state
%   x(17)    - tof_dsm22comet_d [days]  DSM2 -> comet
%   x(18:20) - dv_dsm2           [km/s]  DSM2 delta-v
%
% Forward-propagated design variables: halo waiting tof, injection dv,
%   INJ->DSM1 tof, DSM1 dv, DSM1->flyby tof.
% Backward-propagated design variables: flyby epoch, comet arrival velocity,
%   DSM2->comet tof, DSM2 dv.

c = [];
ceq = [];

%% ===== UNPACK DESIGN VARIABLES =======================================

% --- Forward trajectory (Halo -> flyby) -------------------------------
tof_onHalo_ad     = x(17)*(24*3600)/c_const.Tstar;   % shift along the Halo [days]
dv_inj           = x(1:3)'./1e+3;         % injection delta-v [km/s]
tof_inj2dsm1_d   = x(4);            % TOF injection -> DSM1 [days]
dv_dsm1          = x(5:7)'./10;         % DSM1 delta-v [km/s]
tof_dsm12flyby_d = x(8);            % TOF DSM1 -> flyby [days]

% --- Backward trajectory (comet -> flyby) ----------------------------
epoch_flyby      = x(9)*1e+8;           % lunar flyby epoch [ET s]
S_comet_arr_vel  = x(10:12)';       % comet arrival state [km; km/s]
tof_dsm22comet_d = x(13)*100;           % TOF DSM2 -> comet [days]
dv_dsm2          = x(14:16)'./10;       % DSM2 delta-v [km/s]

% --- Convert TOFs to seconds (for the integrations) -----------------
% tof_onHalo_s       = tof_onHalo_ad   * 86400;
tof_inj2dsm1_s   = tof_inj2dsm1_d   * 86400;
tof_dsm12flyby_s = tof_dsm12flyby_d * 86400;
tof_dsm22comet_s = tof_dsm22comet_d * 86400;

tof_flyby2dsm2_s = epoch_comet_flyby - tof_dsm22comet_s - epoch_flyby;

initial_epoch = epoch_flyby - tof_dsm12flyby_s - tof_inj2dsm1_s;

% S_halo is in synodic non-dim; convert it at the current initial_epoch
S_halo_J2000 = synodic2sun_J2000(S_halo(:)', initial_epoch, c_const.mu);
S_halo_J2000 = S_halo_J2000(:);   % always a 6×1 column

moon_state = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN')';


%% Forward propagation

opt_traj = odeset('AbsTol',1e-13,'RelTol',1e-13);

% % Shift along the Halo -> posso partire tot giorni prima o tot giorni dopo, può decidere l'ottimizzatore
% if abs(tof_onHalo_s) > 3600
%     tspan = [0 tof_onHalo_s];
%     [~, S_waiting] = ode45(@(t,S) NBODY_J2000(t,S,initial_epoch,c_const), tspan, S_halo, opt_traj);
% else
%     S_waiting = S_halo';
% end
% 
% % Injection maneuver -> vincolo sul modulo massimo
% c = [c; norm(dv_inj) - max_dv_inj*1e-3];
% S_injection = S_waiting(end,:);
% S_injection(4:6) = S_injection(4:6) + dv_inj;
% 
% % Injection to DSM1
% tspan = [0 tof_inj2dsm1_s];
% [~, S_preDSM1] = ode45(@(t,S) NBODY_J2000(t,S,initial_epoch,c_const), tspan, S_injection, opt_traj);
% 
% % DSM1 maneuver
% S_postDSM1 = S_preDSM1(end,:);
% S_postDSM1(4:6) = S_postDSM1(4:6) + dv_dsm1;
% 
% % DSM1 to moon
% tspan = [0 tof_dsm12flyby_s];
% [~, S_preFlyby] = ode45(@(t,S) NBODY_J2000(t,S,initial_epoch+tof_inj2dsm1_s,c_const), tspan, S_postDSM1, opt_traj);

% % Shift along the Halo
% if abs(tof_onHalo_ad) > (5*24*3600)/c_const.Tstar
%     tspan = [0 tof_onHalo_s];
%     [~, S_waiting] = ode45(@(t,S) NBODY_J2000(t,S,initial_epoch,c_const), tspan, S_halo_J2000, opt_traj);
% else
%     S_waiting = S_halo_J2000(:)';
% end

% Shift along the Halo
if abs(tof_onHalo_ad) > (3600)/c_const.Tstar
    tspan = [0 tof_onHalo_ad];
    [~, S_partenza] = ode45(@(t,S) CR3BP(t,S,c_const.mu), tspan, S_halo(:), opt_traj);
else
    S_partenza = S_halo(:)';
end

S_halo_J2000 = synodic2sun_J2000(S_partenza(end,:), initial_epoch, c_const.mu);
S_halo_J2000 = S_halo_J2000(:);   % always a 6×1 column

% Injection
% S_injection = S_waiting(end,:);
S_injection = S_halo_J2000';
S_injection(4:6) = S_injection(4:6) + dv_inj;
c = [c; norm(dv_inj)*1e+3 - max_dv_inj];

% Injection -> DSM1 (starts at initial_epoch)
tspan = [0 tof_inj2dsm1_s];
[~, S_preDSM1] = ode45(@(t,S) NBODY_J2000(t,S, initial_epoch, c_const), tspan, S_injection, opt_traj);

% DSM1 -> Moon (starts at initial_epoch + tof_inj2dsm1_s)
S_postDSM1 = S_preDSM1(end,:);
S_postDSM1(4:6) = S_postDSM1(4:6) + dv_dsm1;
tspan = [0 tof_dsm12flyby_s];
[~, S_preFlyby] = ode45(@(t,S) NBODY_J2000(t,S, initial_epoch + tof_inj2dsm1_s, c_const), tspan, S_postDSM1, opt_traj);


ceq = [ceq; (S_preFlyby(end,1:3) - moon_state(1:3))'./1e5];
vinf_arrival = S_preFlyby(end,4:6) - moon_state(4:6);

%% Backward propagation

% Comet to DSM2
S_comet = [target_position, S_comet_arr_vel];
tspan = [0 -tof_dsm22comet_s];
[~, S_postDSM2] = ode45(@(t,S) NBODY_J2000(t,S,epoch_comet_flyby,c_const), tspan, S_comet, opt_traj);

% DSM2
S_preDSM2 = S_postDSM2(end,:);
S_preDSM2(4:6) = S_preDSM2(4:6) - dv_dsm2;

% DSM2 to Moon
tspan = [0 -tof_flyby2dsm2_s];
[~, S_postFlyby] = ode45(@(t,S) NBODY_J2000(t,S,epoch_comet_flyby-tof_dsm22comet_s,c_const), tspan, S_preDSM2, opt_traj);

ceq = [ceq; (S_postFlyby(end,1:3) - moon_state(1:3))'./1e5];
vinf_departure = S_postFlyby(end,4:6) - moon_state(4:6);

%% Constraint on the maximum bending angle

% Maximum deflection angle (lunar gravity assist)
muMoon = 4902.800066;   % km^3/s^2
rMoon  = 1737;           % km
rM     = rMoon + hp;

vinf_avg = (norm(vinf_departure) + norm(vinf_arrival))/2; % NOTE: could be revised (e.g. use the fmincon result?)

delta_max = 2*asin( muMoon ./ (rM * vinf_avg^2 + muMoon) );
%delta_max = deg2rad(90);

cos_theta = dot(vinf_departure, vinf_arrival) / (norm(vinf_departure)*norm(vinf_arrival));
c = [c; (cos(delta_max) - cos_theta)*10];

%% Constraint on the Vinf norm
ceq = [ceq; norm(vinf_departure) - norm(vinf_arrival)];


%% Constraint on the maximum delta-v
total_dv = norm(dv_inj) + norm(dv_dsm1) + norm(dv_dsm2);
c = [c; (total_dv - max_dv*1e-3)*100];


%% A constraint on the minimum Earth flyby (e.g. 10,000 km) could be added here


end
