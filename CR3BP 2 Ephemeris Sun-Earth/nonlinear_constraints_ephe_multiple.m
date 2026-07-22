function [c, ceq] = nonlinear_constraints_ephe_multiple(x, S_halo, max_dv, target_position, epoch_comet_flyby, hp, c_const, max_dv_inj, k_vec, f_nodes)
% nonlinear_constraints_ephe_multiple  Constraints for the direct ephemeris
% optimization with MULTIPLE SHOOTING.
%
% x layout (17 base + 6*sum(k_vec) node STATES):
%   x(1:3)   - dv_inj * 1e+3            [m/s]
%   x(4:6)   - dv_dsm1 * 10             [km/s*10]
%   x(7)     - epoch_flyby  [days before comet]  (epoch = epoch_comet - x*86400)
%   x(8:10)  - S_comet_arr_vel          [km/s]
%   x(11:13) - dv_dsm2 * 10             [km/s*10]
%   x(14)    - t_halo                   [days]
%   x(15)    - epoch_dep (initial_epoch) [days before comet]
%   x(16)    - epoch_dsm1               [days before comet]
%   x(17)    - epoch_dsm2               [days before comet]
%   -- then, in order, 6*ki STATES per segment (synodic non-dim). Node epochs
%      are NOT in x: they are reconstructed from f_nodes (fixed fractions) as
%      T = Ts + f*(Te - Ts).
%     Segment 1 (k1): injection -> pre_DSM1, chronological, last node = pre_DSM1
%     Segment 2 (k2): post_DSM1 -> flyby,    chronological, interior points only
%     Segment 3 (k3): post_flyby -> pre_DSM2, reversed,     first node = pre_DSM2
%     Segment 4 (k4): post_DSM2 -> comet,    reversed,      interior points only
%
% Inputs:
%   x                 - design vector (see layout above)
%   S_halo            - halo orbit states
%   max_dv            - total delta-v budget [m/s]
%   target_position   - target comet position [km]
%   epoch_comet_flyby - comet arrival epoch [ET s]
%   hp                - minimum lunar flyby altitude [km]
%   c_const           - constants struct
%   max_dv_inj        - maximum injection delta-v [m/s]
%   k_vec             - [k1 k2 k3 k4] nodes per segment
%   f_nodes           - fixed node time fractions
%
% Outputs:
%   c   - inequality constraints (<= 0)
%   ceq - equality constraints (= 0): node continuity and closure

c   = [];
ceq = [];

% --- Constraint normalization scales (must match ConstraintTolerance=1e-3
%     in optimization_ephe) ---
pos_tol_km = 10;                      % desired position continuity [km]
pos_tol_km_moon = 1;
vel_tol_ms = 0.1;                     % desired velocity continuity [m/s]
tol_feas   = 1e-3;                    % = ConstraintTolerance
Lscale = pos_tol_km / tol_feas;       % = 1e5  [km]
Lscale_moon = pos_tol_km_moon / tol_feas;
Vscale = (vel_tol_ms*1e-3) / tol_feas;% = 0.1  [km/s]

%% ===== UNPACK BASE DESIGN VARIABLES =================================
dv_inj          = x(1:3) ./ 1e+3;          % [km/s] column
dv_dsm1         = x(4:6) ./ 10;            % [km/s] column
epoch_flyby     = epoch_comet_flyby - x(7) * 86400;     % [ET s]  (x in days before comet)
S_comet_arr_vel = x(8:10).';               % row [km/s]
dv_dsm2         = x(11:13) ./ 10;          % [km/s] column
tof_onHalo_ad   = x(14) * (24*3600) / c_const.Tstar;   % non-dim CR3BP
initial_epoch   = epoch_comet_flyby - x(15) * 86400;    % [ET s] injection epoch (post halo waiting)
epoch_dsm1      = epoch_comet_flyby - x(16) * 86400;    % [ET s] DSM1 epoch
epoch_dsm2      = epoch_comet_flyby - x(17) * 86400;    % [ET s] DSM2 epoch

%% ===== UNPACK MULTIPLE SHOOTING NODES ===============================
k1 = k_vec(1);  k2 = k_vec(2);
k3 = k_vec(3);  k4 = k_vec(4);

% Node states (node epochs are NOT in x)
n_base = 17;
off1 = n_base;          off2 = off1 + 6*k1;
off3 = off2 + 6*k2;     off4 = off3 + 6*k3;
S1_syn = reshape(x(off1+1 : off1+6*k1), 6, k1);
S2_syn = reshape(x(off2+1 : off2+6*k2), 6, k2);
S3_syn = reshape(x(off3+1 : off3+6*k3), 6, k3);
S4_syn = reshape(x(off4+1 : off4+6*k4), 6, k4);

% Node epochs reconstructed from fixed fractions: T = Ts + f*(Te - Ts)
%   Seg1: [dep, dsm1]   Seg2: [dsm1, flyby]
%   Seg3: [flyby, dsm2] Seg4: [dsm2, comet]
% Event nodes come out exact: T1(k1)=epoch_dsm1, T3(1)=epoch_dsm2 (f=1).
T1 = initial_epoch + f_nodes{1}(:) .* (epoch_dsm1        - initial_epoch);
T2 = epoch_dsm1    + f_nodes{2}(:) .* (epoch_flyby       - epoch_dsm1);
T3 = epoch_flyby   + f_nodes{3}(:) .* (epoch_dsm2        - epoch_flyby);
T4 = epoch_dsm2    + f_nodes{4}(:) .* (epoch_comet_flyby - epoch_dsm2);

% Convert synodic -> J2000 at each corresponding epoch
S1 = nodes_syn2J2000(S1_syn, T1, c_const.mu);
S2 = nodes_syn2J2000(S2_syn, T2, c_const.mu);
S3 = nodes_syn2J2000(S3_syn, T3, c_const.mu);
S4 = nodes_syn2J2000(S4_syn, T4, c_const.mu);

%% ===== Moon state at flyby ==========================================
moon_state = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN').';   % 1×6

%% ===== HALO WAITING + INJECTION =====================================
opt_traj = odeset('AbsTol',1e-13,'RelTol',1e-12);

if abs(tof_onHalo_ad) > 3600/c_const.Tstar
    [~, S_partenza] = ode45(@(t,S) CR3BP(t,S,c_const.mu), [0 tof_onHalo_ad], S_halo(:), opt_traj);
    S_halo_end_syn = S_partenza(end,:);
else
    S_halo_end_syn = S_halo(:).';
end

S_halo_J2000 = synodic2sun_J2000(S_halo_end_syn, initial_epoch, c_const.mu);
S_halo_J2000 = S_halo_J2000(:);            % 6×1

% Injection
S_injection = S_halo_J2000;
S_injection(4:6) = S_injection(4:6) + dv_inj;
c = [c; norm(dv_inj)*1e+3/max_dv_inj - 1];

%% ===== SEGMENT 1: forward (injection -> pre_DSM1) =====================
% Sub-arcs: k1 total
%   arc 1: S_injection (at initial_epoch) -> S1(:,1) (at T1(1))
%   arc i (i=2..k1): S1(:,i-1) -> S1(:,i)
% 6D continuity at the end of each sub-arc.

for i = 1:k1
    if i == 1
        S_start = S_injection;   T_start = initial_epoch;
    else
        S_start = S1(:,i-1);     T_start = T1(i-1);
    end
    dt = T1(i) - T_start;
    [~, S_prop] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);
    gap = S_prop(end,:).' - S1(:,i);
    ceq = [ceq; gap(1:3)./Lscale; gap(4:6)./Vscale];
end

% Pre-DSM1 = last node of segment 1; apply the DSM1 maneuver
S_postDSM1 = S1(:,k1);
S_postDSM1(4:6) = S_postDSM1(4:6) + dv_dsm1;
T_postDSM1 = T1(k1);

%% ===== SEGMENT 2: forward (post_DSM1 -> flyby) ========================
% Sub-arcs: k2+1 total (k2 interior nodes)
%   arc 1: S_postDSM1 -> S2(:,1)
%   arc i (i=2..k2): S2(:,i-1) -> S2(:,i)
%   arc k2+1: S2(:,k2) -> flyby (position constraint = Moon)

for i = 1:k2
    if i == 1
        S_start = S_postDSM1;    T_start = T_postDSM1;
    else
        S_start = S2(:,i-1);     T_start = T2(i-1);
    end
    dt = T2(i) - T_start;
    [~, S_prop] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);
    gap = S_prop(end,:).' - S2(:,i);
    ceq = [ceq; gap(1:3)./Lscale; gap(4:6)./Vscale];
end

% Final sub-arc -> Moon
if k2 > 0
    S_start = S2(:,k2);   T_start = T2(k2);
else
    S_start = S_postDSM1; T_start = T_postDSM1;
end
dt = epoch_flyby - T_start;
[~, S_preFlyby] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);

ceq = [ceq; (S_preFlyby(end,1:3) - moon_state(1:3)).' ./ Lscale_moon];
vinf_arrival = S_preFlyby(end,4:6) - moon_state(4:6);

%% ===== SEGMENT 4: backward (comet -> S_postDSM2) ======================
% Nodes in REVERSE order: S4(:,1) is the most recent (near the comet).
% Sub-arcs: k4+1 total
%   arc 1: S_comet (at epoch_comet_flyby) -> S4(:,1) (at T4(1))
%   arc i (i=2..k4): S4(:,i-1) -> S4(:,i)
%   arc k4+1: S4(:,k4) -> S_postDSM2  (= pre_DSM2 + dv_dsm2)

S_comet = [target_position(:); S_comet_arr_vel(:)];

for i = 1:k4
    if i == 1
        S_start = S_comet;     T_start = epoch_comet_flyby;
    else
        S_start = S4(:,i-1);   T_start = T4(i-1);
    end
    dt = T4(i) - T_start;       % negative (backward)
    [~, S_prop] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);
    gap = S_prop(end,:).' - S4(:,i);
    ceq = [ceq; gap(1:3)./Lscale; gap(4:6)./Vscale];
end

% Final sub-arc: last node of segment 4 -> S_postDSM2
if k4 > 0
    S_start = S4(:,k4);    T_start = T4(k4);
else
    S_start = S_comet;     T_start = epoch_comet_flyby;
end
T_end_dsm2 = T3(1);                         % pre_DSM2 epoch (= DSM2 epoch)
dt = T_end_dsm2 - T_start;                  % negative
[~, S_postDSM2_prop] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);

% Continuity with S_postDSM2 = pre_DSM2 + dv_dsm2
S_postDSM2_target = S3(:,1);
S_postDSM2_target(4:6) = S_postDSM2_target(4:6) + dv_dsm2;
gap = S_postDSM2_prop(end,:).' - S_postDSM2_target;
ceq = [ceq; gap(1:3)./Lscale; gap(4:6)./Vscale];

%% ===== SEGMENT 3: backward (pre_DSM2 -> post_flyby) ===================
% Nodes in REVERSE order: S3(:,1) = pre_DSM2 (most recent),
%                         S3(:,k3) = node closest to the flyby.
% Sub-arcs: k3 total
%   arc i (i=1..k3-1): S3(:,i) -> S3(:,i+1)
%   arc k3: S3(:,k3) -> flyby (position constraint = Moon)

for i = 1:k3-1
    S_start = S3(:,i);     T_start = T3(i);
    dt = T3(i+1) - T_start;     % negative
    [~, S_prop] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);
    gap = S_prop(end,:).' - S3(:,i+1);
    ceq = [ceq; gap(1:3)./Lscale; gap(4:6)./Vscale];
end

% Final sub-arc -> Moon (backward)
S_start = S3(:,k3);   T_start = T3(k3);
dt = epoch_flyby - T_start;     % negative
[~, S_postFlyby] = ode45(@(t,S) NBODY_J2000(t,S, T_start, c_const), [0 dt], S_start, opt_traj);

ceq = [ceq; (S_postFlyby(end,1:3) - moon_state(1:3)).' ./ Lscale_moon];
vinf_departure = S_postFlyby(end,4:6) - moon_state(4:6);

%% ===== CONSTRAINT BENDING ANGLE (lunar gravity assist) ==============
muMoon = 4902.800066;
rMoon  = 1737;
rM     = rMoon + hp;

vinf_avg  = (norm(vinf_departure) + norm(vinf_arrival))/2;
delta_max = 2*asin( muMoon ./ (rM * vinf_avg^2 + muMoon) );
cos_theta = dot(vinf_departure, vinf_arrival) / (norm(vinf_departure)*norm(vinf_arrival));
c = [c; (cos(delta_max) - cos_theta)];

%% ===== CONSTRAINT: match |Vinf| in and out =========================
ceq = [ceq; (norm(vinf_departure) - norm(vinf_arrival))];

%% ===== CONSTRAINT: maximum total delta-v ===========================
total_dv = norm(dv_inj) + norm(dv_dsm1) + norm(dv_dsm2);
c = [c; (total_dv - max_dv*1e-3)];

end


% =====================================================================
% Helper: converte k stati sinodico-adim → J2000 (km, km/s) usando
%         l'epoca specifica di ogni nodo
% =====================================================================
function S_J = nodes_syn2J2000(S_syn, T, mu)
    k = size(S_syn, 2);
    S_J = zeros(6, k);
    for i = 1:k
        tmp = synodic2sun_J2000(S_syn(:,i).', T(i), mu);
        S_J(:,i) = tmp(:);
    end
end
