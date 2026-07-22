function plot_guess_pre_flyby(x, c, state_tcm, epoch_tcm, epoch_flyby, state_inj, epoch_inj)
% PLOT_GUESS_PRE_FLYBY  Plot pre-flyby trajectory in Sun-Earth synodic frame.
%
% Propagates the trajectory with full ephemeris (same chain as
% NC_refinement_pre_flyby) and plots in synodic adimensional coordinates.
%
% INPUTS
%   x          [7x1]  design vector:
%                       x(1:3) – injection DV  [km/s]
%                       x(4:6) – DSM1 DV       [km/s]
%                       x(7)   – epoch_dsm1/1e8 [s/1e8]
%   c                  constants struct (needs c.mu)
%   state_tcm  [6x1]  heliocentric J2000 state at TCM1 epoch [km, km/s]
%   epoch_tcm          SPICE ET of TCM1 [s]
%   epoch_flyby        SPICE ET of flyby [s]
%   state_inj  [6x1]  heliocentric J2000 departure state [km, km/s]
%   epoch_inj          SPICE ET of departure [s]
%
% USAGE (initial guess):
%   x0 = [out.injection.J2000_kms; out.dsm1.dv.J2000_kms; out.dsm1.epoch/1e+8];
%   plot_guess_pre_flyby(x0, c, state_tcm, epoch_tcm, epoch_flyby, ...
%                        out.departure.state.J2000, out.departure.epoch)
%
% USAGE (converged solution):
%   plot_guess_pre_flyby(ultima_iter, c, state_tcm, epoch_tcm, epoch_flyby, ...
%                        out.departure.state.J2000, out.departure.epoch)

%% --- Extract design variables ----------------------------------------------
x = x(:);
inj        = x(1:3);
dsm1       = x(4:6);
epoch_dsm1 = x(7) * 1e+8;

state_inj  = state_inj(:);
state_tcm  = state_tcm(:);

%% --- Propagation -----------------------------------------------------------
opt_traj = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
N_pts    = 300;

% Arc 1: injection -> DSM1
s0_arc1 = state_inj;
s0_arc1(4:6) = s0_arc1(4:6) + inj;
tof1   = epoch_dsm1 - epoch_inj;
tspan1 = linspace(0, tof1, N_pts);
[t1, S1] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_inj, c), ...
                 tspan1, s0_arc1, opt_traj);
epochs_arc1 = epoch_inj + t1';

S_pre_dsm1  = S1(end, :);
S_post_dsm1 = S_pre_dsm1';
S_post_dsm1(4:6) = S_post_dsm1(4:6) + dsm1;

% Arc 2: DSM1 -> TCM1  (tof must be positive: TCM > DSM1)
tof2   = epoch_tcm - epoch_dsm1;
tspan2 = linspace(0, tof2, N_pts);
[t2, S2] = ode45(@(t,S) NBODY_J2000_full_ephe(t, S, epoch_dsm1, c), ...
                 tspan2, S_post_dsm1, opt_traj);
epochs_arc2 = epoch_dsm1 + t2';

S_pre_tcm_prop = S2(end, :);

%% --- Convert to synodic adimensional ---------------------------------------
SYN1 = j2000_to_synodic(S1, epochs_arc1, c.mu);
SYN2 = j2000_to_synodic(S2, epochs_arc2, c.mu);

syn_inj_post  = j2000_to_synodic(s0_arc1',       epoch_inj,    c.mu);
syn_pre_dsm1  = j2000_to_synodic(S_pre_dsm1,     epoch_dsm1,   c.mu);
syn_post_dsm1 = j2000_to_synodic(S_post_dsm1',   epoch_dsm1,   c.mu);
syn_pre_tcm   = j2000_to_synodic(S_pre_tcm_prop, epoch_tcm,    c.mu);
syn_tcm_tgt   = j2000_to_synodic(state_tcm',     epoch_tcm,    c.mu);

%% --- Moon orbit in synodic -------------------------------------------------
N_moon = 80;
t_moon = linspace(epoch_inj, epoch_flyby, N_moon);
r_moon_syn = zeros(N_moon, 6);
for i = 1:N_moon
    st = cspice_spkezr('MOON', t_moon(i), 'ECLIPJ2000', 'NONE', 'SUN');
    r_moon_syn(i, :) = j2000_to_synodic(st(1:6)', t_moon(i), c.mu);
end

st_moon_tcm   = cspice_spkezr('MOON', epoch_tcm,   'ECLIPJ2000', 'NONE', 'SUN');
st_moon_flyby = cspice_spkezr('MOON', epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN');
syn_moon_tcm   = j2000_to_synodic(st_moon_tcm',   epoch_tcm,   c.mu);
syn_moon_flyby = j2000_to_synodic(st_moon_flyby', epoch_flyby, c.mu);

%% --- Lagrange points -------------------------------------------------------
mu  = c.mu;
rL  = (mu/3)^(1/3);
L1  = [1 - mu - rL, 0, 0];
L2  = [1 - mu + rL, 0, 0];

%% --- Plot ------------------------------------------------------------------
figure('Name', 'Pre-Flyby Refinement (Synodic Adim)', 'Color', 'w');
hold on; grid on; axis equal;

plot3(SYN1(:,1), SYN1(:,2), SYN1(:,3), ...
      'b-', 'LineWidth', 1.8, 'DisplayName', 'Arc 1: inj \rightarrow DSM1');
plot3(SYN2(:,1), SYN2(:,2), SYN2(:,3), ...
      'r-', 'LineWidth', 1.8, 'DisplayName', 'Arc 2: DSM1 \rightarrow TCM1');

plot3(r_moon_syn(:,1), r_moon_syn(:,2), r_moon_syn(:,3), ...
      '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'DisplayName', 'Moon orbit');

plot3(-mu,   0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 10, 'DisplayName', 'Sun');
plot3(1-mu,  0, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize',  8, 'DisplayName', 'Earth');
plot3(L1(1), 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5,       'DisplayName', 'L1');
plot3(L2(1), 0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5,       'DisplayName', 'L2');

plot3(syn_moon_tcm(1),   syn_moon_tcm(2),   syn_moon_tcm(3), ...
      'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 7, 'DisplayName', 'Moon @ TCM1');
plot3(syn_moon_flyby(1), syn_moon_flyby(2), syn_moon_flyby(3), ...
      'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 8, 'DisplayName', 'Moon @ flyby');

plot3(syn_inj_post(1),  syn_inj_post(2),  syn_inj_post(3), ...
      'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8, 'DisplayName', 'Injection (post DV)');
plot3(syn_pre_dsm1(1),  syn_pre_dsm1(2),  syn_pre_dsm1(3), ...
      'r^', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', 'Pre-DSM1');
plot3(syn_post_dsm1(1), syn_post_dsm1(2), syn_post_dsm1(3), ...
      'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Post-DSM1');
plot3(syn_pre_tcm(1),   syn_pre_tcm(2),   syn_pre_tcm(3), ...
      'rs', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', 'Propagated pre-TCM1');
plot3(syn_tcm_tgt(1),   syn_tcm_tgt(2),   syn_tcm_tgt(3), ...
      'k*', 'MarkerSize', 10, 'DisplayName', 'Target TCM1 state');

xlabel('x [adim]'); ylabel('y [adim]'); zlabel('z [adim]');
title('Pre-Flyby Refinement — Sun-Earth Synodic Adim');
legend('show', 'Location', 'best');
view(0, 90);

%% --- Console summary -------------------------------------------------------
pos_err = S_pre_tcm_prop(1:3) - state_tcm(1:3)';
fprintf('\n--- Pre-TCM1 position mismatch ---\n');
fprintf('  |dr| = %.4e km\n', norm(pos_err));
fprintf('  TOF inj -> DSM1 : %.2f days\n', tof1/86400);
fprintf('  TOF DSM1 -> TCM1: %.2f days\n', tof2/86400);
fprintf('  |DV_inj|  = %.4f m/s\n', norm(inj)*1e3);
fprintf('  |DV_DSM1| = %.4f m/s\n\n', norm(dsm1)*1e3);

end


%% =========================================================================
function syn = j2000_to_synodic(S_J2000, epochs, mu)
% Convert N×6 heliocentric J2000 states to synodic adimensional.
% Inverse of synodic2sun_J2000, including Ldot correction.

    epochs = epochs(:)';
    N = numel(epochs);
    if size(S_J2000, 1) ~= N
        S_J2000 = S_J2000';
    end

    ES   = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');
    rE   = ES(1:3, :);
    vE   = ES(4:6, :);
    L    = sqrt(sum(rE.^2, 1));
    h    = cross(rE, vE);
    hnrm = sqrt(sum(h.^2, 1));
    om   = hnrm ./ L.^2;
    V    = L .* om;
    Ldot = dot(rE, vE, 1) ./ L;

    e1 = rE ./ L;
    e3 = h  ./ hnrm;
    e2 = cross(e3, e1);

    syn = zeros(N, 6);
    for i = 1:N
        R = [e1(:,i), e2(:,i), e3(:,i)];

        r_J = S_J2000(i, 1:3)';
        v_J = S_J2000(i, 4:6)';

        r_rot_dim = R' * r_J;
        r_syn     = r_rot_dim / L(i) - [mu; 0; 0];

        omega_vec = om(i) * e3(:, i);
        v_rot_dim = R' * (v_J - cross(omega_vec, r_J) - r_rot_dim * Ldot(i) / L(i));
        v_syn     = v_rot_dim / V(i);

        syn(i, :) = [r_syn', v_syn'];
    end
end
