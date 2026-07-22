function plot_traj_post_flyby(x_guess, x_sol, c, state_tcm, dv_tcm, ...
                               epoch_tcm, t_tcm2, epoch_flyby, comet_pos, epoch_comet)
% plot_traj_post_flyby  Plot the guess and the solution of the post-flyby
% refinement in the non-dimensional Sun-Earth synodic frame.
%
% The propagation chain is identical to NC_refinement_post_flyby.
% The TCM1->TCM2 arc is fixed (same inputs) and is plotted only once.
% From TCM2 onward the guess (dashed) and the solution (solid) diverge.
%
% Inputs:
%   x_guess [7x1] - initial guess passed to refinement_post_flyby
%   x_sol   [7x1] - converged solution vector
%   c             - constants struct (c.mu, ...)
%   state_tcm     - J2000 state at TCM1 epoch [km, km/s] (pre-maneuver)
%   dv_tcm  [3x1] - TCM1 delta-v [km/s]
%   epoch_tcm     - SPICE ET of TCM1 [s]
%   t_tcm2        - days after flyby when TCM2 is applied [days]
%   epoch_flyby   - SPICE ET of flyby [s]
%   comet_pos     - target comet position [km] (3x1 or 1x3)
%   epoch_comet   - SPICE ET of comet rendezvous [s]
%
% Outputs:
%   (none) - produces a figure
%
% USAGE
%   x0   = [tcm2_guess; out.dsm2.dv.J2000_kms; out.dsm2.epoch/1e+8];
%   plot_traj_post_flyby(x0, ultima_iter_post, c, state_tcm, dv_tcm, ...
%                        epoch_tcm, t_tcm2, epoch_flyby, ...
%                        selected_comet.comet_pos, selected_comet.epoch)

%% --- Setup -----------------------------------------------------------------
opt  = odeset('AbsTol', 1e-13, 'RelTol', 1e-13);
N    = 300;
mu   = c.mu;

state_tcm = state_tcm(:);
dv_tcm    = dv_tcm(:);

epoch_tcm2 = epoch_flyby + t_tcm2 * 86400;

%% --- Arc 0: TCM1 → TCM2 (identical for guess and solution) -----------------
s0_arc0 = state_tcm;
s0_arc0(4:6) = s0_arc0(4:6) + dv_tcm;

tof0   = epoch_tcm2 - epoch_tcm;
tspan0 = linspace(0, tof0, N);
[t0, S0] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm, c), ...
                 tspan0, s0_arc0, opt);
epochs0 = epoch_tcm + t0';

S_pre_tcm2 = S0(end, :);   % same for both

%% --- Propagate guess and solution from TCM2 onward ------------------------
[arcs_g, nodes_g] = propagate_post(x_guess, S_pre_tcm2, epoch_tcm2, epoch_comet, c, opt, N);
[arcs_s, nodes_s] = propagate_post(x_sol,   S_pre_tcm2, epoch_tcm2, epoch_comet, c, opt, N);

%% --- Convert all arcs to synodic adimensional ------------------------------
SYN0     = j2000_to_synodic(S0,          epochs0,                mu);
SYN1_g   = j2000_to_synodic(arcs_g.S1,  arcs_g.epochs1,         mu);
SYN2_g   = j2000_to_synodic(arcs_g.S2,  arcs_g.epochs2,         mu);
SYN1_s   = j2000_to_synodic(arcs_s.S1,  arcs_s.epochs1,         mu);
SYN2_s   = j2000_to_synodic(arcs_s.S2,  arcs_s.epochs2,         mu);

% Key nodes
syn_tcm1_post = j2000_to_synodic(s0_arc0.',         epoch_tcm,          mu);
syn_tcm2      = j2000_to_synodic(S_pre_tcm2,        epoch_tcm2,         mu);

syn_dsm2_pre_g  = j2000_to_synodic(nodes_g.S_pre_dsm2,  nodes_g.epoch_dsm2, mu);
syn_dsm2_post_g = j2000_to_synodic(nodes_g.S_post_dsm2.', nodes_g.epoch_dsm2, mu);
syn_comet_g     = j2000_to_synodic(nodes_g.S_comet.',    epoch_comet,        mu);

syn_dsm2_pre_s  = j2000_to_synodic(nodes_s.S_pre_dsm2,  nodes_s.epoch_dsm2, mu);
syn_dsm2_post_s = j2000_to_synodic(nodes_s.S_post_dsm2.', nodes_s.epoch_dsm2, mu);
syn_comet_s     = j2000_to_synodic(nodes_s.S_comet.',    epoch_comet,        mu);

syn_comet_tgt   = j2000_to_synodic([comet_pos(:); zeros(3,1)].', epoch_comet, mu);

%% --- Moon and comet reference orbit ----------------------------------------
N_bod = 100;
t_bod = linspace(epoch_tcm, epoch_comet, N_bod);
r_moon_syn = zeros(N_bod, 6);
for i = 1:N_bod
    st = cspice_spkezr('MOON', t_bod(i), 'ECLIPJ2000', 'NONE', 'SUN');
    r_moon_syn(i, :) = j2000_to_synodic(st', t_bod(i), mu);
end

st_moon_flyby = cspice_spkezr('MOON', epoch_flyby,  'ECLIPJ2000', 'NONE', 'SUN');
st_moon_tcm2  = cspice_spkezr('MOON', epoch_tcm2,   'ECLIPJ2000', 'NONE', 'SUN');
syn_moon_flyby = j2000_to_synodic(st_moon_flyby', epoch_flyby, mu);
syn_moon_tcm2  = j2000_to_synodic(st_moon_tcm2',  epoch_tcm2,  mu);

%% --- Lagrange points -------------------------------------------------------
rL = (mu/3)^(1/3);
L1 = [1-mu-rL, 0, 0];
L2 = [1-mu+rL, 0, 0];

%% --- Plot ------------------------------------------------------------------
figure('Name', 'Post-Flyby Refinement — Guess vs Solution (Synodic Adim)', ...
       'Color', 'w');
hold on; grid on; axis equal;

% Arc 0: TCM1 -> TCM2 (common)
plot3(SYN0(:,1), SYN0(:,2), SYN0(:,3), ...
      '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Arc TCM1 \rightarrow TCM2 (common)');

% Guess arcs
plot3(SYN1_g(:,1), SYN1_g(:,2), SYN1_g(:,3), ...
      'b--', 'LineWidth', 1.5, 'DisplayName', 'Guess: TCM2 \rightarrow DSM2');
plot3(SYN2_g(:,1), SYN2_g(:,2), SYN2_g(:,3), ...
      'r--', 'LineWidth', 1.5, 'DisplayName', 'Guess: DSM2 \rightarrow comet');

% Solution arcs
plot3(SYN1_s(:,1), SYN1_s(:,2), SYN1_s(:,3), ...
      'b-', 'LineWidth', 2, 'DisplayName', 'Solution: TCM2 \rightarrow DSM2');
plot3(SYN2_s(:,1), SYN2_s(:,2), SYN2_s(:,3), ...
      'r-', 'LineWidth', 2, 'DisplayName', 'Solution: DSM2 \rightarrow comet');

% Moon orbit (reference)
plot3(r_moon_syn(:,1), r_moon_syn(:,2), r_moon_syn(:,3), ...
      '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'Moon orbit');

% Fixed bodies
plot3(-mu,   0, 0, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 10, 'DisplayName', 'Sun');
plot3(1-mu,  0, 0, 'go', 'MarkerFaceColor', 'g', 'MarkerSize',  8, 'DisplayName', 'Earth');
plot3(L1(1), 0, 0, 'k+', 'MarkerSize', 10, 'LineWidth', 1.5,       'DisplayName', 'L1');
plot3(L2(1), 0, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.5,       'DisplayName', 'L2');

% Moon at flyby and TCM2
plot3(syn_moon_flyby(1), syn_moon_flyby(2), syn_moon_flyby(3), ...
      'm^', 'MarkerFaceColor', 'm', 'MarkerSize', 9, 'DisplayName', 'Moon @ flyby');
plot3(syn_moon_tcm2(1),  syn_moon_tcm2(2),  syn_moon_tcm2(3), ...
      'ms', 'MarkerFaceColor', 'm', 'MarkerSize', 7, 'DisplayName', 'Moon @ TCM2');

% TCM1 post-maneuver (start of arc 0)
plot3(syn_tcm1_post(1), syn_tcm1_post(2), syn_tcm1_post(3), ...
      'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'DisplayName', 'TCM1 (post-DV)');

% TCM2 (same for both)
plot3(syn_tcm2(1), syn_tcm2(2), syn_tcm2(3), ...
      'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'DisplayName', 'TCM2');

% DSM2 guess
plot3(syn_dsm2_pre_g(1),  syn_dsm2_pre_g(2),  syn_dsm2_pre_g(3), ...
      'b^', 'MarkerFaceColor', 'w', 'MarkerSize', 7, 'DisplayName', 'Pre-DSM2 (guess)');
plot3(syn_dsm2_post_g(1), syn_dsm2_post_g(2), syn_dsm2_post_g(3), ...
      'bv', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Post-DSM2 (guess)');

% DSM2 solution
plot3(syn_dsm2_pre_s(1),  syn_dsm2_pre_s(2),  syn_dsm2_pre_s(3), ...
      'r^', 'MarkerFaceColor', 'w', 'MarkerSize', 7, 'DisplayName', 'Pre-DSM2 (sol)');
plot3(syn_dsm2_post_s(1), syn_dsm2_post_s(2), syn_dsm2_post_s(3), ...
      'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'DisplayName', 'Post-DSM2 (sol)');

% Comet guess and solution endpoints
plot3(syn_comet_g(1), syn_comet_g(2), syn_comet_g(3), ...
      'b*', 'MarkerSize', 10, 'DisplayName', 'Comet arrival (guess)');
plot3(syn_comet_s(1), syn_comet_s(2), syn_comet_s(3), ...
      'r*', 'MarkerSize', 10, 'DisplayName', 'Comet arrival (sol)');

% Comet target position
plot3(syn_comet_tgt(1), syn_comet_tgt(2), syn_comet_tgt(3), ...
      'k*', 'MarkerSize', 12, 'LineWidth', 1.5, 'DisplayName', 'Comet target');

xlabel('x [adim]'); ylabel('y [adim]'); zlabel('z [adim]');
title('Post-Flyby Refinement — Guess (--) vs Solution (—)  |  Synodic Adim');
legend('show', 'Location', 'best');
view(0, 90);

%% --- Console summary -------------------------------------------------------
fprintf('\n--- Post-flyby refinement summary ---\n');
fprintf('  TCM2  epoch       : %s\n', cspice_et2utc(epoch_tcm2, 'C', 0));
fprintf('  DSM2  epoch guess : %s\n', cspice_et2utc(nodes_g.epoch_dsm2, 'C', 0));
fprintf('  DSM2  epoch sol   : %s\n', cspice_et2utc(nodes_s.epoch_dsm2, 'C', 0));
fprintf('  |DV_TCM2| guess   : %.4f m/s\n', norm(x_guess(1:3))*1e3);
fprintf('  |DV_TCM2| sol     : %.4f m/s\n', norm(x_sol(1:3))*1e3);
fprintf('  |DV_DSM2| guess   : %.4f m/s\n', norm(x_guess(4:6))*1e3);
fprintf('  |DV_DSM2| sol     : %.4f m/s\n', norm(x_sol(4:6))*1e3);
fprintf('  Comet pos err guess: %.4e km\n', norm(nodes_g.S_comet(1:3) - comet_pos(:)'));
fprintf('  Comet pos err sol  : %.4e km\n', norm(nodes_s.S_comet(1:3) - comet_pos(:)'));
fprintf('\n');

end


%% =========================================================================
function [arcs, nodes] = propagate_post(x, S_pre_tcm2, epoch_tcm2, epoch_comet, c, opt, N)
% Propagate arcs 1 and 2 from TCM2 given design vector x.

    tcm2       = x(1:3);
    dsm2       = x(4:6);
    epoch_dsm2 = x(7) * 1e8;

    S_post_tcm2 = S_pre_tcm2(:);
    S_post_tcm2(4:6) = S_post_tcm2(4:6) + tcm2;

    % Arc 1: TCM2 -> DSM2
    tof1   = epoch_dsm2 - epoch_tcm2;
    tspan1 = linspace(0, tof1, N);
    [t1, S1] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm2, c), ...
                     tspan1, S_post_tcm2, opt);
    epochs1 = epoch_tcm2 + t1';

    S_pre_dsm2  = S1(end, :);
    S_post_dsm2 = S_pre_dsm2(:);
    S_post_dsm2(4:6) = S_post_dsm2(4:6) + dsm2;

    % Arc 2: DSM2 -> comet
    tof2   = epoch_comet - epoch_dsm2;
    tspan2 = linspace(0, tof2, N);
    [t2, S2] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_dsm2, c), ...
                     tspan2, S_post_dsm2, opt);
    epochs2 = epoch_dsm2 + t2';

    arcs.S1      = S1;
    arcs.epochs1 = epochs1;
    arcs.S2      = S2;
    arcs.epochs2 = epochs2;

    nodes.epoch_dsm2  = epoch_dsm2;
    nodes.S_pre_dsm2  = S_pre_dsm2;
    nodes.S_post_dsm2 = S_post_dsm2;
    nodes.S_comet     = S2(end, :);
end


%% =========================================================================
function syn = j2000_to_synodic(S_J2000, epochs, mu)
% Convert N×6 heliocentric J2000 states to synodic adimensional.

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
