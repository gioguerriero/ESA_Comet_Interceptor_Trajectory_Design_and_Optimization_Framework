function [x0, f] = build_initial_guess_MS(x0, k_vec, S_traj_syn, epoch_vec, ...
                                     S_pre_dsm1_syn, S_pre_dsm2_syn, ...
                                     epoch_dep, epoch_dsm1, ...
                                     epoch_flyby, epoch_dsm2, epoch_arr)
% build_initial_guess_MS  Build the multiple-shooting initial guess by
% appending the node states of the 4 segments to x0.
%
% Inputs:
%   x0             - initial design-variable vector (column)
%   k_vec          - [k1 k2 k3 k4] number of nodes per segment
%   S_traj_syn     - (N x >=6) propagated trajectory in the synodic frame
%   epoch_vec      - (N x 1) epochs matching S_traj_syn
%   S_pre_dsm1_syn - pre-DSM1 state (6 x 1) - design variable
%   S_pre_dsm2_syn - pre-DSM2 state (6 x 1) - design variable
%   epoch_dep, epoch_dsm1, epoch_flyby, epoch_dsm2 - event epochs
%   epoch_arr      - comet arrival epoch
%
% Outputs:
%   x0 - extended vector with ONLY the node STATES of the 4 segments
%   f  - 1x4 cell, fixed node time fractions (T = Ts + f*(Te-Ts)). Node epochs
%        are NOT design variables; only the 4 event epochs (dep, dsm1, flyby,
%        dsm2) remain in the base of x0.
%
% Segment 1 (injection -> dsm1): k1 chronological nodes, injection excluded,
%                                pre_dsm1 included.
% Segment 2 (dsm1 -> flyby)    : k2 chronological nodes, interior points only.
% Segment 3 (flyby -> dsm2)    : k3 time-REVERSED nodes, post-flyby excluded,
%                                pre_dsm2 included (first).
% Segment 4 (dsm2 -> comet)    : k4 time-REVERSED nodes, interior only.

    k1 = k_vec(1);  k2 = k_vec(2);
    k3 = k_vec(3);  k4 = k_vec(4);

    %epoch_arr = epoch_vec(end);   % epoca arrivo cometa

    %% Segment 1: halo -> dsm1 (chronological) ---------------------------
    tau1 = linspace(epoch_dep, epoch_dsm1, k1 + 1);
    tau1 = tau1(2:end);
    X1 = zeros(6, k1);
    T1 = zeros(k1, 1);
    for i = 1:k1-1
        [X1(:,i), T1(i)] = pick_state(tau1(i), epoch_vec, S_traj_syn);
    end
    X1(:, k1) = S_pre_dsm1_syn(:);
    T1(k1)    = epoch_dsm1;

    %% Segment 2: dsm1 -> flyby (chronological, interior only) ------------
    tau2 = linspace(epoch_dsm1, epoch_flyby, k2 + 2);
    tau2 = tau2(2:end-1);
    X2 = zeros(6, k2);
    T2 = zeros(k2, 1);
    for i = 1:k2
        [X2(:,i), T2(i)] = pick_state(tau2(i), epoch_vec, S_traj_syn);
    end

    %% Segment 3: flyby -> dsm2 (time-reversed) --------------------
    tau3 = linspace(epoch_flyby, epoch_dsm2, k3 + 1);
    tau3 = tau3(2:end);
    tau3 = fliplr(tau3);
    X3 = zeros(6, k3);
    T3 = zeros(k3, 1);
    X3(:, 1) = S_pre_dsm2_syn(:);
    T3(1)    = epoch_dsm2;
    for i = 2:k3
        [X3(:,i), T3(i)] = pick_state(tau3(i), epoch_vec, S_traj_syn);
    end

    %% Segment 4: dsm2 -> comet (reversed, interior only) ----------------
    tau4 = linspace(epoch_dsm2, epoch_arr, k4 + 2);
    tau4 = tau4(2:end-1);
    tau4 = fliplr(tau4);
    X4 = zeros(6, k4);
    T4 = zeros(k4, 1);
    for i = 1:k4
        [X4(:,i), T4(i)] = pick_state(tau4(i), epoch_vec, S_traj_syn);
    end

    %% Fixed node time fractions + append the STATES only ------
    % Each node: T = Ts + f*(Te - Ts), with f = (T_guess - Ts)/(Te - Ts).
    %   Seg1: [dep, dsm1]   Seg2: [dsm1, flyby]
    %   Seg3: [flyby, dsm2] Seg4: [dsm2, comet]
    % Event nodes come out with f = 1 (T1(k1)=dsm1, T3(1)=dsm2), so during
    % reconstruction they land exactly on the boundary.
    f    = cell(1, 4);
    f{1} = (T1 - epoch_dep)   ./ (epoch_dsm1  - epoch_dep);
    f{2} = (T2 - epoch_dsm1)  ./ (epoch_flyby - epoch_dsm1);
    f{3} = (T3 - epoch_flyby) ./ (epoch_dsm2  - epoch_flyby);
    f{4} = (T4 - epoch_dsm2)  ./ (epoch_arr   - epoch_dsm2);

    % Node epochs are NOT design variables: append the states only.
    x0 = [x0(:);
          X1(:);
          X2(:);
          X3(:);
          X4(:)];
end

% =====================================================================
function [x_pick, t_pick] = pick_state(tau, epoch_vec, S_traj_syn)
% pick_state  Return the trajectory state (and epoch) nearest to time tau.
    [~, idx] = min(abs(epoch_vec - tau));
    x_pick = S_traj_syn(idx, 1:6).';
    t_pick = epoch_vec(idx);
end