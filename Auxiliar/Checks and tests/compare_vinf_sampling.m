% COMPARE_VINF_SAMPLING
% Visualizes the difference between rectangular (meshgrid) and
% polar (disk) sampling of V-infinity directions on the deflection cone.
%
% Each point is a unit direction vector, rotated by (fpa, oop) from
% a reference direction, then plotted on the unit sphere.

clear; clc; close all;

delta_max = deg2rad(70);  % example deflection angle
M = 1000;                  % target number of points (all methods aim for ~M)

%% ========================================================
%  METHOD 1: RECTANGULAR GRID (current)
% ========================================================

n_fpa = round(sqrt(M));
n_out = n_fpa;

fpa_vec = linspace(-delta_max, delta_max, n_fpa);
out_vec = linspace(-delta_max, delta_max, n_out);

[FPA_rect, OUT_rect] = meshgrid(fpa_vec, out_vec);
FPA_rect = FPA_rect(:);
OUT_rect = OUT_rect(:);

%% ========================================================
%  METHOD 2: POLAR (DISK) GRID (proposed)
% ========================================================

n_r   = round(sqrt(M));
n_ang = round(M / n_r);

r_vec   = linspace(0, delta_max, n_r);
ang_vec = linspace(0, 2*pi, n_ang+1); ang_vec(end) = [];

[R, ANG] = meshgrid(r_vec, ang_vec);
FPA_disk = R(:) .* cos(ANG(:));
OUT_disk = R(:) .* sin(ANG(:));

%% ========================================================
%  METHOD 3: POLAR (DISK) GRID — UNIFORM AREA (variable angular points)
% ========================================================

n_r_u = round(sqrt(M));
r_vec_u = linspace(0, delta_max, n_r_u);   % equispaced radii



% Distribute ~M points: 1 at center, rest proportional to ring index k
ring_weights = (1:n_r_u-1);                        % weight ~ k (proportional to radius)
pts_per_ring = round((M - 1) * ring_weights / sum(ring_weights));
pts_per_ring = max(pts_per_ring, 3);                % at least 3 per ring

FPA_unif = 0;   % center point
OUT_unif = 0;

for k = 2:n_r_u
    n_ang_k = pts_per_ring(k-1);
    ang_k   = linspace(0, 2*pi, n_ang_k+1); ang_k(end) = [];
    FPA_unif = [FPA_unif; r_vec_u(k) * cos(ang_k')];
    OUT_unif = [OUT_unif; r_vec_u(k) * sin(ang_k')];
end

%% ========================================================
%  CONVERT ANGLES TO 3D UNIT VECTORS
% ========================================================
% Reference direction: [1, 0, 0] (arbitrary, just for visualization)
% fpa  -> rotation in XY plane (around Z)
% oop  -> rotation toward Z (around axis perpendicular to rotated vector)

dir_rect = angles_to_vectors(FPA_rect, OUT_rect);
dir_disk = angles_to_vectors(FPA_disk, OUT_disk);
dir_unif = angles_to_vectors(FPA_unif, OUT_unif);

%% ========================================================
%  PLOT
% ========================================================

figure('Color','w','Position',[100 100 1800 500]);

[xs, ys, zs] = sphere(40);

% --- Method 1: Rectangular ---
subplot(1,3,1);
hold on; grid on; axis equal;
surf(xs, ys, zs, 'FaceAlpha',0.05, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);
quiver3(0,0,0, 1.3,0,0, 0, 'k','LineWidth',2);
scatter3(dir_rect(:,1), dir_rect(:,2), dir_rect(:,3), ...
         30, 'b', 'filled', 'MarkerFaceAlpha',0.7);
title(sprintf('Rectangular (%d pts)', length(FPA_rect)));
xlabel('X'); ylabel('Y'); zlabel('Z');
view(30, 25);

% --- Method 2: Polar (linear r) ---
subplot(1,3,2);
hold on; grid on; axis equal;
surf(xs, ys, zs, 'FaceAlpha',0.05, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);
quiver3(0,0,0, 1.3,0,0, 0, 'k','LineWidth',2);
scatter3(dir_disk(:,1), dir_disk(:,2), dir_disk(:,3), ...
         30, 'r', 'filled', 'MarkerFaceAlpha',0.7);
title(sprintf('Polar linear r (%d pts)', length(FPA_disk)));
xlabel('X'); ylabel('Y'); zlabel('Z');
view(30, 25);

% --- Method 3: Polar (sqrt r — uniform area) ---
subplot(1,3,3);
hold on; grid on; axis equal;
surf(xs, ys, zs, 'FaceAlpha',0.05, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);
quiver3(0,0,0, 1.3,0,0, 0, 'k','LineWidth',2);
scatter3(dir_unif(:,1), dir_unif(:,2), dir_unif(:,3), ...
         30, [0 0.7 0], 'filled', 'MarkerFaceAlpha',0.7);
title(sprintf('Polar uniform (%d pts)', length(FPA_unif)));
xlabel('X'); ylabel('Y'); zlabel('Z');
view(30, 25);

sgtitle(sprintf('V_{\\infty} direction sampling | \\delta_{max} = %.1f°', rad2deg(delta_max)), ...
        'FontWeight','bold','FontSize',14);

%% ========================================================
%  THIRD FIGURE: OVERLAY + ANGLE SPACE
% ========================================================

figure('Color','w','Position',[100 100 1400 600]);

% --- 3D overlay ---
subplot(1,2,1);
hold on; grid on; axis equal;

surf(xs, ys, zs, 'FaceAlpha',0.05, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);
quiver3(0,0,0, 1.3,0,0, 0, 'k','LineWidth',2);

scatter3(dir_rect(:,1), dir_rect(:,2), dir_rect(:,3), ...
         30, 'b', 'filled', 'MarkerFaceAlpha',0.4);
scatter3(dir_disk(:,1), dir_disk(:,2), dir_disk(:,3), ...
         30, 'r', 'filled', 'MarkerFaceAlpha',0.4);
scatter3(dir_unif(:,1), dir_unif(:,2), dir_unif(:,3), ...
         30, [0 0.7 0], 'filled', 'MarkerFaceAlpha',0.4);

legend('','Reference dir','Rectangular','Polar (linear)','Polar (sqrt)','Location','best');
title('3D overlay');
xlabel('X'); ylabel('Y'); zlabel('Z');
view(30, 25);

% --- 2D angle space ---
subplot(1,2,2);
hold on; grid on; axis equal;

th = linspace(0, 2*pi, 200);
plot(rad2deg(delta_max)*cos(th), rad2deg(delta_max)*sin(th), 'k--','LineWidth',1.5);

scatter(rad2deg(FPA_rect), rad2deg(OUT_rect), 30, 'b', 'filled', 'MarkerFaceAlpha',0.4);
scatter(rad2deg(FPA_disk), rad2deg(OUT_disk), 30, 'r', 'filled', 'MarkerFaceAlpha',0.4);
scatter(rad2deg(FPA_unif), rad2deg(OUT_unif), 30, [0 0.7 0], 'filled', 'MarkerFaceAlpha',0.4);

legend('\delta_{max} circle','Rectangular','Polar (linear)','Polar (sqrt)','Location','best');
title('Angle space (fpa vs oop)');
xlabel('fpa [deg]'); ylabel('oop [deg]');

sgtitle('Comparison of sampling methods','FontWeight','bold','FontSize',14);


%% ========================================================
%  FINAL PLOT: 2000 points, delta_max = 70 deg (polar uniform)
% ========================================================

delta_max_final = deg2rad(70);
M_final = 2000;

n_r_f = round(sqrt(M_final));
r_vec_f = linspace(0, delta_max_final, n_r_f);

ring_weights_f = (1:n_r_f-1);
pts_per_ring_f = round((M_final - 1) * ring_weights_f / sum(ring_weights_f));
pts_per_ring_f = max(pts_per_ring_f, 3);

FPA_final = 0;
OUT_final = 0;

for k = 2:n_r_f
    n_ang_k = pts_per_ring_f(k-1);
    ang_k   = linspace(0, 2*pi, n_ang_k+1); ang_k(end) = [];
    FPA_final = [FPA_final; r_vec_f(k) * cos(ang_k')];
    OUT_final = [OUT_final; r_vec_f(k) * sin(ang_k')];
end

dir_final = angles_to_vectors(FPA_final, OUT_final);

figure('Color','w','Position',[100 100 800 700]);
hold on; grid on; axis equal;

h_sphere = surf(xs, ys, zs, 'FaceAlpha',0.15, 'EdgeColor','none', 'FaceColor',[0.8 0.8 1]);
h_ref    = quiver3(0,0,0, 1.3,0,0, 0, 'k','LineWidth',2.5);
h_pts    = scatter3(dir_final(:,1), dir_final(:,2), dir_final(:,3), ...
                    25, 'r', 'filled', 'MarkerFaceAlpha',0.7);

legend([h_ref, h_pts], ...
       {'Reference $v_{\infty}$', 'Sampled $v_{\infty}$ directions'}, ...
       'Interpreter','latex','Location','best','FontSize',12);

xlabel('X'); ylabel('Y'); zlabel('Z');
view(30, 25);

%% ========================================================
%  LOCAL FUNCTION
% ========================================================

function dirs = angles_to_vectors(fpa, oop)
% Converts (fpa, oop) angle pairs to 3D unit direction vectors.
% Reference direction: [1, 0, 0]
% fpa: clockwise rotation in XY plane
% oop: rotation toward +Z

    ref = [1; 0; 0];
    N = length(fpa);
    dirs = zeros(N, 3);

    for k = 1:N
        % In-plane rotation (around Z, clockwise -> negative angle)
        Rz = [ cos(-fpa(k))  -sin(-fpa(k))  0;
               sin(-fpa(k))   cos(-fpa(k))  0;
               0              0              1];

        v_xy = Rz * ref;

        % Out-of-plane rotation axis: cross(v_xy, [0;0;1])
        rot_ax = cross(v_xy, [0;0;1]);
        rot_ax = rot_ax / norm(rot_ax);

        % Rodrigues rotation
        K = [     0        -rot_ax(3)   rot_ax(2);
              rot_ax(3)      0         -rot_ax(1);
             -rot_ax(2)   rot_ax(1)      0       ];

        R_oop = eye(3) + sin(oop(k))*K + (1 - cos(oop(k)))*(K*K);

        v_3d = R_oop * v_xy;
        dirs(k,:) = v_3d';
    end
end
