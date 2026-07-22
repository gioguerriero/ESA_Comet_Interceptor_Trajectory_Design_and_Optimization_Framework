%% TEST OF VINF_ROTATION WITH UPDATED SIGN CONVENTIONS
% Standalone sanity-check script: exercises vinf_rotation with a dummy synodic
% Moon state and visualizes the resulting v-infinity directions.
clear; clc; close all;

%% -----------------------------
% Dummy synodic Moon state
% (only velocity is used)
%% -----------------------------
synodic_moon = [0 0 0  -1 0 0]';   % velocity along +X in synodic frame

%% Parameters
Vinf = 1;                       % scale factor for visualization
fpa  = deg2rad(30);             % in-plane rotation (clockwise)
out_of_plane = deg2rad(20);     % out-of-plane rotation (toward +Z)

%% -----------------------------
% Call the rotation function
%% -----------------------------
v_inf_vec = vinf_rotation(synodic_moon, Vinf, fpa, out_of_plane);

%% -----------------------------
% For visualization, also compute intermediate in-plane rotation
%% -----------------------------
v0 = synodic_moon(4:6);
v0 = v0 / norm(v0);

% XY rotation only (clockwise)
Rz = [ cos(-fpa) -sin(-fpa) 0;
       sin(-fpa)  cos(-fpa) 0;
       0          0         1 ];
v_xy = Rz * v0;

%% -----------------------------
% 3D Plot
%% -----------------------------
figure;
hold on; grid on; axis equal;

% Original velocity direction
quiver3(0,0,0, v0(1), v0(2), v0(3), ...
        'LineWidth', 2, 'MaxHeadSize', 0.5);

% After XY rotation (clockwise)
quiver3(0,0,0, v_xy(1), v_xy(2), v_xy(3), ...
        'LineWidth', 2, 'MaxHeadSize', 0.5);

% Final rotated Vinf (with out-of-plane rotation)
quiver3(0,0,0, v_inf_vec(1), v_inf_vec(2), v_inf_vec(3), ...
        'LineWidth', 2, 'MaxHeadSize', 0.5);

% Reference axes
quiver3(0,0,0, 1,0,0,'k--');
quiver3(0,0,0, 0,1,0,'k--');
quiver3(0,0,0, 0,0,1,'k--');

xlabel('X (synodic)');
ylabel('Y (synodic)');
zlabel('Z (synodic)');

legend( ...
    'Original velocity direction', ...
    'After XY rotation (fpa clockwise)', ...
    'Final V_\infty vector (with out-of-plane)', ...
    'Location', 'best');

title('3D Visualization of V_\infty Rotations (Updated Signs)');

view(35,25);