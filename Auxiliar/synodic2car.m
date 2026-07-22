function [svcar] = synodic2car(SV_sync, time, mu, theta_0)
% synodic2car  Convert a CR3BP synodic state to an inertial Cartesian state.
% Primaries: Sun (primary), Earth (secondary). Inverse of car2synodic.
%
% Inputs:
%   SV_sync - synodic state [1x6]
%   time    - non-dimensional time (rotation angle of the synodic frame)
%   mu      - CR3BP mass parameter
%   theta_0 - initial rotation angle [rad] (defaults to 0 if empty)
%
% Outputs:
%   svcar   - inertial Cartesian state [1x6]

%addpath(genpath('C:\Users\Rita\Dropbox\PHD\Matlab\CranRepo'))

if isempty(theta_0)
    theta_0 = 0;
end

%Rotation
theta = time + theta_0;  % time needs to be adimensionalised
A_t = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];


%     SV_sync=[X_S V_S];
X_S=SV_sync(1:3);
V_S=SV_sync(4:6);

X_car_xFix= (X_S + [mu 0 0]);
V_S=[V_S(1)-X_S(2) V_S(2)+X_S(1) V_S(3)];
V_car_xFix=(V_S + [0 mu 0]);

xcar=(A_t*X_car_xFix')';
vcar=(A_t*V_car_xFix')';

svcar=[xcar vcar];

end

%%
%     v_barycentre_normalized=[YE(i,4)-YE(i,2) YE(i,5)+YE(i,1) YE(i,6)];
%     r_barycentre_normalized=[YE(i,1:3)];
%     v_InvMan = [v_barycentre_normalized(1) v_barycentre_normalized(2)+mu v_barycentre_normalized(3)]*2*pi*AU/365/24/3600;
%     r_InvMan = [r_barycentre_normalized(1)+mu r_barycentre_normalized(2:3)]*AU;
%     theta=TE(i);
%     R=[cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1];
%     v1=(inv(R)*v_InvMan')';
%     r1=(inv(R)*r_InvMan')';