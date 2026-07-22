function [SV_sync] = car2synodic(svcar, time, mu, theta_0)
% car2synodic  Convert an inertial Cartesian state to the CR3BP synodic frame.
% Primaries: Sun (primary), Earth (secondary).
%
% Inputs:
%   svcar   - inertial Cartesian state, [1x3] position or [1x6] state
%   time    - non-dimensional time (rotation angle of the synodic frame)
%   mu      - CR3BP mass parameter
%   theta_0 - initial rotation angle [rad] (defaults to 0 if empty)
%
% Outputs:
%   SV_sync - state in the synodic frame ([1x3] or [1x6] to match input)

%addpath(genpath('C:\Users\Rita\Dropbox\PHD\Matlab\CranRepo'))

if isempty(theta_0)
    theta_0 = 0;
end

%Rotation
theta = time + theta_0;  % time needs to be adimensionalised
A_t = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

xcar=svcar(1:3);

%Renormalize and change center
X_car_xFix = (inv(A_t)*xcar')';
X_S = (X_car_xFix - [mu 0 0]);
% X_S = (inv(A_t)*x')';

if length(svcar)==6
    
    vcar=svcar(4:6);
    V_car_xFix = (inv(A_t)*vcar')';
    V_S=(V_car_xFix - [0 mu 0]);
    V_S=[V_S(1)+X_S(2) V_S(2)-X_S(1) V_S(3)];

    SV_sync=[X_S V_S];
else
    SV_sync=X_S;
end


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