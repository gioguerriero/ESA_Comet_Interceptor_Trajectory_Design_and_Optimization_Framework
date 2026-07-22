function [B_vec, B_T, B_R, B_mag, r_p, S_hat, T_hat, R_hat] = bplane_from_vinf(vinf_in, vinf_out, mu)
% B-Plane parameters from incoming and outgoing v-infinity vectors.
% Implements Algorithm 79 (Vallado, "Fundamentals of Astrodynamics and
% Applications"), valid for gravity-assist flybys.
%
% Inputs:
%   vinf_in   - [3x1] incoming v-infinity vector (km/s), body-relative frame
%   vinf_out  - [3x1] outgoing v-infinity vector (km/s), body-relative frame
%   mu        - gravitational parameter of the flyby body (km^3/s^2)
%
% Outputs:
%   B_vec  - [3x1] B-vector (km)
%   B_T    - B-plane T component (km)
%   B_R    - B-plane R component (km)
%   B_mag  - magnitude of B (km)
%   r_p    - periapsis radius (km)
%   S_hat  - unit vector along incoming asymptote
%   T_hat  - T unit vector (in ecliptic plane, perpendicular to S)
%   R_hat  - R unit vector (normal to S and T)

vinf_in  = vinf_in(:);
vinf_out = vinf_out(:);

vinf_in_mag  = norm(vinf_in);
vinf_out_mag = norm(vinf_out);

% --- Reference frame construction (Algorithm 79) ---
S_hat = vinf_in / vinf_in_mag;

h_hat = cross(vinf_in, vinf_out);
h_hat = h_hat / norm(h_hat);

B_hat = cross(S_hat, h_hat);           % unit B direction

K_hat = [0; 0; 1];                     % ecliptic Z axis

T_hat = cross(S_hat, K_hat);
T_hat = T_hat / norm(T_hat);

R_hat = cross(S_hat, T_hat);

% --- Turning angle and periapsis radius ---
cos_phi = dot(vinf_in, vinf_out) / (vinf_in_mag * vinf_out_mag);
phi     = acos(max(-1, min(1, cos_phi)));

r_p = (mu / vinf_in_mag^2) * (1 / cos((pi - phi) / 2) - 1);

% --- B magnitude and vector ---
B_mag = (mu / vinf_in_mag^2) * sqrt((1 + vinf_in_mag^2 * r_p / mu)^2 - 1);

B_vec = B_mag * B_hat;

B_T = dot(B_vec, T_hat);
B_R = dot(B_vec, R_hat);

end
