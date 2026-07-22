function d_hat = direction(state, alpha, phi)
% direction  Unit delta-v direction in the local velocity (TNH) frame,
% built from an in-plane and an out-of-plane steering angle.
%
% Inputs:
%   state - 6x1 state vector [r; v] (dimensional or non-dimensional)
%   alpha - in-plane steering angle [rad]
%   phi   - out-of-plane elevation angle [rad]
%
% Outputs:
%   d_hat - 3x1 unit direction vector for the delta-v

r = state(1:3);
v = state(4:6);

% Tangential direction (along velocity)
t_hat = v / norm(v);

% Angular momentum direction (out-of-plane)
h = cross(r,v);
h_hat = h / norm(h);

% In-plane normal direction
n_hat = cross(h_hat, t_hat);
n_hat = n_hat / norm(n_hat);

% Direction construction
d_hat = cos(phi)*cos(alpha)*t_hat + ...
        cos(phi)*sin(alpha)*n_hat + ...
        sin(phi)*h_hat;

% Safety normalization (numerical robustness)
d_hat = d_hat / norm(d_hat);

end