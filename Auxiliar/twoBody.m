function dsdt = twoBody(t, S, muReal)
% twoBody  Two-body (Keplerian) equations of motion for ode integration.
%
% Inputs:
%   t       - integration time (unused, kept for the ode interface)
%   S       - state [x y z vx vy vz] in an inertial frame
%   muReal  - gravitational parameter of the central body [km^3/s^2]
%
% Outputs:
%   dsdt    - state derivative [vx vy vz ax ay az]
    r = sqrt(S(1)^2 + S(2)^2 + S(3)^2);
    dsdt = [S(4); S(5); S(6); -muReal*S(1)/r^3; -muReal*S(2)/r^3; -muReal*S(3)/r^3];
end



