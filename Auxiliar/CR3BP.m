function ds = CR3BP(t, s, mu)
% CR3BP  Equations of motion of the Circular Restricted Three-Body Problem.
% Right-hand side for ode45 in the rotating (synodic), non-dimensional frame.
%
% Inputs:
%   t  - integration time (unused, kept for the ode45 interface)
%   s  - state [x y z vx vy vz] in synodic non-dimensional units
%   mu - CR3BP mass parameter of the system
%
% Outputs:
%   ds - state derivative [vx vy vz ax ay az]
    
    % Preallocate the state-derivative vector
    ds = zeros(6,1);

    % Unpack position and velocity components
    x = s(1);
    y = s(2);
    z = s(3);
    vx = s(4);
    vy = s(5);
    vz = s(6);
    
    % Compute the distance from the smaller primary
    r = norm(s(1:3)-[1-mu;0;0]);    
    % Compute the distance from the larger primary
    d = norm(s(1:3)+[mu;0;0]);    
    % Assign velocities along x, y and z axes of the CR3BP rotating-frame
    ds(1) = vx;
    ds(2) = vy;
    ds(3) = vz;
    
    % Assign accelerations along x, y and z axes of the CR3BP rotating-frame,
    % using equations of motion for CR3BP model
    ds(4) = 2*vy + x - (1-mu)/d^3 * (x+mu) - mu*(x-(1-mu))/r^3;    ds(5) = -2*vx + y - (1-mu)*y/d^3 - mu*y/r^3;    ds(6) = -(1-mu)*z/d^3 - mu*z/r^3;end


