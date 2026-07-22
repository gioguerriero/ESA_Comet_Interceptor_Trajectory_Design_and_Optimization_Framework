function ds = CR3BP_STM(t, s, mu, n)
% CR3BP_STM  CR3BP equations of motion augmented with the State Transition
% Matrix. RHS for ode integration of the 6 states plus the 6x6 STM (42 total).
%
% Inputs:
%   t  - integration time (unused, kept for the ode interface)
%   s  - augmented state [6x1 state; 36x1 flattened STM]
%   mu - CR3BP mass parameter of the system
%   n  - mean motion of the secondary (always 1 in the CR3BP; unused)
%
% Outputs:
%   ds - derivative of the augmented state [6x1 state deriv; 36x1 STM deriv]

    % Preallocate the augmented state-derivative vector
    ds = zeros(42,1);

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
    ds(4) = 2*vy + x - (1-mu)/d^3 * (x+mu) - mu*(x-(1-mu))/r^3;    ds(5) = -2*vx + y - (1-mu)*y/d^3 - mu*y/r^3;    ds(6) = -(1-mu)*z/d^3 - mu*z/r^3;
    % % Compute the partial-derivatives of the Psuedo-Potential function - TODO: FILL IN THE EQUATION
    % Uxx = 1 - mu * (1 / r^3 - 3 * (x - 1 + mu)^2 / r^5) - (1 - mu) * (1 / d^3 - 3 * (x + mu)^2 / d^5);
    % Uxy = - 3 * mu * (x - 1 + mu) * y / r^5 - 3 * (1 - mu) * (x + mu) * y / d^5;
    % Uxz = - 3 * mu * (x - 1 + mu) * z / r^5 - 3 * (1 - mu) * (x + mu) * z / d^5;
    % Uyx = Uxy;
    % Uyy = 1 - mu * (1/r^3 - 3 * y^2 / r^5) - (1 - mu) * (1/d^3 - 3 * y^2 / d^5);
    % Uyz = - 3 * y * z * (mu / r^5 + (1 - mu) / d^5);
    % Uzx = Uxz;
    % Uzy = Uyz;
    % Uzz = - mu * (1/r^3 - 3 * z^2 / r^5) - (1 - mu) * (1/d^3 - 3 * z^2 / d^5);
    % 
    % %test dai documenti di guardabasso -> ottengo lo stesso risultato
    % r1 = d;
    % r2 = r;
    % Uxx = -1 + (1 - mu)/r1^3 + mu/r2^3 - (3*(1 - mu)*(x + mu)^2)/r1^5 - (3*mu*(x - 1 + mu)^2)/r2^5;
    % Uyy = -1 + (1 - mu)/r1^3 + mu/r2^3 - (3*(1 - mu)*y^2)/r1^5 - (3*mu*y^2)/r2^5;
    % Uzz = (1 - mu)/r1^3 + mu/r2^3 - (3*(1 - mu)*z^2)/r1^5 - (3*mu*z^2)/r2^5;
    % Uxy = - (3*(1 - mu)*(x + mu)*y)/r1^5 - (3*mu*(x - 1 + mu)*y)/r2^5;
    % Uyx = Uxy;
    % Uxz = - (3*(1 - mu)*(x + mu)*z)/r1^5 - (3*mu*(x - 1 + mu)*z)/r2^5;
    % Uyz = - (3*(1 - mu)*y*z)/r1^5 - (3*mu*y*z)/r2^5;
    % Uzx = Uxz;
    % Uzy = Uyz;

    % Second partials of the pseudo-potential (from the referenced thesis)
    Uxx = 1 - (1 - mu)/d^3 - mu/r^3 + 3*(1 - mu)*(x + mu)^2/d^5 + 3*mu*(x - 1 + mu)^2/r^5;
    Uyy = 1 - (1 - mu)/d^3 - mu/r^3 + 3*(1 - mu)*y^2/d^5 + 3*mu*y^2/r^5;
    Uzz = - (1 - mu)/d^3 - mu/r^3 + 3*(1 - mu)*z^2/d^5 + 3*mu*z^2/r^5;
    Uxy = 3*(1 - mu)*(x + mu)*y/d^5 + 3*mu*(x - 1 + mu)*y/r^5;
    Uxz = 3*(1 - mu)*(x + mu)*z/d^5 + 3*mu*(x - 1 + mu)*z/r^5;
    Uyz = 3*(1 - mu)*y*z/d^5 + 3*mu*y*z/r^5;
    Uyx = Uxy;
    Uzx = Uxz;
    Uzy = Uyz;

    % Assemble the system (Jacobian) matrix A from the pseudo-potential partials
    A=[0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1;
       Uxx Uxy Uxz 0 2 0;
       Uyx Uyy Uyz -2 0 0;
       Uzx Uzy Uzz 0 0 0];

    % Rebuild the 6x6 STM from the flattened state
    phi = reshape(s(7:42), 6, 6);

    % STM derivative: Phi_dot = A * Phi
    phidot = A * phi;

    % Flatten the STM derivative back into the augmented state vector
    ds(7:42) = reshape(phidot, 36, 1);
end



