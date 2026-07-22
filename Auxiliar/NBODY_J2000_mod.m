function dS = NBODY_J2000_mod(t, S, initial_epoch, eps)
%NBODY_J2000_MOD  Homotopy dynamical model: CR3BP → full ephemeris.
%
%   dS = NBODY_J2000_mod(t, S, initial_epoch, eps)
%
%   State:   S = [x; y; z; vx; vy; vz]  in Sun-centred ECLIPJ2000 [km, km/s]
%   t             – integration time from sub-segment start [s]
%   initial_epoch – epoch ET [s] at t = 0
%   eps           – homotopy weight
%                     eps = 0 → Sun + Earth gravity (CR3BP-like, no Moon)
%                     eps = 1 → Sun + Earth + Moon  (full ephemeris)
%
%   Physical model (Sun-centred, non-barycentric frame):
%
%     a = a_Sun  +  a_Earth  +  eps * a_Moon
%
%   where each term includes the *indirect* (frame-correction) component:
%     a_body = mu_body * [ (r_body - r) / |r_body - r|^3
%                        - r_body      / |r_body|^3       ]
%
%   NOTE — why the rotating-frame approach was dropped:
%   The original implementation transformed to the Sun-Earth rotating frame,
%   evaluated the nondimensional CR3BP EOM, then converted back. That
%   approach had two bugs:
%     (1) Velocity normalization used v/L instead of v/(L*omega), making the
%         Coriolis contribution off by a factor of omega (~2e-7).
%     (2) The centrifugal term in the rotating→inertial back-conversion had
%         the wrong sign (+omega×(omega×r) instead of −).
%   Both errors are moot here: in an inertial frame the CR3BP acceleration IS
%   simply Sun + Earth gravity — no rotating frame algebra needed.

    % Standard gravitational parameters [km^3/s^2]
    % (using mu directly is more accurate than G*m)
    mu_sun   = 1.32712440018e11;
    mu_earth = 3.986004418e5;
    mu_moon  = 4.902800066e3;

    % Current epoch
    et = initial_epoch + t;

    % Spacecraft state
    r = S(1:3);
    v = S(4:6);

    % ---- Earth ephemeris (always needed) --------------------------------
    st_e  = cspice_spkezr('EARTH', et, 'ECLIPJ2000', 'NONE', 'SUN');
    r_e   = st_e(1:3);

    d_sc  = norm(r);
    d_e   = norm(r_e - r);
    d_e0  = norm(r_e);

    % Sun direct
    a_sun   = -mu_sun * r / d_sc^3;

    % Earth: direct + indirect
    a_earth = mu_earth * ((r_e - r) / d_e^3  -  r_e / d_e0^3);

    % ---- Moon (skip SPICE call when eps = 0) ----------------------------
    if eps > 0
        st_m = cspice_spkezr('MOON', et, 'ECLIPJ2000', 'NONE', 'SUN');
        r_m  = st_m(1:3);
        d_m  = norm(r_m - r);
        d_m0 = norm(r_m);
        a_moon = mu_moon * ((r_m - r) / d_m^3  -  r_m / d_m0^3);
    else
        a_moon = zeros(3,1);
    end

    % ---- Output ---------------------------------------------------------
    dS = [v; a_sun + a_earth + eps * a_moon];
end
