function dS = NBODY_J2000_full_ephe(t, S, initial_epoch, c)
% NBODY_J2000_full_ephe  Spacecraft acceleration under Sun + Earth + Moon
% point-mass gravity from SPICE ephemerides (Sun-centered ECLIPJ2000).
%
% Inputs:
%   t             - time elapsed since initial_epoch [s]
%   S             - Sun-centered spacecraft state [x y z vx vy vz] [km, km/s]
%   initial_epoch - reference epoch (SPICE ET) [s]
%   c             - constants struct (fields G, mEarth, mSun, mMoon)
%
% Outputs:
%   dS            - state derivative [vx vy vz ax ay az]

dS = zeros(6,1);

%% ------------------ Constants ------------------
G = c.G; % km^3/kg/s^2
m_earth = c.mEarth;
m_sun   = c.mSun;
m_moon = c.mMoon;

%% ------------------ Time ------------------
time = initial_epoch + t;

observer = 'SUN';
frame = 'ECLIPJ2000';

%% ------------------ Ephemerides ------------------
[posvel_moon, ~]  = cspice_spkezr('MOON',  time, frame, 'NONE', observer);
[posvel_earth, ~] = cspice_spkezr('EARTH', time, frame, 'NONE', observer);

r_moon = posvel_moon(1:3);
r_earth = posvel_earth(1:3);
v_earth = posvel_earth(4:6);

r_sc = S(1:3);
v_sc = S(4:6);

%% ------------------ N-body ephemeris acceleration ------------------

rsi = norm(r_sc);
rsm = norm(r_moon);
rse = norm(r_earth);

rim = norm(r_moon - r_sc);
rie = norm(r_earth - r_sc);

% a_ephem = ...
%     -G*m_sun * r_sc / rsi^3 + ...
%      G*m_earth * ((r_earth - r_sc)/rie^3 - r_earth/rse^3);

a_ephem = ...
    -G*m_sun * r_sc / rsi^3 + ...
     G*m_earth * ((r_earth - r_sc)/rie^3 - r_earth/rse^3) + ...
     G*m_moon  * ((r_moon  - r_sc)/rim^3 - r_moon /rsm^3);


%% ------------------ Output ------------------

dS(1:3) = v_sc;
dS(4:6) = a_ephem;

end