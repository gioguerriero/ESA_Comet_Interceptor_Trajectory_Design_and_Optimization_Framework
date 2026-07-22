function [S_J2000, scales] = synodic2sun_J2000(S_syn, epochs, mu)
%======================================================================
% SYNODIC2SUN_J2000
%
% Converts CR3BP rotating-frame adimensional states (Sun-Earth
% barycenter-centered, L-scaled, V-scaled) into Sun-centered J2000
% inertial dimensional states [km, km/s], using the INSTANTANEOUS
% Sun-Earth ephemeris at each epoch (so length and velocity scales
% are epoch-dependent).
%
% A single batch SPICE call is made for all epochs (efficient).
%
% INPUTS:
%   S_syn  - N×6 matrix of rotating-frame CR3BP states (adim)
%            [x y z vx vy vz], barycenter at origin, x-axis toward Earth
%   epochs - N×1 (or 1×N) ephemeris times [ET s]
%   mu     - Sun-Earth mass parameter (c.mu)
%
% OUTPUTS:
%   S_J2000 - N×6 Sun-centered J2000 states [km, km/s]
%   scales  - struct with per-epoch geometric info used for the
%             conversion. Useful downstream for ΔV conversions:
%               .L       1×N  instantaneous Sun-Earth distance [km]
%               .V       1×N  instantaneous velocity scale  [km/s]
%                             (= L * omega)
%               .omega   1×N  instantaneous angular rate    [rad/s]
%               .R       3×3×N rotation matrices (cols = e1,e2,e3)
%                             mapping rotating-frame vectors into
%                             Sun-centered J2000 (inertial).
%
% CONVENTION:
%   e1 = r_Earth / |r_Earth|        (Sun→Earth line)
%   e3 = (r × v)_Earth / |...|      (orbit normal)
%   e2 = e3 × e1                    (completes right-handed triad)
%
% POSITION:
%   r_syn_dim = ([x+mu; y; z]) * L_now
%   r_J2000   = R * r_syn_dim
%
% VELOCITY (with Coriolis correction for rotating→inertial):
%   v_syn_dim = [vx; vy; vz] * V_now
%   v_J2000   = R * v_syn_dim + omega_vec × r_J2000
%   omega_vec = e3 * omega_now
%======================================================================

    epochs = epochs(:)';            % force 1×N
    N = numel(epochs);

    if size(S_syn,1) ~= N
        error('synodic2sun_J2000:sizeMismatch', ...
              'S_syn must have one row per epoch (got %d rows, %d epochs).', ...
              size(S_syn,1), N);
    end

    % ---- Batch SPICE: Earth state w.r.t. Sun in J2000 ---------------
    ES  = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');  % 6×N
    rE  = ES(1:3,:);
    vE  = ES(4:6,:);

    L_all   = sqrt(sum(rE.^2, 1));             % 1×N [km]
    h_all   = cross(rE, vE);                   % 3×N
    h_nrm   = sqrt(sum(h_all.^2, 1));          % 1×N
    om_all  = h_nrm ./ L_all.^2;               % 1×N [rad/s]
    V_all   = L_all .* om_all;                 % 1×N [km/s]
    % Length-rate dL/dt of the Sun-Earth distance (velocity-scale correction)
    Ldot_all = dot(rE, vE, 1) ./ L_all;        % 1×N [km/s]

    e1 = rE ./ L_all;                          % 3×N
    e3 = h_all ./ h_nrm;                       % 3×N
    e2 = cross(e3, e1);                        % 3×N

    % ---- Per-epoch conversion ---------------------------------------
    S_J2000 = zeros(N, 6);
    Rmats   = zeros(3, 3, N);

    for i = 1:N
        R = [e1(:,i), e2(:,i), e3(:,i)];
        Rmats(:,:,i) = R;

        r_rot = [S_syn(i,1) + mu; S_syn(i,2); S_syn(i,3)] * L_all(i);
        v_rot = [S_syn(i,4);      S_syn(i,5); S_syn(i,6)] * V_all(i);

        r_J = R * r_rot;
        v_J = R * v_rot + cross(e3(:,i) * om_all(i), r_J);

        % Include the length-rate (Ldot) contribution in the inertial velocity
        r_body_dim = [S_syn(i,1) + mu; S_syn(i,2); S_syn(i,3)] * Ldot_all(i);
        v_J = R * v_rot + cross(e3(:,i) * om_all(i), r_J) + R * r_body_dim;

        S_J2000(i,:) = [r_J.' v_J.'];
    end

    scales.L     = L_all;
    scales.V     = V_all;
    scales.omega = om_all;
    scales.R     = Rmats;
end
