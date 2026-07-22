function S_syn = flyby_periapsis_state(vinf_in, vinf_out, moon_state, epoch, muMoon, mu, c_const)
% flyby_periapsis_state  Reconstruct the flyby-hyperbola periapsis state
% (non-dimensional synodic) from patched-conics v-infinity vectors.
%
% The patched-conics data only give pre/post states at the Moon centre; this
% builds a first approximation of the true periapsis state from that geometry:
%     r_p     from bplane_from_vinf
%     h_hat   = vinf_in x vinf_out normalised           (plane normal)
%     v_p_hat = (vinf_in_hat + vinf_out_hat) normalised (periapsis velocity, bisector)
%     p_hat   = v_p_hat x h_hat                          (periapsis position, outer side)
%     v_p     = sqrt(vinf^2 + 2*muMoon/r_p)
% The heliocentric Moon state is then added and converted to synodic.
%
% Inputs:
%   vinf_in, vinf_out - 3x1 Moon-relative v-infinity, ECLIPJ2000 [km/s]
%   moon_state        - Moon state w.r.t. Sun, ECLIPJ2000 [km, km/s]
%   epoch             - flyby epoch [ET s]
%   muMoon            - Moon gravitational parameter [km^3/s^2]
%   mu                - Sun-Earth mass parameter
%   c_const           - constants struct (uses rMoon for the altitude check)
%
% Outputs:
%   S_syn             - 1x6 periapsis state, non-dimensional synodic

    vinf_in  = vinf_in(:);
    vinf_out = vinf_out(:);

    % Periapsis radius from the B-plane geometry
    [~, ~, ~, ~, r_p] = bplane_from_vinf(vinf_in, vinf_out, muMoon);

    vin_hat  = vinf_in  / norm(vinf_in);
    vout_hat = vinf_out / norm(vinf_out);

    h_vec = cross(vinf_in, vinf_out);
    h_hat = h_vec / norm(h_vec);

    v_p_hat = (vin_hat + vout_hat);
    v_p_hat = v_p_hat / norm(v_p_hat);           % periapsis velocity direction

    p_hat = cross(v_p_hat, h_hat);               % periapsis position direction (p_hat x v_p_hat = h_hat)

    vinf = 0.5 * (norm(vinf_in) + norm(vinf_out));   % mean |v_inf| (equal in patched conics)
    v_p  = sqrt(vinf^2 + 2*muMoon / r_p);            % periapsis speed [km/s]

    r_rel = r_p * p_hat;                          % Moon-relative position [km]
    v_rel = v_p * v_p_hat;                         % Moon-relative velocity [km/s]

    % Periapsis altitude check
    if isfield(c_const, 'rMoon')
        fprintf('  [flyby] r_p = %.1f km  (altitudine = %.1f km)\n', ...
                r_p, r_p - c_const.rMoon);
        if r_p < c_const.rMoon
            warning('flyby_periapsis_state: periasse SOTTO la superficie lunare (r_p < rMoon).');
        end
    end

    % Heliocentric ECLIPJ2000 periapsis state
    S_J2000 = moon_state(:).' + [r_rel.' , v_rel.'];   % 1×6 [km, km/s]

    % Convert to non-dimensional synodic (exact inverse of synodic2sun_J2000)
    S_syn = sun_J2000_to_synodic(S_J2000, epoch, mu);
end
