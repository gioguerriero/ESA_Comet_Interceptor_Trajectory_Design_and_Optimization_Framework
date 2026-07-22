function S_syn = sun_J2000_to_synodic(S_J2000, epochs, mu)
% sun_J2000_to_synodic  Convert Sun-centred ECLIPJ2000 states to the
% non-dimensional Sun-Earth synodic frame. Exact inverse of synodic2sun_J2000
% (includes the Ldot term and instantaneous Sun-Earth ephemerides per epoch).
%
% Inputs:
%   S_J2000 - Nx6 (or 1x6) Sun-centred ECLIPJ2000 states [km, km/s]
%   epochs  - Nx1 (or scalar) epochs [ET s]
%   mu      - Sun-Earth mass parameter
%
% Outputs:
%   S_syn   - Nx6 non-dimensional synodic states (Sun-Earth barycentre)

    epochs = epochs(:).';
    N      = numel(epochs);

    % Instantaneous Sun-Earth geometry at each epoch (rotating frame reference)
    ES       = cspice_spkezr('EARTH', epochs, 'ECLIPJ2000', 'NONE', 'SUN');
    rE       = ES(1:3,:);
    vE       = ES(4:6,:);

    % Length, angular rate and their time derivative used for scaling
    L_all    = sqrt(sum(rE.^2, 1));
    h_all    = cross(rE, vE);
    h_nrm    = sqrt(sum(h_all.^2, 1));
    om_all   = h_nrm ./ L_all.^2;
    V_all    = L_all .* om_all;
    Ldot_all = dot(rE, vE, 1) ./ L_all;

    % Synodic frame axes: e1 Sun->Earth, e3 orbit normal, e2 completes the triad
    e1 = rE ./ L_all;
    e3 = h_all ./ h_nrm;
    e2 = cross(e3, e1);

    S_J2000 = reshape(S_J2000, N, 6);
    S_syn   = zeros(N, 6);

    for i = 1:N
        R   = [e1(:,i), e2(:,i), e3(:,i)];
        r_J = S_J2000(i, 1:3).';
        v_J = S_J2000(i, 4:6).';

        % Rotate and scale position, then shift origin to the barycentre
        r_syn_adim_mu = R' * r_J / L_all(i);
        r_syn         = r_syn_adim_mu - [mu; 0; 0];

        % Velocity: remove frame rotation (Coriolis) and length-rate (Ldot) terms
        omega_vec = e3(:,i) * om_all(i);
        v_rot_km  = R' * (v_J - cross(omega_vec, r_J) - r_J * Ldot_all(i)/L_all(i));
        v_syn     = v_rot_km / V_all(i);

        S_syn(i,:) = [r_syn.' v_syn.'];
    end
end
