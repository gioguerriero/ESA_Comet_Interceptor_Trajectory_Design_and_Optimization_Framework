function vinf_rotated = vinf_rotation(synodic_moon, Vinf, fpa, out_of_plane)

    % vinf_rotation  Build a v-infinity vector by rotating the Moon's synodic
    % velocity direction in-plane (fpa) and out-of-plane, then scaling by Vinf.
    %
    % Inputs:
    %   synodic_moon - Moon synodic state [1x6]; its velocity sets the reference direction
    %   Vinf         - v-infinity magnitude to apply
    %   fpa          - in-plane angle [rad]; fpa > 0 rotates clockwise in the XY plane
    %   out_of_plane - out-of-plane angle [rad]; > 0 rotates toward +Z
    %
    % Outputs:
    %   vinf_rotated - 3x1 v-infinity vector in synodic coordinates

    %% Normalize initial direction
    v_direction = synodic_moon(4:6) ./ norm(synodic_moon(4:6));

    %% Clockwise in-plane rotation: Rz(+theta) is counterclockwise, so use Rz(-fpa)
    Rz_fpa = [ cos(-fpa)  -sin(-fpa)   0;
               sin(-fpa)   cos(-fpa)   0;
               0           0           1 ];

    v_rot_xy = Rz_fpa * v_direction;

    %% Out-of-plane rotation via Rodrigues' formula; axis cross(v,k) makes
    %  out_of_plane > 0 tilt the vector toward +Z
    k = [0; 0; 1];

    rot_axis = cross(v_rot_xy, k);
    rot_axis = rot_axis / norm(rot_axis);

    % Skew-symmetric matrix of the rotation axis
    K = [     0              -rot_axis(3)   rot_axis(2);
          rot_axis(3)            0         -rot_axis(1);
         -rot_axis(2)       rot_axis(1)         0        ];

    R_out = eye(3) ...
          + sin(out_of_plane)*K ...
          + (1 - cos(out_of_plane))*(K*K);

    v_rot_3D = R_out * v_rot_xy;

    %% Apply Vinf magnitude
    vinf_rotated = Vinf * v_rot_3D;

end