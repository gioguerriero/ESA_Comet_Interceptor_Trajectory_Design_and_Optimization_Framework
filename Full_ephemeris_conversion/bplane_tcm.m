function [dv_tcm, dv_mag, info] = bplane_tcm(state_tcm, epoch_tcm, epoch_flyby, ...
                                               flyby_body, BT_target, BR_target, ...
                                               mu_body, c, options)
% BPLANE_TCM  Compute a TCM that targets a desired B-plane point.
%
% The function propagates the spacecraft from the TCM epoch to the flyby
% epoch, computes B-plane parameters (Vallado Algorithm 78), and solves
% iteratively for the delta-v that drives (B_T, B_R) to the desired target
% using a Newton scheme with a numerical Jacobian.
%
% INPUTS
%   state_tcm   [6x1]  heliocentric ECLIPJ2000 state at TCM epoch [km, km/s]
%   epoch_tcm          SPICE ET at TCM point [s]
%   epoch_flyby        SPICE ET at flyby periapsis [s]
%   flyby_body         SPICE name of the flyby body (e.g. '67P', 'MARS') -> MOON
%   BT_target          desired B_T [km]
%   BR_target          desired B_R [km]
%   mu_body            GM of flyby body [km^3/s^2]
%   options (optional struct)
%       .dv_pert    finite-difference step [km/s]   default 1e-4
%       .max_iter   max Newton iterations           default 15
%       .tol        convergence threshold [km]      default 1e-2
%       .ode_opts   odeset options struct
%       .verbose    print iteration log (true/false) default true
%
% OUTPUTS
%   dv_tcm  [3x1]  TCM delta-v in heliocentric ECLIPJ2000 [km/s]
%   dv_mag         TCM magnitude [km/s]
%   info    struct
%       .converged      logical
%       .iter           iterations performed
%       .BT_initial     B_T before TCM [km]
%       .BR_initial     B_R before TCM [km]
%       .BT_final       B_T after TCM  [km]
%       .BR_final       B_R after TCM  [km]
%       .residual       [ΔB_T, ΔB_R] at convergence [km]
%       .J              last Jacobian (2x3)
%       .r_peri         periapsis distance from flyby body centre [km]
%       .h_flyby        flyby altitude above body surface [km]

%% --- Options -----------------------------------------------------------
if nargin < 9 || isempty(options), options = struct(); end
dv_pert  = getopt(options, 'dv_pert',  1e-4);
max_iter = getopt(options, 'max_iter', 15);
tol      = getopt(options, 'tol',      1e-4);
verbose  = getopt(options, 'verbose',  true);
ode_opts = getopt(options, 'ode_opts', odeset('AbsTol',1e-10,'RelTol',1e-10));

state_tcm = state_tcm(:);
tof = epoch_flyby - epoch_tcm;   % [s]

%% --- Nominal B-plane (Δv = 0) -----------------------------------------
[BT0, BR0] = compute_bplane(state_tcm, zeros(3,1), epoch_tcm, epoch_flyby, ...
                             flyby_body, mu_body, c, ode_opts);
info.BT_initial = BT0;
info.BR_initial = BR0;

if verbose
    fprintf('\n=== B-Plane TCM Targeting ===\n');
    fprintf('  TOF to flyby:  %.2f days\n', tof/86400);
    fprintf('  Target:  B_T = %.4f km,  B_R = %.4f km\n', BT_target, BR_target);
    fprintf('  Initial: B_T = %.4f km,  B_R = %.4f km\n', BT0, BR0);
    fprintf('  Error:   ΔB_T = %.4f km, ΔB_R = %.4f km\n', ...
            BT_target-BT0, BR_target-BR0);
    fprintf('-------------------------------\n');
end

%% --- Newton iterations -------------------------------------------------
dv = zeros(3,1);

for iter = 1:max_iter

    % Current B-plane at accumulated Δv
    [BT_cur, BR_cur] = compute_bplane(state_tcm, dv, epoch_tcm, epoch_flyby, ...
                                      flyby_body, mu_body, c, ode_opts);
    err = [BT_target - BT_cur;
           BR_target - BR_cur];

    if verbose
        fprintf('  iter %2d:  |err| = %.4e km   |Δv| = %.4e km/s\n', ...
                iter, norm(err), norm(dv));
    end

    if norm(err) < tol
        break
    end

    % Numerical Jacobian  J(2x3):  J_ij = ∂B_param_i / ∂Δv_j
    J = zeros(2, 3);
    for j = 1:3
        dv_p = dv;  dv_p(j) = dv_p(j) + dv_pert;
        dv_m = dv;  dv_m(j) = dv_m(j) - dv_pert;

        [BTp, BRp] = compute_bplane(state_tcm, dv_p, epoch_tcm, epoch_flyby, ...
                                    flyby_body, mu_body, c, ode_opts);
        [BTm, BRm] = compute_bplane(state_tcm, dv_m, epoch_tcm, epoch_flyby, ...
                                    flyby_body, mu_body, c, ode_opts);

        J(1,j) = (BTp - BTm) / (2*dv_pert);
        J(2,j) = (BRp - BRm) / (2*dv_pert);
    end

    % Minimum-norm correction: Δv_step = J^T (J J^T)^{-1} err
    dv_step = J' * ((J * J') \ err);
    dv = dv + dv_step;
end

%% --- Final evaluation --------------------------------------------------
[BT_fin, BR_fin, r_peri] = compute_bplane(state_tcm, dv, epoch_tcm, epoch_flyby, ...
                                           flyby_body, mu_body, c, ode_opts);
residual = [BT_target - BT_fin; BR_target - BR_fin];

% Flyby altitude above body surface
radii    = cspice_bodvrd(flyby_body, 'RADII', 3);
R_body   = radii(1);   % mean equatorial radius [km]
R_body  = c.rMoon;
h_flyby  = r_peri - R_body;

info.converged  = norm(residual) < tol;
info.iter       = iter;
info.BT_final   = BT_fin;
info.BR_final   = BR_fin;
info.residual   = residual;
info.J          = J;
info.r_peri     = r_peri;
info.h_flyby    = h_flyby;

dv_tcm = dv;
dv_mag = norm(dv);

if verbose
    fprintf('-------------------------------\n');
    fprintf('  Final:   B_T = %.4f km,  B_R = %.4f km\n', BT_fin, BR_fin);
    fprintf('  Residual: %.4e km\n', norm(residual));
    fprintf('  Periapsis: r_peri = %.4f km,  h_flyby = %.4f km\n', r_peri, h_flyby);
    fprintf('  TCM Δv = [%.6f  %.6f  %.6f] km/s\n', dv(1), dv(2), dv(3));
    fprintf('  |Δv|   = %.6f km/s  (%.4f m/s)\n', dv_mag, dv_mag*1e3);
    if info.converged
        fprintf('  STATUS: CONVERGED in %d iterations\n\n', iter);
    else
        fprintf('  STATUS: NOT CONVERGED (residual = %.4e km)\n\n', norm(residual));
    end
end

end % main function

%% =========================================================================
function [BT, BR, r_peri] = compute_bplane(state_tcm, dv, epoch_tcm, epoch_flyby, ...
                                            flyby_body, mu_body, c, ode_opts)
% Propagate state_tcm + dv from epoch_tcm to epoch_flyby, then compute
% B-plane parameters via Vallado Algorithm 78.

    % Apply TCM
    s0 = state_tcm;
    s0(4:6) = s0(4:6) + dv;

    % Propagate heliocentric ECLIPJ2000
    tof = epoch_flyby - epoch_tcm;
    [~, S] = ode45(@(t,s) NBODY_J2000_full_ephe(t, s, epoch_tcm, c), ...
                   [0, tof], s0, ode_opts);
    s_flyby = S(end,:)';   % heliocentric state at flyby epoch

    % Flyby body state from SPICE
    st_body = cspice_spkezr(flyby_body, epoch_flyby, 'ECLIPJ2000', 'NONE', 'SUN');
    r_body  = st_body(1:3);
    v_body  = st_body(4:6);

    % Relative (body-centred) state
    r_rel = s_flyby(1:3) - r_body;
    v_rel = s_flyby(4:6) - v_body;

    r_mag = norm(r_rel);
    v_mag = norm(v_rel);

    % Hyperbolic excess speed squared
    vinf2 = v_mag^2 - 2*mu_body/r_mag;
    if vinf2 <= 0
        warning('compute_bplane: v_inf^2 <= 0 at flyby epoch — spacecraft may be below escape speed relative to body.');
        vinf2 = abs(vinf2);   % fallback
    end

    % --- Vallado Algorithm 78 -------------------------------------------
    % Angular momentum unit vector
    h_vec = cross(r_rel, v_rel);
    h_hat = h_vec / norm(h_vec);

    % Eccentricity vector
    e_vec = ((v_mag^2 - mu_body/r_mag)*r_rel - dot(r_rel,v_rel)*v_rel) / mu_body;
    e_mag = norm(e_vec);
    e_hat = e_vec / e_mag;

    % Hyperbola semi-major and semi-minor axes (a < 0 for hyperbola)
    a = -mu_body / vinf2;
    b = -a * sqrt(e_mag^2 - 1);   % B = b (semi-minor axis of hyperbola)

    % Asymptote direction
    phi_s = acos(min(1, 1/e_mag));
    hxe   = cross(h_hat, e_hat);
    S_hat = e_hat*cos(phi_s) + hxe/norm(hxe)*sin(phi_s);

    % B-plane axes
    K_hat = [0; 0; 1];
    T_vec = cross(S_hat, K_hat);
    if norm(T_vec) < 1e-10
        K_hat = [0; 1; 0];   % fallback if S parallel to Z
        T_vec = cross(S_hat, K_hat);
    end
    T_hat = T_vec / norm(T_vec);
    R_hat = cross(S_hat, T_hat);

    % B vector
    B_hat = cross(S_hat, h_hat);
    B_vec = b * B_hat;

    BT = dot(B_vec, T_hat);
    BR = dot(B_vec, R_hat);

    % Periapsis distance: r_p = a*(1 - e)  [a < 0, e > 1 → r_p > 0]
    r_peri = a * (1 - e_mag);

end % compute_bplane

%% =========================================================================
function val = getopt(s, field, default)
% getopt  Return s.(field) if present, otherwise the supplied default.
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
