function info = verify_flyby_geometry(out_cr3bp, c_const)
% verify_flyby_geometry  Visualize and verify the lunar flyby geometry.
%
% Reconstructs the flyby-hyperbola periapsis state (same logic as
% write_python_inputs.m) and plots it in NON-DIMENSIONAL Moon-relative
% coordinates (lengths / Lstar, velocities / Vstar), showing:
%   - the Moon at natural size (non-dim rMoon radius)
%   - the periapsis point
%   - the relative trajectory (hyperbola, 2-body propagated about the Moon)
%   - the relative velocities vinf_in (incoming), vinf_out (outgoing) and v_periapsis
%   - the B-plane with axes T_hat / R_hat, the B vector, and the incoming
%     asymptote piercing it at the aim point.
%
% It also prints a numerical check: propagating the hyperbola toward the two
% asymptotes must reproduce vinf_in / vinf_out.
%
% Inputs:
%   out_cr3bp - variables_organizer output (flyby vinf_in/out, moon_state, epoch)
%   c_const   - constants struct
%
% Outputs:
%   info - struct with the reconstructed geometry and verification residuals
%
%   Details:
%     out_cr3bp - variables_organizer output (uses out.flyby.vinf_in/out,
%                 moon_state, epoch)
%     c_const   - struct with Lstar, Tstar, rMoon and muMoon (or G+mMoon)
%     info      - struct with r_p, h_flyby, vinf magnitudes, turn angle, B-plane,
%                 and the verification errors at the asymptotes.

%% ===== CONSTANTS / UNITS =================================================
Lstar = c_const.Lstar;            % [km]
Tstar = c_const.Tstar;            % [s]
Vstar = Lstar / Tstar;            % [km/s]
rMoon = c_const.rMoon;            % [km]

if isfield(c_const, 'muMoon')
    muMoon = c_const.muMoon;
else
    muMoon = c_const.G * c_const.mMoon;
end

vinf_in  = out_cr3bp.flyby.vinf_in(:);
vinf_out = out_cr3bp.flyby.vinf_out(:);

%% ===== HYPERBOLA GEOMETRY (same as flyby_periapsis_state) ==============
[B_vec, B_T, B_R, B_mag, r_p, S_hat, T_hat, R_hat] = ...
    bplane_from_vinf(vinf_in, vinf_out, muMoon);

vin_hat  = vinf_in  / norm(vinf_in);
vout_hat = vinf_out / norm(vinf_out);

h_hat   = cross(vinf_in, vinf_out);  h_hat = h_hat / norm(h_hat);
v_p_hat = (vin_hat + vout_hat);      v_p_hat = v_p_hat / norm(v_p_hat);
p_hat   = cross(v_p_hat, h_hat);     % periapsis position (outer side)

vinf = 0.5 * (norm(vinf_in) + norm(vinf_out));
v_p  = sqrt(vinf^2 + 2*muMoon / r_p);

r_peri_vec = r_p * p_hat;            % [km]  Moon-relative
v_peri_vec = v_p * v_p_hat;          % [km/s]

turn_angle = acos(max(-1,min(1, dot(vin_hat, vout_hat))));   % [rad]

%% ===== RELATIVE TRAJECTORY PROPAGATION (2-body about the Moon) ====
r_target = max(66100, 20*r_p);       % ~ SOI lunare [km]
t_guess  = 3 * r_target / vinf;      % [s] abbondante: l'evento ferma prima

opts = odeset('AbsTol',1e-9,'RelTol',1e-10, ...
              'Events', @(t,s) soi_event(t, s, r_target));

s0 = [r_peri_vec; v_peri_vec];

[~, S_fwd, ~, Sf_end, ~] = ode45(@(t,s) twoBody(t,s,muMoon), [0  t_guess], s0, opts);
[~, S_bwd, ~, Sb_end, ~] = ode45(@(t,s) twoBody(t,s,muMoon), [0 -t_guess], s0, opts);

% Full trajectory: incoming branch (reversed) + outgoing branch
traj_km = [flipud(S_bwd(:,1:3)); S_fwd(:,1:3)];   % N×3 [km]

%% ===== ASYMPTOTE CHECK ================================================
% At the branch ends the velocity should tend to vinf_out (fwd) and vinf_in (bwd)
if ~isempty(Sf_end), v_out_far = Sf_end(end,4:6).'; else, v_out_far = S_fwd(end,4:6).'; end
if ~isempty(Sb_end), v_in_far  = Sb_end(end,4:6).'; else, v_in_far  = S_bwd(end,4:6).'; end

ang_err_out = acosd(max(-1,min(1, dot(v_out_far,vinf_out)/(norm(v_out_far)*norm(vinf_out)))));
ang_err_in  = acosd(max(-1,min(1, dot(v_in_far, vinf_in) /(norm(v_in_far) *norm(vinf_in)))));

%% ===== NON-DIMENSIONALIZATION ============================================
ad        = 1/Lstar;                 % km   -> adim
traj_ad   = traj_km * ad;
rperi_ad  = r_peri_vec * ad;
rMoon_ad  = rMoon * ad;
Bvec_ad   = B_vec * ad;
r_tgt_ad  = r_target * ad;

% Velocity arrow length (visual scale, NOT to scale with positions)
L_arr = 0.45 * r_tgt_ad;
arr_in  = vin_hat  * L_arr;          % incoming direction
arr_out = vout_hat * L_arr;
arr_vp  = v_p_hat  * L_arr;

%% ===== PLOT =============================================================
figure('Color','w','Name','Flyby geometry (adim, Moon-relative)'); hold on; grid on;

% --- Moon (natural size) ---
[xs,ys,zs] = sphere(40);
surf(xs*rMoon_ad, ys*rMoon_ad, zs*rMoon_ad, ...
     'FaceColor',[0.6 0.6 0.62],'EdgeColor','none','FaceAlpha',0.95);

% --- Relative trajectory ---
plot3(traj_ad(:,1), traj_ad(:,2), traj_ad(:,3), 'b-', 'LineWidth',1.6);

% --- Periapsis ---
plot3(rperi_ad(1), rperi_ad(2), rperi_ad(3), 'r.', 'MarkerSize',26);

% --- Relative velocities (quiver, visual scale only) ---
% vinf_in at the incoming end, vinf_out at the outgoing end, v_p at periapsis
p_in  = traj_ad(1,:).';
p_out = traj_ad(end,:).';
quiver3(p_in(1),p_in(2),p_in(3),      arr_in(1),arr_in(2),arr_in(3), 0, ...
        'Color',[0 0.5 0],'LineWidth',2,'MaxHeadSize',0.6);
quiver3(p_out(1),p_out(2),p_out(3),   arr_out(1),arr_out(2),arr_out(3), 0, ...
        'Color',[0.85 0.3 0],'LineWidth',2,'MaxHeadSize',0.6);
quiver3(rperi_ad(1),rperi_ad(2),rperi_ad(3), arr_vp(1),arr_vp(2),arr_vp(3), 0, ...
        'Color',[0.5 0 0.6],'LineWidth',1.6,'MaxHeadSize',0.6);

% --- B-plane (perpendicular to S_hat, through the Moon centre) ---
hs = 1.3 * B_mag * ad;               % patch half-side [adim]
c1 = ( hs*T_hat + hs*R_hat) * 1;
c2 = ( hs*T_hat - hs*R_hat) * 1;
c3 = (-hs*T_hat - hs*R_hat) * 1;
c4 = (-hs*T_hat + hs*R_hat) * 1;
Bpatch = [c1, c2, c3, c4];           % 3×4
patch('XData',Bpatch(1,:),'YData',Bpatch(2,:),'ZData',Bpatch(3,:), ...
      'FaceColor',[0.2 0.4 0.9],'FaceAlpha',0.12,'EdgeColor',[0.2 0.4 0.9]);

% T_hat / R_hat axes in the B-plane
quiver3(0,0,0, hs*T_hat(1),hs*T_hat(2),hs*T_hat(3), 0, ...
        'Color',[0.2 0.4 0.9],'LineWidth',1.2,'LineStyle','-','MaxHeadSize',0.5);
quiver3(0,0,0, hs*R_hat(1),hs*R_hat(2),hs*R_hat(3), 0, ...
        'Color',[0.2 0.4 0.9],'LineWidth',1.2,'LineStyle','-','MaxHeadSize',0.5);
text(hs*T_hat(1),hs*T_hat(2),hs*T_hat(3),'  T','Color',[0.2 0.4 0.9]);
text(hs*R_hat(1),hs*R_hat(2),hs*R_hat(3),'  R','Color',[0.2 0.4 0.9]);

% B vector (centre -> aim point)
quiver3(0,0,0, Bvec_ad(1),Bvec_ad(2),Bvec_ad(3), 0, ...
        'Color','k','LineWidth',2,'MaxHeadSize',0.5);
text(Bvec_ad(1),Bvec_ad(2),Bvec_ad(3),'  B','Color','k','FontWeight','bold');

% Incoming asymptote: line || S_hat through the aim point (B_vec)
asy = [Bvec_ad - S_hat*r_tgt_ad, Bvec_ad + S_hat*r_tgt_ad];
plot3(asy(1,:),asy(2,:),asy(3,:),'--','Color',[0.4 0.4 0.4],'LineWidth',1);

axis equal; view(35,20);
xlabel('x (adim)'); ylabel('y (adim)'); zlabel('z (adim)');
title(sprintf('Lunar flyby - r_p = %.0f km (h = %.0f km),  v_\\infty = %.3f km/s,  turn = %.1f deg', ...
              r_p, r_p-rMoon, vinf, rad2deg(turn_angle)));
legend({'Moon','rel. trajectory','periapsis','v_{\infty,in}','v_{\infty,out}', ...
        'v_{periapsis}','B-plane','T','R','B','incoming asymptote'}, ...
        'Location','eastoutside');
hold off;

%% ===== PRINT SUMMARY / CHECK ======================================
fprintf('\n========== VERIFY FLYBY GEOMETRY ==========\n');
fprintf('  r_p              : %.2f km   (h = %.2f km)\n', r_p, r_p-rMoon);
fprintf('  |vinf_in|        : %.5f km/s\n', norm(vinf_in));
fprintf('  |vinf_out|       : %.5f km/s\n', norm(vinf_out));
fprintf('  v_periapsis      : %.5f km/s\n', v_p);
fprintf('  turn angle       : %.4f deg\n', rad2deg(turn_angle));
fprintf('  B magnitude      : %.2f km  (B_T=%.2f, B_R=%.2f)\n', B_mag, B_T, B_R);
fprintf('  --- asymptote check (should be ~0) ---\n');
fprintf('  err. dir. vinf_out : %.4e deg\n', ang_err_out);
fprintf('  err. dir. vinf_in  : %.4e deg\n', ang_err_in);
if r_p < rMoon
    fprintf(2,'  WARNING: periapsis BELOW the lunar surface!\n');
end
fprintf('===========================================\n\n');

%% ===== OUTPUT ===========================================================
info.r_p            = r_p;
info.h_flyby        = r_p - rMoon;
info.vinf_in_kms    = norm(vinf_in);
info.vinf_out_kms   = norm(vinf_out);
info.v_peri_kms     = v_p;
info.turn_angle_deg = rad2deg(turn_angle);
info.B_mag = B_mag;  info.B_T = B_T;  info.B_R = B_R;
info.ang_err_in_deg  = ang_err_in;
info.ang_err_out_deg = ang_err_out;
info.r_peri_vec_km   = r_peri_vec.';
info.v_peri_vec_kms  = v_peri_vec.';

end


%% =========================================================================
function [value, isterminal, direction] = soi_event(~, s, r_target)
% Ferma l'integrazione quando il raggio Moon-relative raggiunge r_target.
    value      = norm(s(1:3)) - r_target;
    isterminal = 1;
    direction  = 0;
end
