% RUN_FLUTTER_ANALYSIS
% Top-level script: coupled structural (FEM beam) + quasi-steady strip-theory
% aerodynamic flutter / vibration analysis of a finite-span NACA 2412 wing.
%
% Mirrors aeroelastic_naca2412.tex section by section. Edit the PARAMETERS
% block below for a specific wing; everything downstream is generic.

clear; clc; close all;

%% ------------------------- PARAMETERS -----------------------------------

% --- Geometry ---
L_span   = 1.5;        % semi-span [m]  (root to tip)
chord    = 0.25;       % constant chord [m] (or use a function handle @(y) ... )
AR       = (2*L_span)^2 / (2*L_span*chord);   % aspect ratio of full wing (b^2/S)
e_s      = 0.9;         % Oswald span efficiency factor

% --- Structural properties (uniform along span; replace with handles @(y)... for taper) ---
EI       = 450;         % bending stiffness [N*m^2]
GJ       = 180;         % torsional stiffness [N*m^2]
m_per_L  = 0.35;        % mass per unit length [kg/m]
I_theta  = 0.0018;      % mass moment of inertia per unit length about EA [kg*m]

% --- Elastic axis / aerodynamic center geometry ---
ea_loc_frac = 0.35;     % elastic axis location, fraction of chord from LE
a_ea        = 2*ea_loc_frac - 1;   % EA location from mid-chord, used in pitch-rate term

% --- Flight condition ---
rho      = 1.225;       % air density [kg/m^3]
U_vec    = linspace(1, 60, 200);   % airspeed sweep [m/s]

% --- Discretization ---
n_el = 12;              % number of beam elements along the span

% --- Modal tracking ---
n_modes_keep = [];       % [] = retain all oscillatory modes for complete search

% --- Structural damping (Rayleigh damping, Cs = alpha*Ms + beta*Ks) ---
% A small nonzero beta is used by default. The baseline aerodynamic model
% (see aero_strip_matrices.m) has no torsional aerodynamic damping term, so
% with EXACTLY zero structural damping the torsion-dominated mode can show a
% very weak instability even at low airspeed -- an expected consequence of
% that omission, not a bug. This small beta stands in for it; replace with
% measured/estimated structural damping, or add a proper torsional
% aerodynamic damping term, for quantitative work.
alpha_rayleigh = 0.0;
beta_rayleigh  = 2e-5;

%% ------------------------- BUILD MODELS ---------------------------------

aero = naca2412_data();

[Ms, Ks, node_y, dof_map, free_dofs] = beam_fem_matrices(L_span, n_el, EI, GJ, ...
                                                          m_per_L, I_theta);

Cs = alpha_rayleigh*Ms + beta_rayleigh*Ks;   % full-size; reduced inside flutter_sweep

%% ------------------------- FLUTTER / VIBRATION SWEEP --------------------

[U_vec, freq_Hz, damping, ~, divergence_indicator, max_growth, ...
                       dominant_eigenvalue, real_eigenvalues] = ...
    flutter_sweep(Ms, Ks, free_dofs, node_y, dof_map, ...
                                chord, ea_loc_frac, a_ea, AR, e_s, aero, rho, ...
                                U_vec, n_modes_keep, Cs);

%% ------------------------- FLUTTER SPEED DETECTION -----------------------

% Use the complete state-space spectrum, including real eigenvalues. This
% avoids missing an instability simply because its oscillatory pair vanished.
unstable_idx = find(max_growth >= 0, 1, 'first');
U_flutter = NaN;
if ~isempty(unstable_idx)
    if unstable_idx == 1
        U_flutter = U_vec(1);
    else
        U1 = U_vec(unstable_idx-1); g1 = max_growth(unstable_idx-1);
        U2 = U_vec(unstable_idx);   g2 = max_growth(unstable_idx);
        U_flutter = U1 + (0-g1)*(U2-U1)/(g2-g1);
    end
end

%% ------------------------- DIVERGENCE DETECTION --------------------------

% A sign change in det(Kr+Sa)/det(Kr) brackets static divergence. This is
% checked separately from the dynamic Re(s)=0 flutter crossing.
div_idx = find(divergence_indicator(1:end-1).*divergence_indicator(2:end) <= 0, 1, 'first');
U_divergence = NaN;
if ~isempty(div_idx)
    d1 = divergence_indicator(div_idx);
    d2 = divergence_indicator(div_idx+1);
    U1 = U_vec(div_idx); U2 = U_vec(div_idx+1);
    U_divergence = U1 + (0-d1)*(U2-U1)/(d2-d1);
    fprintf('Estimated static-divergence speed: U_D = %.2f m/s\n', U_divergence);
else
    fprintf('No static-divergence determinant crossing found in the swept speed range.\n');
end

figure('Name','Static divergence indicator');
plot(U_vec, divergence_indicator, 'LineWidth', 1.5); grid on; hold on;
yline(0, 'k--');
if ~isnan(U_divergence)
    xline(U_divergence, 'r--', sprintf('U_D = %.1f m/s', U_divergence));
end
xlabel('Airspeed U [m/s]');
ylabel('det(K_r + S_{a,r}) / det(K_r)');
title('Static-divergence indicator');

if ~isnan(U_flutter)
    fprintf('First complete-spectrum instability: U_F = %.2f m/s\n', U_flutter);
    fprintf('Dominant eigenvalue at first unstable grid point: %.6g%+.6gi 1/s\n', ...
            real(dominant_eigenvalue(unstable_idx)), ...
            imag(dominant_eigenvalue(unstable_idx)));
else
    fprintf('No flutter crossing found in the swept speed range [%.1f, %.1f] m/s.\n', ...
            U_vec(1), U_vec(end));
end

%% ------------------------- PLOTS -----------------------------------------

figure('Name','V-f diagram');
plot(U_vec, freq_Hz.', 'LineWidth', 1.5); grid on;
xlabel('Airspeed U [m/s]'); ylabel('Frequency [Hz]');
title('V-f diagram: modal frequency vs. airspeed');
legend(arrayfun(@(k) sprintf('mode %d',k), 1:size(freq_Hz,1), 'UniformOutput', false));

figure('Name','V-g diagram');
plot(U_vec, damping.', 'LineWidth', 1.5); grid on; hold on;
yline(0, 'k--');
if ~isnan(U_flutter)
    xline(U_flutter, 'r--', sprintf('U_F = %.1f m/s', U_flutter));
end
xlabel('Airspeed U [m/s]'); ylabel('Damping, Re(s) [1/s]');
title('V-g diagram: modal damping vs. airspeed');
legend(arrayfun(@(k) sprintf('mode %d',k), 1:size(damping,1), 'UniformOutput', false));

figure('Name','Real eigenvalues');
plot(U_vec, real_eigenvalues.', 'LineWidth', 1.0); grid on; hold on;
yline(0, 'k--');
if ~isnan(U_flutter)
    xline(U_flutter, 'r--', sprintf('U_F = %.1f m/s', U_flutter));
end
xlabel('Airspeed U [m/s]'); ylabel('Real eigenvalue [1/s]');
title('Real/overdamped eigenvalues');

% Report the actual variation in the stored histories. This distinguishes a
% genuinely speed-independent result from a plotting-orientation problem.
freq_span = zeros(size(freq_Hz,1),1);
growth_span = zeros(size(damping,1),1);
for kdiag = 1:size(freq_Hz,1)
    fk = freq_Hz(kdiag,isfinite(freq_Hz(kdiag,:)));
    gk = damping(kdiag,isfinite(damping(kdiag,:)));
    if ~isempty(fk); freq_span(kdiag) = max(fk)-min(fk); end
    if ~isempty(gk); growth_span(kdiag) = max(gk)-min(gk); end
end
fprintf('Largest oscillatory frequency variation: %.6g Hz\n', max(freq_span));
fprintf('Largest tracked growth-rate variation: %.6g 1/s\n', max(growth_span));
fprintf('Dominant growth-rate range: %.6g to %.6g 1/s\n', ...
        min(max_growth), max(max_growth));

%% ------------------------- TIP VERTICAL VIBRATION AT A GIVEN SPEED -------
% Optional: time-domain forced response to an initial disturbance, well
% below flutter speed, to see the finite-span wing's vertical (plunge-
% dominated) vibration decay/growth directly.

U_check = 0.5 * (U_flutter * ~isnan(U_flutter) + U_vec(end) * isnan(U_flutter));
if isnan(U_flutter); U_check = 0.6*U_vec(end); end

[Da_chk, Sa_chk, ~] = aero_strip_matrices(node_y, dof_map, chord, ea_loc_frac, a_ea, ...
                                           AR, e_s, aero, U_check, rho);
Mr = Ms(free_dofs, free_dofs);
Kr = Ks(free_dofs, free_dofs) + Sa_chk(free_dofs, free_dofs);
Cr = Cs(free_dofs, free_dofs) + Da_chk(free_dofs, free_dofs);

n = numel(free_dofs);
odefun = @(t,z) [z(n+1:end); Mr\(-Cr*z(n+1:end) - Kr*z(1:n))];

z0 = zeros(2*n,1);
% perturb the tip transverse DOF (last node's 'w' DOF, first entry in its block)
tip_w_local = find(free_dofs == dof_map(end,1));
z0(tip_w_local) = 0.01;   % 1 cm initial tip displacement

[t, Z] = ode45(odefun, [0 2], z0);

figure('Name','Tip vertical vibration time history');
plot(t, Z(:, tip_w_local), 'LineWidth', 1.5); grid on;
xlabel('time [s]'); ylabel('Tip vertical displacement [m]');
title(sprintf('Tip plunge response at U = %.1f m/s', U_check));

%% ------------------------- FORCED GUST RESPONSE ---------------------------
% Prescribed sinusoidal vertical gust. The resulting lift and pitching
% moment are assembled as an external generalized force and integrated in
% time, rather than using an initial displacement only.
gust_amplitude = 1.0;       % vertical gust velocity amplitude [m/s]
gust_frequency = 2.0;       % gust frequency [Hz]
gust_tspan = linspace(0, 5, 1001);
gust_velocity = @(t) gust_amplitude*sin(2*pi*gust_frequency*t);

[tg, Zg] = gust_response(Ms, Ks, Cs, free_dofs, node_y, dof_map, ...
                         chord, ea_loc_frac, AR, e_s, aero, U_check, rho, ...
                         gust_tspan, gust_velocity, []);
tip_w_local = find(free_dofs == dof_map(end,1));
figure('Name','Forced gust response');
plot(tg, Zg(:,tip_w_local), 'LineWidth', 1.5); grid on;
xlabel('time [s]'); ylabel('Tip vertical displacement [m]');
title(sprintf('Forced tip response: %.1f m/s sinusoidal gust at U = %.1f m/s', ...
             gust_frequency, U_check));
