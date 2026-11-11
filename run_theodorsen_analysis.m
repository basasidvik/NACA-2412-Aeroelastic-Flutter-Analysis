% RUN_THEODORSEN_ANALYSIS
% Example runner for the separate circulatory Theodorsen extension.
%
% This script intentionally leaves run_flutter_analysis.m unchanged. Adjust
% the parameters below consistently with the quasi-steady example.

clear; clc; close all;

L_span = 1.5;
chord = 0.25;
AR = (2*L_span)^2/(2*L_span*chord);
e_s = 0.9;
EI = 450;
GJ = 180;
m_per_L = 0.35;
I_theta = 0.0018;
ea_loc_frac = 0.35;
a_ea = 2*ea_loc_frac - 1;
rho = 1.225;
U_vec = linspace(1,60,120);
n_el = 12;
n_modes_keep = [];       % retain all available oscillatory branches
alpha_rayleigh = 0.0;
beta_rayleigh = 2e-5;

aero = naca2412_data();
[Ms,Ks,node_y,dof_map,free_dofs] = beam_fem_matrices( ...
    L_span,n_el,EI,GJ,m_per_L,I_theta);
Cs = alpha_rayleigh*Ms + beta_rayleigh*Ks;

[U_vec,freq_Hz,growth_rate,converged,~,max_growth, ...
    dominant_eigenvalue,real_eigenvalues,spectrum_real_parts] = ...
    theodorsen_flutter_sweep( ...
    Ms,Ks,free_dofs,node_y,dof_map,chord,ea_loc_frac,a_ea, ...
    AR,e_s,aero,rho,U_vec,n_modes_keep,Cs,'jones');

first_unstable = find(max_growth >= 0,1,'first');
if isempty(first_unstable)
    fprintf('No complete-spectrum Theodorsen-extension instability in the swept range.\n');
else
    fprintf('First complete-spectrum Theodorsen-extension instability near %.2f m/s.\n', ...
            U_vec(first_unstable));
    fprintf('Dominant eigenvalue: %.6g%+.6gi 1/s\n', ...
            real(dominant_eigenvalue(first_unstable)), ...
            imag(dominant_eigenvalue(first_unstable)));
end
fprintf('Converged mode solutions: %.1f%%\n', ...
        100*nnz(converged)/numel(converged));

figure('Name','Theodorsen V-f');
plot(U_vec,freq_Hz.','LineWidth',1.3); grid on;
xlabel('Airspeed U [m/s]'); ylabel('Frequency [Hz]');
title('Theodorsen circulatory extension: V-f');

figure('Name','Theodorsen V-g');
plot(U_vec,growth_rate.','LineWidth',1.3); grid on; hold on;
yline(0,'k--');
xlabel('Airspeed U [m/s]'); ylabel('Growth rate Re(s) [1/s]');
title('Theodorsen circulatory extension: V-g');

figure('Name','Theodorsen real eigenvalues');
plot(U_vec,real_eigenvalues.','LineWidth',1.0); grid on; hold on;
yline(0,'k--');
xlabel('Airspeed U [m/s]'); ylabel('Real eigenvalue [1/s]');
title('Theodorsen extension: purely real/overdamped eigenvalues');

figure('Name','Theodorsen complete-spectrum real parts');
plot(U_vec,spectrum_real_parts.','LineWidth',0.8); grid on; hold on;
yline(0,'k--');
xlabel('Airspeed U [m/s]'); ylabel('Re(s) [1/s]');
title('Theodorsen extension: candidate real parts nearest zero');

if ~any(isfinite(real_eigenvalues(:)))
    fprintf(['No purely real Theodorsen roots were found. The separate ', ...
             'spectrum-real-parts plot shows damping/growth of all candidate roots.\n']);
end

%% Harmonic displacement comparison: quasi-steady versus Theodorsen/Sears
if ~isempty(first_unstable)
    U_compare = max(U_vec(1),0.8*U_vec(first_unstable));
else
    U_compare = 0.6*U_vec(end);
end
gust_amplitude = 1.0;
gust_frequency = 2.0;
amp_qs = nan(size(U_vec));
amp_unsteady = nan(size(U_vec));
for ku = 1:numel(U_vec)
    [~,~,amp_qs(ku)] = quasi_harmonic_gust_response( ...
        Ms,Ks,Cs,free_dofs,node_y,dof_map,chord,ea_loc_frac, ...
        AR,e_s,aero,U_vec(ku),rho,gust_amplitude,gust_frequency);
    [~,~,amp_unsteady(ku)] = sears_gust_response( ...
        Ms,Ks,Cs,free_dofs,node_y,dof_map,chord,ea_loc_frac,a_ea, ...
        AR,e_s,aero,U_vec(ku),rho,gust_amplitude,gust_frequency,'jones');
end

figure('Name','Gust displacement comparison');
plot(U_vec,amp_qs,'k--','LineWidth',1.4); hold on;
plot(U_vec,amp_unsteady,'b-','LineWidth',1.4); grid on;
xlabel('Airspeed U [m/s]'); ylabel('Tip displacement amplitude [m]');
legend('Quasi-steady','Theodorsen + Sears','Location','best');
title(sprintf('Harmonic gust displacement comparison at %.1f Hz',gust_frequency));

[q_qs,~,~] = quasi_harmonic_gust_response( ...
    Ms,Ks,Cs,free_dofs,node_y,dof_map,chord,ea_loc_frac, ...
    AR,e_s,aero,U_compare,rho,gust_amplitude,gust_frequency);
[q_unsteady,~,~] = sears_gust_response( ...
    Ms,Ks,Cs,free_dofs,node_y,dof_map,chord,ea_loc_frac,a_ea, ...
    AR,e_s,aero,U_compare,rho,gust_amplitude,gust_frequency,'jones');
t_compare = linspace(0,5,1001);
tip_w_compare = find(free_dofs == dof_map(end,1),1,'first');
omega_compare = 2*pi*gust_frequency;
z_qs = real(q_qs(tip_w_compare)*exp(1i*omega_compare*t_compare));
z_unsteady = real(q_unsteady(tip_w_compare)*exp(1i*omega_compare*t_compare));
figure('Name','Gust displacement time comparison');
plot(t_compare,z_qs,'k--','LineWidth',1.3); hold on;
plot(t_compare,z_unsteady,'b-','LineWidth',1.3); grid on;
xlabel('time [s]'); ylabel('Tip vertical displacement [m]');
legend('Quasi-steady','Theodorsen + Sears','Location','best');
title(sprintf('Harmonic gust response at U = %.2f m/s',U_compare));

