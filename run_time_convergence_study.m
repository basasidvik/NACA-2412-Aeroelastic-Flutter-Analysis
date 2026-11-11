% RUN_TIME_CONVERGENCE_STUDY
% Compare 4-, 8-, and 12-mode truncated FEM gust responses.
%
% The modal bases are calculated from the full FEM matrices. Each reduced
% system is integrated with ode15s, and the resulting tip displacement
% histories are plotted side-by-side and overlaid for comparison.

clear; clc; close all;

%% Model parameters
L_span = 1.5;
chord = 0.25;
AR = (2*L_span)^2/(2*L_span*chord);
e_s = 0.9;
EI = 450;
GJ = 180;
m_per_L = 0.35;
I_theta = 0.0018;
ea_loc_frac = 0.35;
a_ea = 2*ea_loc_frac-1;
rho = 1.225;
n_el = 12;
U_check = 12.0;                 % below the default quasi-steady flutter speed
gust_amplitude = 1.0;           % vertical gust velocity [m/s]
gust_frequency = 2.0;           % gust frequency [Hz]
t_span = linspace(0,5,1001);
mode_counts = [4 8 12];

alpha_rayleigh = 0.0;
beta_rayleigh = 2e-5;

%% Build full FEM model
aero = naca2412_data();
[Ms,Ks,node_y,dof_map,free_dofs] = beam_fem_matrices( ...
    L_span,n_el,EI,GJ,m_per_L,I_theta);
Cs = alpha_rayleigh*Ms + beta_rayleigh*Ks;
gust_velocity = @(t) gust_amplitude*sin(2*pi*gust_frequency*t);

tip_histories = cell(numel(mode_counts),1);
metrics = cell(numel(mode_counts),1);
modal_frequencies = cell(numel(mode_counts),1);

for im = 1:numel(mode_counts)
    nm = mode_counts(im);
    [t,tip_histories{im},~,~,modal_frequencies{im},metrics{im}] = ...
        modal_gust_response(Ms,Ks,Cs,free_dofs,node_y,dof_map, ...
                            chord,ea_loc_frac,a_ea,AR,e_s,aero, ...
                            U_check,rho,t_span,gust_velocity,nm,[]);
    fprintf('%2d modes: max |tip| = %.6g m, RMS = %.6g m\n', ...
            nm,metrics{im}.max_abs_tip,metrics{im}.rms_tip);
end

%% Side-by-side responses
figure('Name','Modal convergence side-by-side');
for im = 1:numel(mode_counts)
    subplot(numel(mode_counts),1,im);
    plot(t,tip_histories{im},'LineWidth',1.2); grid on;
    ylabel(sprintf('%d modes: w [m]',mode_counts(im)));
    title(sprintf('ode15s response using %d retained modes',mode_counts(im)));
end
xlabel('time [s]');

%% Overlay and convergence error
figure('Name','Modal convergence overlay');
plot(t,tip_histories{1},'--','LineWidth',1.1); hold on;
plot(t,tip_histories{2},'-','LineWidth',1.1);
plot(t,tip_histories{3},'-.','LineWidth',1.1); grid on;
xlabel('time [s]'); ylabel('Tip vertical displacement [m]');
legend('4 modes','8 modes','12 modes','Location','best');
title(sprintf('Modal convergence at U = %.1f m/s',U_check));

reference = tip_histories{end};
relative_error = zeros(numel(mode_counts)-1,1);
for im = 1:numel(mode_counts)-1
    relative_error(im) = norm(tip_histories{im}-reference) / ...
                         max(norm(reference),eps);
end
figure('Name','Modal convergence error');
bar(mode_counts(1:end-1),relative_error); grid on;
xlabel('Retained modes'); ylabel('Relative L2 error versus 12 modes');
title('Time-response modal convergence');
