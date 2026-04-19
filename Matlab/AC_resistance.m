%% AC Resistance Comparison: Copper vs Aluminium
% Compares AC resistance as a function of frequency using the skin effect model
% Based on equations from Pyrhonen 2014 Design of Machines

clear all; close all; clc;

%% Physical Constants
mu_0 = 4 * pi * 1e-7;  % Permeability of free space (H/m)

%% Material Parameters
rho_copper = 1.724e-8;     % Copper resistivity (Ω·m)

% Fallback values if JSON not found
h_slot = 21e-3;        % Stator slot height (m)
h_bar = 21e-3;         % Conductor bar height (m)
w_bar = 2.5e-3;        % Conductor bar width (m)
l_stack = 250e-3;      % Stack length (m)
rho_aluminum = 2.65e-8; % Aluminium resistivity (Ω·m)

%% Calculate cross-sectional area and DC resistance
A_bar = w_bar * h_bar;  % Cross-sectional area of conductor bar

% DC Resistance: R_DC = ρ * L / A
R_DC_copper = rho_copper * l_stack / A_bar;
R_DC_aluminum = rho_aluminum * l_stack / A_bar;

fprintf('DC Resistance (Copper): %.6e Ω\n', R_DC_copper);
fprintf('DC Resistance (Aluminium): %.6e Ω\n', R_DC_aluminum);

%% Frequency range (0 to 500 Hz)
f_start = 1e-3;  % Start from 1 mHz to avoid singularities
f_end = 500;
n_points = 5000;
f = linspace(f_start, f_end, n_points);

%% Calculate AC resistance factor k_r for both materials
% k_r = ξ * (sinh(2ξ) + sin(2ξ)) / (cosh(2ξ) - cos(2ξ))
% ξ = h_slot * sqrt(π * μ₀ * f / ρ)

% For copper
xi_copper = h_slot * sqrt(pi * mu_0 * f / rho_copper);
k_r_copper = xi_copper .* (sinh(2*xi_copper) + sin(2*xi_copper)) ./ ...
             (cosh(2*xi_copper) - cos(2*xi_copper));

% For aluminium
xi_aluminum = h_slot * sqrt(pi * mu_0 * f / rho_aluminum);
k_r_aluminum = xi_aluminum .* (sinh(2*xi_aluminum) + sin(2*xi_aluminum)) ./ ...
               (cosh(2*xi_aluminum) - cos(2*xi_aluminum));

%% Calculate AC resistance: R_AC = k_r * R_DC
R_AC_copper = k_r_copper * R_DC_copper;
R_AC_aluminum = k_r_aluminum * R_DC_aluminum;

%% Create plots
figure1 = figure;
matlab_blue = [0 0.4470 0.7410];
matlab_orange = [0.8500 0.3250 0.0980];

% Left subplot: AC Resistance Factor k_r vs Frequency
subplot(1, 2, 1);
hold on;
h1 = plot(f, k_r_copper, '-', 'Color', matlab_blue, 'LineWidth', 1.5, 'DisplayName', 'Copper (\rho = 1.724\times10^{-8} \Omega\cdotm)');
h2 = plot(f, k_r_aluminum, '-', 'Color', matlab_orange, 'LineWidth', 1.5, 'DisplayName', 'Aluminium (\rho = 2.65\times10^{-8} \Omega\cdotm)');
xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('AC Resistance Factor, k_r');
grid;
set(gca, 'Position', [0.08 0.32 0.34 0.61]);
hold off;

% Right subplot: AC Resistance vs Frequency
subplot(1, 2, 2);
hold on;
plot(f, R_AC_copper*1e6, '-', 'Color', matlab_blue, 'LineWidth', 1.5, 'DisplayName', 'Copper (\rho = 1.724\times10^{-8} \Omega\cdotm)');
plot(f, R_AC_aluminum*1e6, '-', 'Color', matlab_orange, 'LineWidth', 1.5, 'DisplayName', 'Aluminium (\rho = 2.65\times10^{-8} \Omega\cdotm)');
xlim([0 500]);
xlabel('Frequency (Hz)');
ylabel('AC Resistance (\mu\Omega)');
grid;
set(gca, 'Position', [0.60 0.32 0.34 0.61]);
hold off;

% Figure-level shared legend below both subplots (1x2 layout)
ax_lgd = axes('Parent', figure1, 'Position', [0 0 1 1], 'Visible', 'off');
hold(ax_lgd, 'on');
hl1 = plot(ax_lgd, nan, nan, '-', 'Color', matlab_blue, 'LineWidth', 1.5);
hl2 = plot(ax_lgd, nan, nan, '-', 'Color', matlab_orange, 'LineWidth', 1.5);
lgd = legend(ax_lgd, [hl1 hl2], {'Copper (\rho = 1.724\times10^{-8} \Omega\cdotm)', 'Aluminium (\rho = 2.65\times10^{-8} \Omega\cdotm)'});
set(lgd, 'Orientation', 'horizontal', 'Location', 'southoutside', 'NumColumns', 2, 'Box', 'on');
hold(ax_lgd, 'off');

% Apply figure formatting style (double width, same height)
set(findall(figure1, 'Type', 'text'), 'FontSize', 9, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
set(findall(figure1, '-property', 'FontSize'), 'FontSize', 9, 'FontWeight', 'normal', 'FontName', 'Times New Roman');
set(figure1, 'Units', 'centimeters', 'PaperUnits', 'centimeters', 'Position', [0 0 17 5.2], 'PaperSize', [17 5.2]);

%% Print summary information
fprintf('\n=== Summary ===\n');
fprintf('Slot height (h_slot): %.2f mm\n', h_slot*1e3);
fprintf('Conductor bar dimensions: %.2f mm × %.2f mm\n', w_bar*1e3, h_bar*1e3);
fprintf('Conductor bar cross-section: %.2f mm²\n', A_bar*1e6);
fprintf('Stack length (l_stack): %.2f mm\n', l_stack*1e3);

fprintf('\nDC Resistance (Copper): %.4f μΩ\n', R_DC_copper*1e6);
fprintf('DC Resistance (Aluminium): %.4f μΩ\n', R_DC_aluminum*1e6);
fprintf('Ratio (Al/Cu): %.3f\n', R_DC_aluminum/R_DC_copper);

fprintf('\nAt 500 Hz:\n');
fprintf('  k_r (Copper): %.4f\n', k_r_copper(end));
fprintf('  k_r (Aluminium): %.4f\n', k_r_aluminum(end));
fprintf('  R_AC (Copper): %.4f μΩ\n', R_AC_copper(end)*1e6);
fprintf('  R_AC (Aluminium): %.4f μΩ\n', R_AC_aluminum(end)*1e6);

%% Save the figure
output_filename = 'ac_resistance_comparison.png';
saveas(gcf, output_filename);
fprintf('\nPlot saved to: %s\n', fullfile(pwd, output_filename));
