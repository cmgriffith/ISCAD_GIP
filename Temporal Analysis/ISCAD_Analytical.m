% ISCAD_Analytical - Analytical calculations for ISCAD motor
% Loads parameters from JSON file into 'mparams' struct and performs analytical calculations

clear; clc;

% Load parameters from JSON file into workspace struct 'mparams'
mparams = load_json_parameters('../ISCAD_parameters.json');

% Constants
mu_0 = 4*pi*1e-7; % Permeability of free space [H/m]
f_r = 0; % Rotor Speed
mparams.p = 2;
Vs_rms = 0.03; % Desired stator voltage input (RMS)

%% RESISTANCES
% Use mparams.variable_name syntax
h_slot = 0.01:0.001:0.025; 
Rs_DC = mparams.rho_s * mparams.l_stack ./ (h_slot * mparams.w_slot(1));
f = linspace(0, 500, 50);
figure
hold on
for j = 1:length(h_slot)
    if f < 1
        K_Rs = 1;
    else
        zeta_Rs = sqrt((pi * mu_0 * f) ./ mparams.rho_s) .* h_slot(j);
        K_Rs = zeta_Rs .* (sinh(2*zeta_Rs) + sin(2*zeta_Rs)) ./ (cosh(2*zeta_Rs) - cos(2*zeta_Rs));
    end
    plot(f, K_Rs)
end
hold off
fprintf('Stator DC resistance = %.2e Ohm\n', Rs_DC);
Rs_AC = Rs_DC .* K_Rs;
figure
plot(f, K_Rs)

%%

% As rotor bars are not defined yet, we will strat by setting the total
% cross sectional area of the bars to be equal to the stator bars to have
% equal reistance.
A_rotorbar = mparams.h_slot * mparams.w_slot(1) * mparams.Qs / mparams.Qr;
Rr_DC = mparams.rho_s * mparams.l_stack / A_rotorbar;
zeta_Rr = sqrt((pi*mu_0*mparams.f)/mparams.rho_s) * (sqrt(4*A_rotorbar/(pi))); % Assumes circular cross section of rotor bars
K_Rr = zeta_Rr .* (sinh(2*zeta_Rr) + sin(2*zeta_Rr)) ./ (cosh(2*zeta_Rr) - cos(2*zeta_Rr));

%% SELF INDUCTANCE
Ls_s = pi/2 * mu_0 * (mparams.r * mparams.l_stack / mparams.ag);
fprintf('Self inductance = %.3f H\n', Ls_s);

%% LEAKAGE INDUCTANCE
% For DC current in rectangular slots
f_sw = linspace(100, 40000, 200);
lambda_slot = (mparams.h_slot / (3*mparams.w_slot(1))) + (mparams.h_so / mparams.w_so);
Ls_sigma_DC = mu_0 * mparams.l_stack * lambda_slot;
zeta_Ls = sqrt(pi*mu_0*f_sw/mparams.rho_s)*mparams.h_slot;
K_Ls = (3./(2*zeta_Ls)).*(sinh(2*zeta_Ls)-sin(2*zeta_Ls))./(cosh(2*zeta_Ls)-cos(2*zeta_Ls));
Ls_sigma_AC = mu_0 * mparams.l_stack * lambda_slot * K_Ls;

fprintf('Leakage inductance DC = %.3f H\n', Ls_sigma_DC);

figure
plot(f_sw, Ls_sigma_AC)


%% MAIN INDUCTANCE
sum_val = 0;
for k = 1:round(mparams.Qs/2/mparams.p)
    sum_val = sum_val + cos(mparams.p * (k-1) * 2*pi / mparams.Qs)^2;
end
Ls_m = Ls_s * sum_val;

fprintf('Main inductance = %.3f H\n', Ls_m);
