clear; clc;

% Load parameters from JSON file into workspace struct 'mparams'
mparams = load_json_parameters('../TEST_parameters.json');
mu_0 = 4*pi*1e-7; % Permeability of free space [H/m]

%% Inverter / Electrical Parameters
f_e = 500;        % Electrical fundamental frequency [Hz]
k_AC = 3;        % AC Resistance Factor (rough approximation)
f_sw = 25e3;     % Switching frequency [Hz]
Td = 2e-6;       % Dead time [s]
T_s = 1/(20*f_sw); % Sampling period [s]
Vdc = 48;        % DC link voltage [V]
I_ph = 10;      % RMS phase current [A]
Cdc = 500e-6;    % DC link capacitance [F]

%% Resistances
Rs_DC = mparams.rho_s * mparams.l_bar ./ (mparams.h_bar * mparams.w_bar(1));
R_ph_AC = Rs_DC * k_AC;

% Inductance
L_ph = 0.25e-3;  % Phase inductance [H]

Z = R_ph_AC + 1i*2*pi*f_e*L_ph;

Z_mag = abs(Z);
phi = angle(Z);

V_ph_rms = I_ph * Z_mag;
mod_index = (V_ph_rms * 2 * sqrt(2)) / Vdc;

fprintf('Phase Resistance (AC): %.4f Ohms\n', R_ph_AC);
fprintf('Phase Impedance Magnitude: %.4f Ohms\n', Z_mag);
fprintf('Phase Impedance Angle: %.2f degrees\n', rad2deg(phi));
fprintf('Required Phase Voltage (RMS): %.2f V\n', V_ph_rms);
fprintf('Modulation Index: %.4f\n', mod_index);


%% LC Resonator Analysis
