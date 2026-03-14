clear; clc;

% Load parameters from JSON file into workspace struct 'mparams'
mparams = load_json_parameters('../ISCAD_parameters.json');

% Constants
f_e = 50; % Electrical Fundamental Speed
f_sw = 25e3; % Switch Frequency
Td = 2e-6; % Deadtime
T_s = 1./(20*f_sw); % Sampling Frequency
Vdc = 48; % DC link voltage (V)
I_ph = 140; % RMS phase current (A)
Cdc = 500e-6; % DC link capacitance (F)

L_ph = 0.25e-3; % Phase Inductance
R_ph = 0.8e-3; % Phase Resistance @ max speed
Z = R_ph + 1i*2*pi*f_e*L_ph;

Z_mag = abs(Z);           % |Z| = sqrt(R^2 + X_L^2)
phi = angle(Z);           % angle = atan(X_L/R) [radians]

V_ph_rms = I_ph * Z_mag;
m = (V_ph_rms * 2 * sqrt(2)) / Vdc;
