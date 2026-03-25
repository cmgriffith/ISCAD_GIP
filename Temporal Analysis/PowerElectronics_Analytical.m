clear; clc;

% Load parameters from JSON file into workspace struct 'mparams'
mparams = load_json_parameters('../TEST_parameters.json');
mu_0 = 4*pi*1e-7; % Permeability of free space [H/m]

%% controlled variables
Vdc = 48;   % DC link voltage [V]
Iph = 10;   % phase current (RMS?)
f_e;    % fundamental frequency
f_sw;

%% operating point dependent variables
R_AC_ph; % ANALYTICAL FROM DRM..
L_ph; % FEMM NUMERICAL lookup table L_ph(Iph, f)

%% PCB design dependent variables
% passive parameters
C_snub
C_DClink = 500e-6;    % DC link capacitance [F]

% IPTG011N08NM5 device parameters
R_DSon = 1e-3;
R_G = 2;
% at V_DS = 48V
C_iss = 13000e-12; % = C_GD + C_GS
C_oss = 1800e-12; % = C_DS + C_GD
C_rss = 70e-12; % = C_GD
C_GD = C_rss;
C_GS = C_iss - C_GD;
C_DS = C_oss - C_GD;
t_don = 35e-9; % rise delay
t_r = 31e-9; % rise
t_doff = 82e-9; % fall delay
t_f = 30e-9; % fall
Td = 2e-6;  % Dead time [s]
T_s = 1/(20*f_sw); % Sampling period [s]

% PCB parasitics
R_pcb % bulk parasitic resistance ANAYLTICAL FROM PCB DESIGN
L_pcb % bulk parasitic inductance ANAYLTICAL FROM PCB DESIGN
R_snub % parasitic snubber resistance
L_snub % parasitic snubber inductance
R_DClink % DC link connection parasitic
L_DClink % DC link connection parasitic
R_gatedrive = % gate trace resistance
L_gatedrive = % gate trace inductance





%% Load
% Resistances
Rs_DC = mparams.rho_s * mparams.l_bar ./ (mparams.h_bar * mparams.w_bar(1));
R_ph_AC = Rs_DC * k_AC;

% Inductances
% INTERPOLATE LOOKUP TABLE FROM FEMM 



%% Electronics
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
