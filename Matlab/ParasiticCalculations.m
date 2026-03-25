l_bus = 90e-3; % DC Bus Bar Length on Half Bridge
t_bus = 0.6e-3; % DC Bus Bar Thickness
w_bus = 12e-3; % DC Bus Bar Width
d_seperation = 12e-3; % DC Bus Bar Seperation
d_datum = 30e-3; %w_bus+d_seperation; % To match with definition of d in Rubley paper
mu_0 = 4*pi*1e-7; % Permeability of free space [H/m]

L_p = (mu_0/(2*pi))*l_bus*(log(2*l_bus/(w_bus+t_bus))+0.5+(2/9)*((w_bus+t_bus)/l_bus));
M_p = (mu_0/(2*pi))*l_bus*(log((1/d_datum)+sqrt(1+(l_bus^2/d_datum^2))) ...
      - sqrt(1+(d_datum^2/l_bus^2)) + (d_datum/l_bus));
L_CON = L_p-M_p;

fprintf('Bus Bar Inductance (L_LCON,L_HCON): %.2f nH\n', L_CON*10^9); % Convert to nano-Henries

%% Parasitic Inductance of Power Loop
% From B. Sun, K. L. Jorgensen, Z. Zhang, and M. A. E. Andersen, 
% "Research of Power Loop Layout and Parasitic Inductance in GaN Transistor Implementation," 
% IEEE Transactions on Industry Applications, vol. 57, no. 2, pp. 1677–1687, 
% Mar. 2021, doi: 10.1109/TIA.2020.3048641.

% With new Vertical Layout
y_powerloop = d_seperation; % Assuming power loop exists in vertical layout and is via'd right on bus edges
t_PCB = 1.5e-3; % PCB thickness
w_PCB = 8.7e-3; % Width of MOSFET trace (need to clarify if this is correct!)

L_PL = mu_0*t_PCB*y_powerloop*((0.27)/(1-0.74*exp(-0.45*(t_PCB/w_PCB))))/(w_PCB);
fprintf('Power Loop Inductance (L_PL): %.2f nH\n', L_PL*10^9); % Convert to nano-Henries

Coss = 2000e-12;
rho_Cu = 1.724e-8; % Copper resistivity @ 20C
f_ring = 1/(2*pi*sqrt(L_PL*Coss));
delta_ringing = sqrt(rho_Cu/(pi*f_ring*mu_0));
fprintf('Power Loop Ringing Frequency (f_ring): %.2f MHz\n', f_ring/10^6);
fprintf('Copper Skin Depth at Ringing Frequency (delta_ringing): %.2f microns\n', delta_ringing*10^6); % Convert to nano-Henries





