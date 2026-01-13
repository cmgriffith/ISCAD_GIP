function [Ls_s, Ls_sigma_DC, Ls, Ls_m_arr] = calculate_inductance()
% CALCULATE_INDUCTANCE Calculate stator inductances from parameters
%
% This MATLAB implementation matches the calculations in ISCAD_inductance.py
%
% Outputs:
%   Ls_s         - Self inductance [H]
%   Ls_sigma_DC  - Leakage inductance DC [H]
%   Ls_m         - Main inductance [H]
%   Ls_m_arr     - Main inductance array for different slot counts

    % Set Parameters
    Qs = 60; % Number of stator slots
    p = 2; % Number of polls
    r_bore = 75e-3; % Inner stator radius
    l_stack = 250e-3; % Stack length
    ag = 0.5e-3; % What is this again?
    b_slot = 4e-3; % Slot Width
    b_so = 16e-3; % Slot opening width
    h_slot = 1e-3; % slot height
    h_so = 1e-3; % Slot opening height
    mu_0 = 4*pi*1e-7;  % permeability of free space [H/m]
    
    %% SELF INDUCTANCE
    Ls_s = (pi/2) * mu_0 * (r_bore * l_stack / ag);
    fprintf('Self inductance = %.3f H\n', Ls_s);
    
    %% LEAKAGE INDUCTANCE
    lambda_slot = (h_slot / (3*b_slot)) + (h_so / b_so);
    Ls_sigma_DC = mu_0 * l_stack * lambda_slot;
    fprintf('Leakage inductance DC = %.3f H\n', Ls_sigma_DC);
    
    %% MAIN INDUCTANCE
    Ls_arr = zeros(Qs, 1);
    
    % Main inductances for range of Qs
    for slot = 1:Qs
        Qs_loop = slot;
        m_loop = Qs_loop / p;
        sum_val = 0;
        
        for k = 1:round(m_loop/2)
            x = cos(p * (k-1) * (2*pi / Qs_loop))^2;
            sum_val = sum_val + x;
        end
        
        Ls = Ls_s * sum_val;
        Ls_arr(slot) = Ls;
    end
    
    Ls = Ls_arr(end);  % Final value for full Qs
    fprintf('Main inductance = %.3f H\n', Ls);
    
    %% PLOT
    figure;
    bar(0:Qs-1, Ls_arr);
    xlabel('Q_s (slot count)');
    ylabel('L_m [H] main inductance');
    title('Main inductance vs Q_s');
    hold on;
    yline(Ls_s, '--r', 'LineWidth', 1.5, 'DisplayName', 'L_s self inductance');
    legend('show');
    hold off;
    
end
