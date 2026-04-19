%% DC Link Capacitor RMS Current - N-phase Inverter
% Based on Schulte et al. (2004), PWM mode only
% Reproduces normalised RMS capacitor current vs modulation index
% for varying phase numbers (cf. Patzak Fig. 3)
 
clear; close all;
 
%% Sweep parameters
m_vec   = linspace(0.1, 1, 100);   % modulation index range
phi     = 0;                         % power factor angle (0 = unity PF)
I_peak  = 1;                         % normalised peak phase current
N_vals  = [3, 15, 30, 60];           % phase numbers to sweep
 
colors  = {'#185FA5','#D85A30','#1D9E75','#7F77DD'};
 
figure; hold on; grid on;
 
for n = 1:length(N_vals)
    N = N_vals(n);
    I_C  = zeros(1, length(m_vec));
    I_DC = zeros(1, length(m_vec));
 
    for i = 1:length(m_vec)
        [I_C(i), I_DC(i)] = dc_link_cap_rms(N, m_vec(i), phi, I_peak);
    end
 
    % Normalise by DC link average current (matches Patzak Fig. 3 y-axis)
    plot(m_vec, I_C ./ I_DC, 'Color', colors{n}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('N = %d', N));
end
 
xlabel('Modulation index m');
ylabel('I_{C,rms} / I_{DC}');
title('Normalised DC link capacitor RMS current (unity power factor)');
legend('Location','northeast');
ylim([0 2]);
 
 
%% Function: dc_link_cap_rms
% Implements Schulte et al. PWM mode analytical method
% Returns I_C_rms and I_link_avg for one operating point
 
function [I_C_rms, I_link_avg] = dc_link_cap_rms(N, m, phi, I_peak)
 
    % --- Integration interval centre (Schulte eq. 3) ---
    % Odd N: centre at pi
    % Even N: centre at pi + pi/(2N)
    if mod(N, 2) == 0
        t_centre = pi + pi/(2*N);
    else
        t_centre = pi;
    end
 
    n_pts = 10000;
    t  = linspace(t_centre - pi/(2*N), t_centre + pi/(2*N), n_pts);
    dt = t(2) - t(1);
 
    k = (0:N-1)';   % phase indices, column vector for broadcasting
 
    % --- Relative conducting times (Schulte eq. 2) ---
    % a is N x n_pts matrix: row k = phase k, column j = time point j
    a = 0.5 * (1 + m * sin(t - 2*k*pi/N));
 
    % --- Phase currents (N x n_pts) ---
    i_phase = I_peak * sin(t - 2*k*pi/N - phi);
 
    % --- Sort RTCs in descending order at each time point ---
    % a_sorted(k,j) = k-th largest RTC at time j
    % sort_idx(k,j) = which phase has the k-th largest RTC at time j
    [a_sorted, sort_idx] = sort(a, 1, 'descend');
 
    % --- Reorder phase currents to match sorted RTCs ---
    % Vectorised linear indexing: converts (row,col) pairs to flat indices
    col_offset = N * repmat(0:n_pts-1, N, 1);
    i_sorted   = i_phase(sort_idx + col_offset);
 
    % --- Differences between adjacent sorted RTCs ---
    % delta_a(k,j) = fraction of switching period where exactly k switches
    % are simultaneously conducting (k = 1..N)
    % Append a row of zeros so delta_a(N) = a_sorted(N) - 0
    a_aug   = [a_sorted; zeros(1, n_pts)];
    delta_a = a_aug(1:N,:) - a_aug(2:N+1,:);
 
    % --- Cumulative phase current sums ---
    % i_cumsum(k,j) = sum of top-k phase currents at time j
    % This is the DC link current when exactly k upper switches conduct
    i_cumsum = cumsum(i_sorted, 1);
 
    % --- DC link current integrands (Schulte eqs. 4, 5) ---
    % At each time j, the switching period is divided into sub-intervals
    % of duration delta_a(k,j), each carrying current i_cumsum(k,j)
    i_avg_integrand  = sum(delta_a .* i_cumsum,    1);   % for avg
    i_rms2_integrand = sum(delta_a .* i_cumsum.^2, 1);   % for rms²
 
    % --- Integrate over fundamental interval pi/N (Schulte eqs. 6, 7) ---
    % Normalisation factor = 1/(pi/N) = N/pi
    I_link_avg  =        (N/pi) * sum(i_avg_integrand)  * dt;
    I_link_rms2 =        (N/pi) * sum(i_rms2_integrand) * dt;
 
    % --- Capacitor RMS current (Schulte eq. 1) ---
    % I_C² = I_rms² - I_avg²
    % max(0,...) guards against small negative floating point errors
    I_C_rms = sqrt(max(0, I_link_rms2 - I_link_avg^2));
 
end