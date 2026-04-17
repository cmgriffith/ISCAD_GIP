%% DPT Waveform Analysis — Paralleled MOSFET Current Imbalance
% Imports CH1–CH4 from 3 DPT test folders, deduces per-MOSFET drain currents
% via rogowski coil subtraction, and plots I_D and V_GS during switch-off
% (35 µs) and switch-on (40 µs) events over a 1150 ns window.

%% Section 1: Paths and Folder Definitions

base_path = 'C:\Users\ossia\OneDrive - University of Bristol\grp-GRP Group 1002 - Documents\DPT Results\Tests at Peak Current';

folders = { ...
    '35-5-5DPT_48VDC_logic_VGS1_IDtotal_VDS_36uHinductor', ...
    '35-5-5DPT_48VDC_logic_VGS2_ID2and3_VDS_36uHinductor', ...
    '35-5-5DPT_48VDC_logic_VGS3_ID3_VDS3_36uHinductor'};

load_ch = @(folder_idx, ch_num) load_csv( ...
    fullfile(base_path, folders{folder_idx}, sprintf('CH%d.CSV', ch_num)));

%% Section 2: Import All Channels

[time, PWM]     = load_ch(1, 1);   % shared time axis (trigger-aligned) + PWM

[~,    VGS1]    = load_ch(1, 2);   % V_GS MOSFET 1
[~,    VGS2]    = load_ch(2, 2);   % V_GS MOSFET 2
[~,    VGS3]    = load_ch(3, 2);   % V_GS MOSFET 3

[~,    I_total] = load_ch(1, 3);   % I_D1 + I_D2 + I_D3 (rogowski at all 3)
[~,    I_23]    = load_ch(2, 3);   % I_D2 + I_D3 (rogowski at MOSFETs 2 & 3)
[~,    I_3]     = load_ch(3, 3);   % I_D3 (rogowski at MOSFET 3 only)

[~,    VDS1]    = load_ch(1, 4);   % V_DS MOSFET 1 (reserved for future analysis)
[~,    VDS2]    = load_ch(2, 4);   % V_DS MOSFET 2
[~,    VDS3]    = load_ch(3, 4);   % V_DS MOSFET 3

%% Section 3: Per-MOSFET Drain Current Deduction

I_D3 = I_3;
I_D2 = I_23    - I_D3;
I_D1 = I_total - I_23;

%% Section 4: Sanity Check — deduction residual (comment out after verification)
figure('Name', 'Deduction Residual');
plot(time * 1e6, I_D1 + I_D2 + I_D3 - I_total);
xlabel('Time (\mus)'); ylabel('Residual (A)');
title('I_{D1}+I_{D2}+I_{D3} - I_{total} (should be ~0)');

%% Section 5: Switching-Event Plots

% Consistent MOSFET colours: blue / orange-red / green
colors = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]};

% Full sample window overview
figure('Name', 'Full Sample Window — Drain Currents', 'NumberTitle', 'off');
hold on;
plot(time * 1e6, I_D1, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 1');
plot(time * 1e6, I_D2, 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 2');
plot(time * 1e6, I_D3, 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 3');
hold off;
xlabel('Time (\mus)');
ylabel('Drain Current I_D (A)');
title('Per-MOSFET Drain Currents — Full Sample Window');
legend('Location', 'best');
grid on;

plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, 35, 'Switch-Off', colors);
plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, 40, 'Switch-On',  colors);

%% -------------------------------------------------------------------------
%  Local Functions
%% -------------------------------------------------------------------------

function [t, v] = load_csv(filepath)
    data = readmatrix(filepath, 'NumHeaderLines', 1, 'Delimiter', ',');
    t = data(:, 1);
    v = data(:, 2);
end

function plot_event(t, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', label, 'NumberTitle', 'off');
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, label, 'FontWeight', 'bold');

    % --- Top tile: per-MOSFET drain currents ---
    nexttile;
    hold on;
    plot(t_ns, I_D1(mask), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 1');
    plot(t_ns, I_D2(mask), 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 2');
    plot(t_ns, I_D3(mask), 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 3');
    hold off;
    ylabel('Drain Current I_D (A)');
    legend('Location', 'best');
    xlim([0, 1150]);
    grid on;

    % --- Bottom tile: gate-source voltages ---
    nexttile;
    hold on;
    plot(t_ns, VGS1(mask), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 1');
    plot(t_ns, VGS2(mask), 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 2');
    plot(t_ns, VGS3(mask), 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 3');
    hold off;
    xlabel('Time (ns)');
    ylabel('V_{GS} (V)');
    legend('Location', 'best');
    xlim([0, 1150]);
    grid on;
end
