%% DPT Waveform Analysis — Paralleled MOSFET Current Imbalance
% Imports CH1–CH4 from 3 DPT test folders, deduces per-MOSFET drain currents
% via rogowski coil subtraction, and plots I_D, I_G and V_GS during switch-off
% and switch-on events over a 1150 ns window.

%% Timing & Circuit Parameters

T1  = 35;   % first pulse duration (µs) — switch-off event occurs at T1
T2  = 5;    % dead time duration (µs)
T3  = 5;    % second pulse duration (µs) — switch-on event occurs at T1+T2
R_G = 5.1;  % gate resistance (Ohms)

% Consistent MOSFET colours: blue / orange-red / green
colors = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]};

%% Section 1: Paths and Folder Definitions

base_path = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Bristol', 'grp-GRP Group 1002 - Documents', 'DPT Results', 'Tests at Peak Current');

folders = { ...
    '35-5-5DPT_48VDC_logic_VGS1_IDtotal_VDS_36uHinductor', ...
    '35-5-5DPT_48VDC_logic_VGS2_ID2and3_VDS_36uHinductor', ...
    '35-5-5DPT_48VDC_logic_VGS3_ID3_VDS3_36uHinductor'};

asym_base = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Bristol', 'grp-GRP Group 1002 - Documents', 'DPT Results', 'Tests at Peak Current', 'Investigating Asymmetry');

asym_folders = { ...
    '35-5-5DPT_48VDC_logic_VGS1_VGRL1_VGRH1', ...
    '35-5-5DPT_48VDC_logic_VGS2_VGRL2_VGRH2', ...
    '35-5-5DPT_48VDC_logic_VGS3_VGRL3_VGRH3'};

load_ch   = @(f, ch) load_csv(fullfile(base_path,  folders{f},      sprintf('CH%d.CSV', ch)));
load_asym = @(f, ch) load_csv(fullfile(asym_base,  asym_folders{f}, sprintf('CH%d.CSV', ch)));

%% Section 2: Import Rogowski & VGS/VDS Channels

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

%% Section 3: Import Gate Current Channels

[~,    VGS1_new] = load_asym(1, 2);   % V_GS MOSFET 1 (re-measured)
[~,    V_b1]     = load_asym(1, 3);   % voltage before gate resistor, MOSFET 1
[~,    V_a1]     = load_asym(1, 4);   % voltage after gate resistor, MOSFET 1

[~,    VGS2_new] = load_asym(2, 2);
[~,    V_b2]     = load_asym(2, 3);
[~,    V_a2]     = load_asym(2, 4);

[~,    VGS3_new] = load_asym(3, 2);
[~,    V_b3]     = load_asym(3, 3);
[~,    V_a3]     = load_asym(3, 4);

% Voltage drop across each gate resistor: CH3 - CH4
V_RG1 = V_b1 - V_a1;
V_RG2 = V_b2 - V_a2;
V_RG3 = V_b3 - V_a3;

% Gate current via Ohm's law: I_G = V_RG / R_G
I_G1 = V_RG1 / R_G;
I_G2 = V_RG2 / R_G;
I_G3 = V_RG3 / R_G;

%% Section 4: Gate Current Plots

% Full sample window
figure('Name', 'Full Sample Window — Gate Currents', 'NumberTitle', 'off');
hold on;
plot(time * 1e6, I_G1, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 1');
plot(time * 1e6, I_G2, 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 2');
plot(time * 1e6, I_G3, 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 3');
hold off;
xlabel('Time (\mus)');
ylabel('Gate Current I_G (A)');
title('Per-MOSFET Gate Currents — Full Sample Window');
legend('Location', 'best');
grid on;

% Switching event windows
plot_ig_event(time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, T1,     'Switch-Off', colors);
plot_ig_event(time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, T1+T2, 'Switch-On',  colors);

% VGS comparison: original vs re-measured
plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                 time, VGS1_new, VGS2_new, VGS3_new, ...
                 T1,    'Switch-Off', colors);
plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                 time, VGS1_new, VGS2_new, VGS3_new, ...
                 T1+T2, 'Switch-On',  colors);

%% Section 5: Rogowski Coil Correction [EDIT AS NEEDED]
% All signals share the same timebase — add or subtract I_G terms directly.
% Default: no correction applied.

I_total_corr = I_total;
I_23_corr    = I_23;
I_3_corr     = I_3;

% Uncorrected rogowski signals with gate currents overlaid — switching event windows
plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1,     'Switch-Off', colors);
plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1+T2,  'Switch-On',  colors);

% Corrected rogowski signals — switching event windows
plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1,     'Switch-Off', colors);
plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1+T2,  'Switch-On',  colors);

%% Section 6: Per-MOSFET Drain Current Deduction

I_D3 = I_3_corr;
I_D2 = I_23_corr    - I_D3;
I_D1 = I_total_corr - I_23_corr;

%% Section 8: Drain Current Plots

% Raw rogowski coil measurements — full sample window
figure('Name', 'Full Sample Window — Raw Rogowski Measurements', 'NumberTitle', 'off');
hold on;
plot(time * 1e6, I_total, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'I_{total} (all 3 MOSFETs)');
plot(time * 1e6, I_23,    'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'I_{D2+D3}');
plot(time * 1e6, I_3,     'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'I_{D3}');
hold off;
xlabel('Time (\mus)');
ylabel('Current (A)');
title('Raw Rogowski Coil Measurements — Full Sample Window');
legend('Location', 'best');
grid on;

% Deduced per-MOSFET drain currents — full sample window
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

% Switching event windows
plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, T1,     'Switch-Off', colors);
plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, T1+T2,  'Switch-On',  colors);

%% -------------------------------------------------------------------------
%  Local Functions
%% -------------------------------------------------------------------------

function [t, v] = load_csv(filepath)
    data = readmatrix(filepath, 'NumHeaderLines', 1, 'Delimiter', ',');
    t = data(:, 1);
    v = data(:, 2);
end

function plot_raw_and_ig_event(t, I_total, I_23, I_3, I_G1, I_G2, I_G3, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', ['Raw Rogowski & Gate Currents — ' label], 'NumberTitle', 'off');
    hold on;
    plot(t_ns, I_total(mask), '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'I_{total}');
    plot(t_ns, I_23(mask),    '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'I_{23}');
    plot(t_ns, I_3(mask),     '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'I_{3}');
    plot(t_ns, I_G1(mask),    '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'I_{G1}');
    plot(t_ns, I_G2(mask),    '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'I_{G2}');
    plot(t_ns, I_G3(mask),    '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'I_{G3}');
    hold off;
    xlabel('Time (ns)');
    ylabel('Current (A)');
    title(['Raw Rogowski & Gate Currents — ' label]);
    legend('Location', 'best');
    xlim([0, 1150]);
    grid on;
end

function plot_corr_event(t, I_total_corr, I_23_corr, I_3_corr, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', ['Corrected Rogowski — ' label], 'NumberTitle', 'off');
    hold on;
    plot(t_ns, I_total_corr(mask), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'I_{total,corr}');
    plot(t_ns, I_23_corr(mask),    'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'I_{23,corr}');
    plot(t_ns, I_3_corr(mask),     'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'I_{3,corr}');
    hold off;
    xlabel('Time (ns)');
    ylabel('Current (A)');
    title(['Corrected Rogowski Signals — ' label]);
    legend('Location', 'best');
    xlim([0, 1150]);
    grid on;
end

function plot_vgs_compare(t_orig, VGS1_orig, VGS2_orig, VGS3_orig, ...
                          t_new,  VGS1_new,  VGS2_new,  VGS3_new, ...
                          t_event_us, label, colors)
    t_start  = t_event_us * 1e-6;
    t_end    = t_start + 1150e-9;
    mk_orig  = t_orig >= t_start & t_orig <= t_end;
    mk_new   = t_new  >= t_start & t_new  <= t_end;
    t_ns_o   = (t_orig(mk_orig) - t_start) * 1e9;
    t_ns_n   = (t_new(mk_new)   - t_start) * 1e9;

    vgs_orig = {VGS1_orig, VGS2_orig, VGS3_orig};
    vgs_new  = {VGS1_new,  VGS2_new,  VGS3_new};
    labels   = {'MOSFET 1', 'MOSFET 2', 'MOSFET 3'};

    figure('Name', ['VGS Comparison — ' label], 'NumberTitle', 'off');
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['V_{GS} Comparison — ' label], 'FontWeight', 'bold');

    for k = 1:3
        nexttile;
        hold on;
        plot(t_ns_o, vgs_orig{k}(mk_orig), '--', 'Color', colors{k}, ...
             'LineWidth', 1.2, 'DisplayName', 'Original');
        plot(t_ns_n, vgs_new{k}(mk_new),   '-',  'Color', colors{k}, ...
             'LineWidth', 1.5, 'DisplayName', 'Re-measured');
        hold off;
        title(labels{k});
        ylabel('V_{GS} (V)');
        legend('Location', 'best');
        xlim([0, 1150]);
        grid on;
        if k == 3, xlabel('Time (ns)'); end
    end
end

function plot_ig_event(t, I_G1, I_G2, I_G3, VGS1, VGS2, VGS3, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', ['Gate Current — ' label], 'NumberTitle', 'off');
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['Gate Current — ' label], 'FontWeight', 'bold');

    nexttile;
    hold on;
    plot(t_ns, I_G1(mask), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 1');
    plot(t_ns, I_G2(mask), 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 2');
    plot(t_ns, I_G3(mask), 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'MOSFET 3');
    hold off;
    ylabel('Gate Current I_G (A)');
    legend('Location', 'best');
    xlim([0, 1150]);
    grid on;

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

function plot_event(t, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', label, 'NumberTitle', 'off');
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, label, 'FontWeight', 'bold');

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
