%% DPT Waveform Analysis - Paralleled MOSFET Current Imbalance
% Imports CH1–CH4 from 3 DPT test folders, deduces per-MOSFET drain currents
% via rogowski coil subtraction, and plots I_D, I_G and V_GS during switch-off
% and switch-on events over a 1150 ns window.

%% Timeseries alignment strategy
% DPT pulsetrain is identical between empirical and simulation, but absolute timestamps are not. 
% empirical timestamps are fixed, since they correlate to the physically meaningful logic trigger 3.6 V rising edge
% LTspice timestamps will have an offset applied IMMEDIATELY after import, such that the timeseries are consistent for further anaylysis
% offset = ...
time_sim_offset = 9.9581e-6

%% DPT Timing & Circuit Parameters

T1  = 35;   % first pulse duration (µs) — switch-off event occurs at T1
T2  = 5;    % dead time duration (µs)
T3  = 5;    % second pulse duration (µs) — switch-on event occurs at T1+T2
R_G = 5.1;  % gate resistance (Ohms)

% Consistent MOSFET colours: blue / orange-red / green
colors = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]};


%% Plot Control Flags
% Set to true/false to enable/disable each figure group

plt_gate_current_full    = 0;   % "Full Sample Window - Gate Currents"
plt_gate_current_events  = 1;   % "Gate Current - Switch-Off/On"
plt_vgs_compare          = 0;   % "VGS Comparison - Switch-Off/On"
plt_raw_rogowski_events  = 0;   % "Raw Rogowski & Gate Currents - Switch-Off/On"
plt_corr_rogowski_events = 0;   % "Corrected Rogowski - Switch-Off/On"
plt_raw_rogowski_full    = 0;   % "Full Sample Window - Raw Rogowski Measurements"
plt_drain_current_full   = 0;   % "Full Sample Window - Drain Currents"
plt_drain_current_events = 1;   % "Switch-Off/On" drain current + VGS

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

%% Section 2: EMPIRICAL Import Rogowski & VGS/VDS Channels

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

%% Section 3: EMPIRICAL Import Gate Current Channels

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

%% Section 4: SIMULATION Import LTspice Channels

ltspice_csv = fullfile(base_path, 'LTspiceExport_35-5-5_48VDC_36uH.csv');
sim_raw = readtable(ltspice_csv);

time_sim   = sim_raw{:,1} - time_sim_offset;   % offset: shift LTspice timestamps by -10 us
VDS2_sim   = sim_raw{:,2};
VDS3_sim   = sim_raw{:,3};
PWLdriver_sim = sim_raw{:,4};
VGS1_sim   = sim_raw{:,5};
VGS2_sim   = sim_raw{:,6};
VGS3_sim   = sim_raw{:,7};
VDS1_sim   = sim_raw{:,8};
I_total_sim = sim_raw{:,9};
ID1_sim   = sim_raw{:,10};
ID2_sim   = sim_raw{:,11};
ID3_sim   = sim_raw{:,12};
IG1_sim   = sim_raw{:,13};
IG2_sim   = sim_raw{:,14};
IG3_sim   = sim_raw{:,15};

%% Section 5: EMPIRICAL Gate Current Plots

if plt_gate_current_full
    figure('Name', 'Full Sample Window - Gate Currents', 'NumberTitle', 'off');
    hold on;
    plot(time * 1e6, I_G1, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 1');
    plot(time * 1e6, I_G2, 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 2');
    plot(time * 1e6, I_G3, 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 3');
    hold off;
    xlabel('Time (\mus)');
    ylabel('Gate Current I_G (A)');
    title('Per-MOSFET Gate Currents - Full Sample Window');
    legend('Location', 'best');
    grid on;
end

if plt_gate_current_events
    plot_ig_event(time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, ...
                  time_sim, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                  T1,     'Switch-Off', colors, timewindow);
    plot_ig_event(time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, ...
                  time_sim, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                  T1+T2,  'Switch-On',  colors, timewindow);
end

if plt_vgs_compare
    plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                     time, VGS1_new, VGS2_new, VGS3_new, ...
                     T1,    'Switch-Off', colors, timewindow);
    plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                     time, VGS1_new, VGS2_new, VGS3_new, ...
                     T1+T2, 'Switch-On',  colors, timewindow);
end

%% Section 6: Rogowski Coil Correction [EDIT AS NEEDED]
% All signals share the same timebase - add or subtract I_G terms directly.
% Default: no correction applied.

I_total_corr = I_total;
I_23_corr    = I_23;
I_3_corr     = I_3;

if plt_raw_rogowski_events
    plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1,     'Switch-Off', colors, timewindow);
    plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1+T2,  'Switch-On',  colors, timewindow);
end

if plt_corr_rogowski_events
    plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1,     'Switch-Off', colors, timewindow);
    plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1+T2,  'Switch-On',  colors, timewindow);
end

%% Section 7: Per-MOSFET Drain Current Deduction

I_D3 = I_3_corr;
I_D2 = I_23_corr    - I_D3;
I_D1 = I_total_corr - I_23_corr;

%% Section 8: Drain Current Plots

if plt_raw_rogowski_full
    figure('Name', 'Full Sample Window - Raw Rogowski Measurements', 'NumberTitle', 'off');
    hold on;
    plot(time * 1e6, I_total, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'I_{total} (all 3 MOSFETs)');
    plot(time * 1e6, I_23,    'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'I_{D2+D3}');
    plot(time * 1e6, I_3,     'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'I_{D3}');
    hold off;
    xlabel('Time (\mus)');
    ylabel('Current (A)');
    title('Raw Rogowski Coil Measurements - Full Sample Window');
    legend('Location', 'best');
    grid on;
end

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

function plot_raw_and_ig_event(t, I_total, I_23, I_3, I_G1, I_G2, I_G3, t_event_us, label, colors, timewindow)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', ['Raw Rogowski & Gate Currents - ' label], 'NumberTitle', 'off');
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
    title(['Raw Rogowski & Gate Currents - ' label]);
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end

function plot_corr_event(t, I_total_corr, I_23_corr, I_3_corr, t_event_us, label, colors, timewindow)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + 1150e-9;
    mask    = t >= t_start & t <= t_end;
    t_ns    = (t(mask) - t_start) * 1e9;

    figure('Name', ['Corrected Rogowski - ' label], 'NumberTitle', 'off');
    hold on;
    plot(t_ns, I_total_corr(mask), 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'I_{total,corr}');
    plot(t_ns, I_23_corr(mask),    'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'I_{23,corr}');
    plot(t_ns, I_3_corr(mask),     'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'I_{3,corr}');
    hold off;
    xlabel('Time (ns)');
    ylabel('Current (A)');
    title(['Corrected Rogowski Signals - ' label]);
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end

function plot_vgs_compare(t_orig, VGS1_orig, VGS2_orig, VGS3_orig, ...
                          t_new,  VGS1_new,  VGS2_new,  VGS3_new, ...
                          t_event_us, label, colors, timewindow)
    t_start  = t_event_us * 1e-6;
    t_end    = t_start + 1150e-9;
    mk_orig  = t_orig >= t_start & t_orig <= t_end;
    mk_new   = t_new  >= t_start & t_new  <= t_end;
    t_ns_o   = (t_orig(mk_orig) - t_start) * 1e9;
    t_ns_n   = (t_new(mk_new)   - t_start) * 1e9;

    vgs_orig = {VGS1_orig, VGS2_orig, VGS3_orig};
    vgs_new  = {VGS1_new,  VGS2_new,  VGS3_new};
    labels   = {'MOSFET 1', 'MOSFET 2', 'MOSFET 3'};

    figure('Name', ['VGS Comparison - ' label], 'NumberTitle', 'off');
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['V_{GS} Comparison - ' label], 'FontWeight', 'bold');

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
        xlim([0, timewindow]);
        grid on;
        if k == 3, xlabel('Time (ns)'); end
    end
end

function plot_ig_event(t, I_G1, I_G2, I_G3, VGS1, VGS2, VGS3, ...
                       t_sim, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                       t_event_us, label, colors, timewindow)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + timewindow * 1e-9;

    mask     = t     >= t_start & t     <= t_end;
    mask_sim = t_sim >= t_start & t_sim <= t_end;
    t_ns     = (t(mask)         - t_start) * 1e9;
    t_ns_sim = (t_sim(mask_sim) - t_start) * 1e9;

    figure('Name', ['Gate V/I - ' label], 'NumberTitle', 'off');
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ['Gate V/I - ' label], 'FontWeight', 'bold');

    nexttile;
    hold on;
    plot(t_ns,     VGS1(mask),         '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
    plot(t_ns,     VGS2(mask),         '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
    plot(t_ns,     VGS3(mask),         '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
    plot(t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
    plot(t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
    plot(t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    hold off;
    ylabel('V_{GS} (V)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;

    nexttile;
    hold on;
    plot(t_ns,     I_G1(mask),         '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
    plot(t_ns,     I_G2(mask),         '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
    plot(t_ns,     I_G3(mask),         '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
    plot(t_ns_sim, IG1_sim(mask_sim),  '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
    plot(t_ns_sim, IG2_sim(mask_sim),  '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
    plot(t_ns_sim, IG3_sim(mask_sim),  '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    hold off;
    xlabel('Time (ns)');
    ylabel('I_G (A)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end

function plot_event(t, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, t_event_us, label, colors)
    t_start = t_event_us * 1e-6;
    t_end   = t_start + timewindow * 1e-9;

    mask     = t     >= t_start & t     <= t_end;
    mask_sim = t_sim >= t_start & t_sim <= t_end;
    t_ns     = (t(mask)         - t_start) * 1e9;
    t_ns_sim = (t_sim(mask_sim) - t_start) * 1e9;

    figure('Name', label, 'NumberTitle', 'off');
    tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, label, 'FontWeight', 'bold');

    nexttile;
    hold on;
    plot(t_ns,     VGS1(mask),         '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
    plot(t_ns,     VGS2(mask),         '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
    plot(t_ns,     VGS3(mask),         '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
    plot(t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
    plot(t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
    plot(t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    hold off;
    ylabel('V_{GS} (V)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;

    nexttile;
    hold on;
    plot(t_ns,     I_D1(mask),         '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
    plot(t_ns,     I_D2(mask),         '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
    plot(t_ns,     I_D3(mask),         '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
    plot(t_ns_sim, ID1_sim(mask_sim),  '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
    plot(t_ns_sim, ID2_sim(mask_sim),  '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
    plot(t_ns_sim, ID3_sim(mask_sim),  '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    hold off;
    xlabel('Time (ns)');
    ylabel('Drain Current I_D (A)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end

function plot_vds_id_report(time, VDS, I_D, T1, T2, off_pts, on_pts)
% Report-quality double-column figure: V_DS (left y) and I_{D,total} (right y)
% for turn-off and turn-on events side by side.

    t_off = T1        * 1e-6;
    t_on  = (T1 + T2) * 1e-6;
    win   = 1e-6;

    mk_off = time >= t_off & time <= (t_off + win);
    mk_on  = time >= t_on  & time <= (t_on  + win);
    ns_off = @(t) (t - t_off) * 1e9;
    ns_on  = @(t) (t - t_on)  * 1e9;
    xlims  = [0, win * 1e9];

    c_vds = [0.20 0.40 0.80];
    c_id  = [0.80 0.30 0.10];
    fs    = 8;    % font size (pt) — matches IEEE body text
    lw    = 1.2;  % line width

    fig = figure('Name', 'VDS & I_Dtotal — Report', 'NumberTitle', 'off');
    fig.Units    = 'centimeters';
    fig.Position = [2 2 18.5 7.5];   % double-column width × compact height

    tl = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % ---- Tile 1: Turn-Off ----
    ax1 = nexttile(1);
    yyaxis left;  hold on;
    plot(ns_off(time(mk_off)), VDS(mk_off), 'Color', c_vds, 'LineWidth', lw);
    yline(0.10 * off_pts.VDS_ss, '--', 'Color', c_vds, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    yline(0.90 * off_pts.VDS_ss, '--', 'Color', c_vds, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    if ~isnan(off_pts.t_VDS_10)
        scatter(ns_off(off_pts.t_VDS_10), sample_at(time, VDS, off_pts.t_VDS_10), 30, c_vds, 'filled', 'HandleVisibility', 'off');
    end
    if ~isnan(off_pts.t_VDS_90)
        scatter(ns_off(off_pts.t_VDS_90), sample_at(time, VDS, off_pts.t_VDS_90), 30, c_vds, 'filled', 'HandleVisibility', 'off');
    end
    ylabel('V_{DS} (V)', 'FontSize', fs);
    yyaxis right;  hold on;
    plot(ns_off(time(mk_off)), I_D(mk_off), 'Color', c_id, 'LineWidth', lw);
    yline(0.90 * off_pts.I_ss, '--', 'Color', c_id, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    yline(0.10 * off_pts.I_ss, '--', 'Color', c_id, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    if ~isnan(off_pts.t_I_90)
        scatter(ns_off(off_pts.t_I_90), sample_at(time, I_D, off_pts.t_I_90), 30, c_id, 'filled', 'HandleVisibility', 'off');
    end
    if ~isnan(off_pts.t_I_10)
        scatter(ns_off(off_pts.t_I_10), sample_at(time, I_D, off_pts.t_I_10), 30, c_id, 'filled', 'HandleVisibility', 'off');
    end
    ylabel('I_{D} (A)', 'FontSize', fs);
    ax1.YAxis(1).Color = c_vds;  ax1.YAxis(2).Color = c_id;
    xlim(xlims);  grid on;  ax1.FontSize = fs;
    xlabel('Time (ns)', 'FontSize', fs);
    title(sprintf('Turn-Off  —  t_d = %.0f ns,  t_f = %.0f ns', off_pts.t_d_ns, off_pts.t_sw_ns), ...
          'FontSize', fs, 'FontWeight', 'normal');

    % ---- Tile 2: Turn-On ----
    ax2 = nexttile(2);
    yyaxis left;  hold on;
    plot(ns_on(time(mk_on)), VDS(mk_on), 'Color', c_vds, 'LineWidth', lw);
    yline(0.90 * on_pts.VDS_ss, '--', 'Color', c_vds, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    yline(0.10 * on_pts.VDS_ss, '--', 'Color', c_vds, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    if ~isnan(on_pts.t_VDS_90)
        scatter(ns_on(on_pts.t_VDS_90), sample_at(time, VDS, on_pts.t_VDS_90), 30, c_vds, 'filled', 'HandleVisibility', 'off');
    end
    if ~isnan(on_pts.t_VDS_10)
        scatter(ns_on(on_pts.t_VDS_10), sample_at(time, VDS, on_pts.t_VDS_10), 30, c_vds, 'filled', 'HandleVisibility', 'off');
    end
    ylabel('V_{DS} (V)', 'FontSize', fs);
    yyaxis right;  hold on;
    plot(ns_on(time(mk_on)), I_D(mk_on), 'Color', c_id, 'LineWidth', lw);
    yline(0.10 * on_pts.I_ss, '--', 'Color', c_id, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    yline(0.90 * on_pts.I_ss, '--', 'Color', c_id, 'LineWidth', 0.7, 'HandleVisibility', 'off');
    if ~isnan(on_pts.t_I_10)
        scatter(ns_on(on_pts.t_I_10), sample_at(time, I_D, on_pts.t_I_10), 30, c_id, 'filled', 'HandleVisibility', 'off');
    end
    if ~isnan(on_pts.t_I_90)
        scatter(ns_on(on_pts.t_I_90), sample_at(time, I_D, on_pts.t_I_90), 30, c_id, 'filled', 'HandleVisibility', 'off');
    end
    ylabel('I_{D} (A)', 'FontSize', fs);
    ax2.YAxis(1).Color = c_vds;  ax2.YAxis(2).Color = c_id;
    xlim(xlims);  grid on;  ax2.FontSize = fs;
    xlabel('Time (ns)', 'FontSize', fs);
    title(sprintf('Turn-On  —  t_d = %.0f ns,  t_r = %.0f ns', on_pts.t_d_ns, on_pts.t_sw_ns), ...
          'FontSize', fs, 'FontWeight', 'normal');

    % ---- Equalise y-scales across both tiles then zero-align ----
    % Read auto-scaled limits from each tile
    yyaxis(ax1, 'left');   yl1L = ylim(ax1);
    yyaxis(ax1, 'right');  yl1R = ylim(ax1);
    yyaxis(ax2, 'left');   yl2L = ylim(ax2);
    yyaxis(ax2, 'right');  yl2R = ylim(ax2);

    % Union: widest range that covers both tiles
    ylL = [min(yl1L(1), yl2L(1)), max(yl1L(2), yl2L(2))];
    ylR = [min(yl1R(1), yl2R(1)), max(yl1R(2), yl2R(2))];

    % Apply union to ax1, zero-align, then read back the final limits
    yyaxis(ax1, 'left');  ylim(ax1, ylL);
    yyaxis(ax1, 'right'); ylim(ax1, ylR);
    align_yyaxis_zeros(ax1);
    yyaxis(ax1, 'left');  ylL_f = ylim(ax1);
    yyaxis(ax1, 'right'); ylR_f = ylim(ax1);

    % Copy identical limits to ax2 (zero alignment is automatic since limits match ax1)
    yyaxis(ax2, 'left');  ylim(ax2, ylL_f);
    yyaxis(ax2, 'right'); ylim(ax2, ylR_f);
end
