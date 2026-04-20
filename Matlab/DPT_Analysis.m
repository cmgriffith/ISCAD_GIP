%% DPT Waveform Analysis — Paralleled MOSFET Current Imbalance
% Imports CH1–CH4 from 3 DPT test folders, deduces per-MOSFET drain currents
% via rogowski coil subtraction, and plots I_D, I_G and V_GS during switch-off
% and switch-on events over a 1150 ns window.

%% Timeseries alignment strategy
% DPT pulsetrain periods are identical between empirical and simulation, but absolute timestamps are not. 
% EMPIRICAL timestamps stay fixed, since they correlate to the physically meaningful logic trigger on a 3.6V T1 rising edge
% SIMULATION timestamps are different, and simulation doesn't model logic signal or the gate driver chip,
% hence simulation data will be aligned with 1% changes in the average device gate voltage, such that gate voltage waveforms are consistently timed

close all

%% Timing & Circuit Parameters

T1  = 35;   % first pulse duration (µs) — switch-off event occurs at T1
T2  = 5;    % dead time duration (µs)
T3  = 5;    % second pulse duration (µs) — switch-on event occurs at T1+T2
R_G = 5.1;  % gate resistance (Ohms)

sim_timeoffset = 10e-6 % LTspice allocates 10us before T1 begins

% Consistent MOSFET colours: blue / orange-red / green
colors = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]};


%% DPT Timing & Circuit Parameters

T1  = 35;   % first pulse duration (µs) - switch-off event occurs at T1
T2  = 5;    % dead time duration (µs)
T3  = 5;    % second pulse duration (µs) - switch-on event occurs at T1+T2
R_G = 5.1;  % gate resistance (Ohms)
timewindow = 1150;

% Consistent MOSFET colours: blue / orange-red / green
colors = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]};


%% Plot Control Flags
% Set to true/false to enable/disable each figure group

plt_include_simulation   = 0;   % overlay LTspice simulation series on all plots

plt_gate_current_full    = 1;   % "Full Sample Window - Gate Currents"
plt_gate_current_events  = 1;   % "Gate Current - Switch-Off/On"
plt_vgs_compare          = 1;   % "VGS Comparison - Switch-Off/On"
plt_raw_rogowski_events  = 0;   % "Raw Rogowski & Gate Currents - Switch-Off/On"
plt_corr_rogowski_events = 1;   % "Corrected Rogowski - Switch-Off/On"
plt_raw_rogowski_full    = 1;   % "Full Sample Window - Raw Rogowski Measurements"
plt_drain_current_full   = 1;   % "Full Sample Window - Drain Currents"
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


sim_path = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Bristol', 'grp-GRP Group 1002 - Documents', 'DPT Results', 'LTspice exports');
% sim_file = 'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_ANALYTICAL_MODEL_allsignals_prePWLtune.csv'
sim_file = 'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_ANALYTICAL_MODEL_allsignals.csv'

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

%% Section 4: SIMULATION import LTspice simulation data

% ltspice_csv = fullfile(sim_path, 'LTspiceExport_35-5-5_48VDC_36uH.csv');
ltspice_csv = fullfile(sim_path, sim_file);
sim_raw = readtable(ltspice_csv);

time_sim   = sim_raw{:,1} - sim_timeoffset;   % offset: shift LTspice timestamps by -10 us
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
I_load_sim = sim_raw{:,13};
IG1_sim   = sim_raw{:,14};
IG2_sim   = sim_raw{:,15};
IG3_sim   = sim_raw{:,16};


%% Section 4b: Simulation Event-Specific Alignment

avg_ns = 20e-9;   % averaging window width

VGS_emp_mean = (VGS1 + VGS2 + VGS3) / 3;
VGS_sim_mean = (VGS1_sim + VGS2_sim + VGS3_sim) / 3;

% Switch-Off: gate falls, event starts at T1
t_off = T1 * 1e-6;

mask_hi_e   = time >= t_off               & time     <= t_off + avg_ns;
mask_lo_e   = time >= t_off + 2e-6 - avg_ns & time   <= t_off + 2e-6;
V_hi_e      = mean(VGS_emp_mean(mask_hi_e));
V_lo_e      = mean(VGS_emp_mean(mask_lo_e));
thresh_e    = V_hi_e - 0.01 * (V_hi_e - V_lo_e);
mask_win_e  = time >= t_off & time <= t_off + 2e-6;
[~, ix]     = min(abs(VGS_emp_mean(mask_win_e) - thresh_e));
t_win       = time(mask_win_e);
t_emp_1pct_off = t_win(ix);

mask_hi_s   = time_sim >= t_off               & time_sim <= t_off + avg_ns;
mask_lo_s   = time_sim >= t_off + 2e-6 - avg_ns & time_sim <= t_off + 2e-6;
V_hi_s      = mean(VGS_sim_mean(mask_hi_s));
V_lo_s      = mean(VGS_sim_mean(mask_lo_s));
thresh_s    = V_hi_s - 0.01 * (V_hi_s - V_lo_s);
mask_win_s  = time_sim >= t_off & time_sim <= t_off + 2e-6;
[~, ix]     = min(abs(VGS_sim_mean(mask_win_s) - thresh_s));
t_win_s     = time_sim(mask_win_s);
t_sim_1pct_off = t_win_s(ix);

sim_timeoffset_off = t_sim_1pct_off - t_emp_1pct_off

% Switch-On: gate rises, event starts at T1+T2
t_on = (T1 + T2) * 1e-6;

mask_lo_e   = time >= t_on               & time     <= t_on + avg_ns;
mask_hi_e   = time >= t_on + 2e-6 - avg_ns & time   <= t_on + 2e-6;
V_lo_e      = mean(VGS_emp_mean(mask_lo_e));
V_hi_e      = mean(VGS_emp_mean(mask_hi_e));
thresh_e    = V_lo_e + 0.01 * (V_hi_e - V_lo_e);
mask_win_e  = time >= t_on & time <= t_on + 2e-6;
[~, ix]     = min(abs(VGS_emp_mean(mask_win_e) - thresh_e));
t_win       = time(mask_win_e);
t_emp_1pct_on = t_win(ix);

mask_lo_s   = time_sim >= t_on               & time_sim <= t_on + avg_ns;
mask_hi_s   = time_sim >= t_on + 2e-6 - avg_ns & time_sim <= t_on + 2e-6;
V_lo_s      = mean(VGS_sim_mean(mask_lo_s));
V_hi_s      = mean(VGS_sim_mean(mask_hi_s));
thresh_s    = V_lo_s + 0.01 * (V_hi_s - V_lo_s);
mask_win_s  = time_sim >= t_on & time_sim <= t_on + 2e-6;
[~, ix]     = min(abs(VGS_sim_mean(mask_win_s) - thresh_s));
t_win_s     = time_sim(mask_win_s);
t_sim_1pct_on = t_win_s(ix);

sim_timeoffset_on = t_sim_1pct_on - t_emp_1pct_on


%% Section 5: Gate Current Plots

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
                  time_sim - sim_timeoffset_off, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                  T1,     'Switch-Off', colors, timewindow, plt_include_simulation);
    plot_ig_event(time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, ...
                  time_sim - sim_timeoffset_on, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                  T1+T2,  'Switch-On',  colors, timewindow, plt_include_simulation);
end

if plt_vgs_compare
    plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                     time, VGS1_new, VGS2_new, VGS3_new, ...
                     time_sim - sim_timeoffset_off, VGS1_sim, VGS2_sim, VGS3_sim, ...
                     T1,    'Switch-Off', colors, timewindow, plt_include_simulation);
    plot_vgs_compare(time, VGS1,     VGS2,     VGS3, ...
                     time, VGS1_new, VGS2_new, VGS3_new, ...
                     time_sim - sim_timeoffset_on, VGS1_sim, VGS2_sim, VGS3_sim, ...
                     T1+T2, 'Switch-On',  colors, timewindow, plt_include_simulation);
end


%% Section 6: Rogowski Coil Correction [EDIT AS NEEDED]
% All signals share the same timebase - add or subtract I_G terms directly.
% Default: no correction applied.

I_total_corr = I_total;
I_23_corr    = I_23;
I_3_corr     = I_3;

if plt_raw_rogowski_events
    plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1,     'Switch-Off', colors);
    plot_raw_and_ig_event(time, I_total, I_23, I_3, I_G1, I_G2, I_G3, T1+T2,  'Switch-On',  colors);
end

if plt_corr_rogowski_events
    plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1,     'Switch-Off', colors);
    plot_corr_event(time, I_total_corr, I_23_corr, I_3_corr, T1+T2,  'Switch-On',  colors);
end


%% Section 6: Per-MOSFET Drain Current Deduction

I_D3 = I_3_corr;
I_D2 = I_23_corr    - I_D3;
I_D1 = I_total_corr - I_23_corr;


%% Section 8: Drain Current Plots

if plt_raw_rogowski_full
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
end

if plt_drain_current_full
    figure('Name', 'Full Sample Window — Drain Currents', 'NumberTitle', 'off');
    hold on;
    plot(time     * 1e6, I_D1,    '-',  'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 empirical');
    plot(time     * 1e6, I_D2,    '-',  'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 empirical');
    plot(time     * 1e6, I_D3,    '-',  'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 empirical');
    if plt_include_simulation
        plot(time_sim * 1e6, ID1_sim, '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
        plot(time_sim * 1e6, ID2_sim, '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
        plot(time_sim * 1e6, ID3_sim, '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    end
    hold off;
    xlabel('Time (\mus)');
    ylabel('Drain Current I_D (A)');
    title('Per-MOSFET Drain Currents — Full Sample Window');
    legend('Location', 'best');
    grid on;
end

if plt_drain_current_events
    plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
               time_sim - sim_timeoffset_off, ID1_sim, ID2_sim, ID3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
               T1,     'Switch-Off', colors, timewindow, plt_include_simulation);
    plot_event(time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
               time_sim - sim_timeoffset_on, ID1_sim, ID2_sim, ID3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
               T1+T2,  'Switch-On',  colors, timewindow, plt_include_simulation);
end

%% Section 9: Per MOSFET Switching Metric Measurement
% Applies the identical methodology used in DPT_SwitchingMetrics (%% Process
% Each T1) to each individual MOSFET at the fixed T1 = 35 µs defined above.
%
% Assumptions:
%   - VDS is identical across all three MOSFETs  →  VDS1 used throughout
%   - Per-MOSFET gate voltage references: VGS1, VGS2, VGS3
%   - Per-MOSFET drain currents:          I_D1, I_D2, I_D3
%   - Search windows: 1000 ns after each switching event (same as DPT_SwitchingMetrics)

VGS_ss_m9 = 15;   % nominal gate drive voltage (V) — local to this section

% Search windows
w_off1 = T1        * 1e-6;   w_off2 = (T1 + 1)      * 1e-6;
w_on1  = (T1 + T2) * 1e-6;   w_on2  = (T1 + T2 + 1) * 1e-6;

% VDS steady-state references — shared (VDS1 used for all MOSFETs)
VDS_off_ss_m9 = sample_at(time, VDS1, w_off2);   % VDS at T1 + 1 µs
VDS_on_ss_m9  = sample_at(time, VDS1, w_on1);    % VDS at T1 + T2

% Shared VDS threshold crossings — identical for all three MOSFETs
t_VDS_10up = find_crossing(time, VDS1, 0.10 * VDS_off_ss_m9, +1, w_off1, w_off2);
t_VDS_90up = find_crossing(time, VDS1, 0.90 * VDS_off_ss_m9, +1, w_off1, w_off2);
t_VDS_90dn = find_crossing(time, VDS1, 0.90 * VDS_on_ss_m9,  -1, w_on1,  w_on2);
t_VDS_10dn = find_crossing(time, VDS1, 0.10 * VDS_on_ss_m9,  -1, w_on1,  w_on2);

dv_dt_off_m9 = (0.90 - 0.10) * VDS_off_ss_m9 / ((t_VDS_90up - t_VDS_10up) * 1e9);
dv_dt_on_m9  = (0.90 - 0.10) * VDS_on_ss_m9  / ((t_VDS_10dn - t_VDS_90dn) * 1e9);

VGS_all = {VGS1, VGS2, VGS3};
I_D_all = {I_D1, I_D2, I_D3};

for k = 1:3
    VGSk = VGS_all{k};
    I_Dk = I_D_all{k};

    % Per-MOSFET steady-state current references
    I_ss_off = sample_at(time, I_Dk, w_off1);   % I_D at T1
    I_ss_on  = sample_at(time, I_Dk, w_on2);    % I_D at T1 + T2 + 1 µs

    % Turn-off: VGS and I_D crossings
    t_VGS_90 = find_crossing(time, VGSk, 0.90 * VGS_ss_m9, -1, w_off1, w_off2);
    t_I_90dn = find_crossing(time, I_Dk, 0.90 * I_ss_off,  -1, w_off1, w_off2);
    t_I_10dn = find_crossing(time, I_Dk, 0.10 * I_ss_off,  -1, w_off1, w_off2);

    t_d_off   = (t_VDS_10up - t_VGS_90)            * 1e9;   % ns
    t_f       = (t_VDS_90up - t_VDS_10up)           * 1e9;   % ns
    di_dt_off = (0.10 - 0.90) * I_ss_off / ((t_I_10dn - t_I_90dn) * 1e9);   % A/ns
    E_off     = energy_integral(time, VDS1, I_Dk, t_VDS_10up, t_I_10dn) * 1e6; % µJ

    % Turn-on: VGS and I_D crossings
    t_VGS_10 = find_crossing(time, VGSk, 0.10 * VGS_ss_m9, +1, w_on1, w_on2);
    t_I_10up = find_crossing(time, I_Dk, 0.10 * I_ss_on,   +1, w_on1, w_on2);
    t_I_90up = find_crossing(time, I_Dk, 0.90 * I_ss_on,   +1, w_on1, w_on2);

    t_d_on   = (t_VDS_90dn - t_VGS_10)           * 1e9;   % ns
    t_r      = (t_VDS_10dn - t_VDS_90dn)          * 1e9;   % ns
    di_dt_on = (0.90 - 0.10) * I_ss_on / ((t_I_90up - t_I_10up) * 1e9);    % A/ns
    E_on     = energy_integral(time, VDS1, I_Dk, t_I_10up, t_VDS_10dn) * 1e6; % µJ

    m9(k).MOSFET    = k;
    m9(k).I_ss_off  = I_ss_off;
    m9(k).I_ss_on   = I_ss_on;
    m9(k).t_VGS_90  = t_VGS_90;
    m9(k).t_VGS_10  = t_VGS_10;
    m9(k).t_I_90dn  = t_I_90dn;
    m9(k).t_I_10dn  = t_I_10dn;
    m9(k).t_I_10up  = t_I_10up;
    m9(k).t_I_90up  = t_I_90up;
    m9(k).t_d_off   = t_d_off;
    m9(k).t_f       = t_f;
    m9(k).E_off     = E_off;
    m9(k).di_dt_off = di_dt_off;
    m9(k).dv_dt_off = dv_dt_off_m9;   % shared — same for all MOSFETs
    m9(k).t_d_on    = t_d_on;
    m9(k).t_r       = t_r;
    m9(k).E_on      = E_on;
    m9(k).di_dt_on  = di_dt_on;
    m9(k).dv_dt_on  = dv_dt_on_m9;    % shared — same for all MOSFETs
end

T9 = table([m9.MOSFET]', [m9.I_ss_off]', [m9.I_ss_on]', ...
    [m9.t_d_off]', [m9.t_f]',  [m9.E_off]',  abs([m9.di_dt_off]'), [m9.dv_dt_off]', ...
    [m9.t_d_on]',  [m9.t_r]',  [m9.E_on]',   [m9.di_dt_on]',       [m9.dv_dt_on]', ...
    'VariableNames', {'MOSFET', 'I_ss_off_A', 'I_ss_on_A', ...
        't_d_off_ns', 't_f_ns',  'E_off_uJ',  'di_dt_off_Ans', 'dv_dt_off_Vns', ...
        't_d_on_ns',  't_r_ns',  'E_on_uJ',   'di_dt_on_Ans',  'dv_dt_on_Vns'});

fprintf('\n=== Per-MOSFET Switching Metrics  (T1=%d µs, T2=%d µs, 48 VDC) ===\n', T1, T2);
fprintf('    VDS_off_ss = %.1f V    VDS_on_ss = %.1f V\n', VDS_off_ss_m9, VDS_on_ss_m9);
disp(T9);

% --- Diagnostic plot: one 3×2 figure per MOSFET ---
for k = 1:3
    pts9.VDS_off_ss = VDS_off_ss_m9;
    pts9.VDS_on_ss  = VDS_on_ss_m9;
    pts9.t_VDS_10up = t_VDS_10up;
    pts9.t_VDS_90up = t_VDS_90up;
    pts9.t_VDS_90dn = t_VDS_90dn;
    pts9.t_VDS_10dn = t_VDS_10dn;
    pts9.I_ss_off   = m9(k).I_ss_off;
    pts9.I_ss_on    = m9(k).I_ss_on;
    pts9.t_VGS_90   = m9(k).t_VGS_90;
    pts9.t_VGS_10   = m9(k).t_VGS_10;
    pts9.t_I_90dn   = m9(k).t_I_90dn;
    pts9.t_I_10dn   = m9(k).t_I_10dn;
    pts9.t_I_10up   = m9(k).t_I_10up;
    pts9.t_I_90up   = m9(k).t_I_90up;
    pts9.t_d_off    = m9(k).t_d_off;
    pts9.t_f        = m9(k).t_f;
    pts9.dv_dt_off  = m9(k).dv_dt_off;
    pts9.di_dt_off  = m9(k).di_dt_off;
    pts9.t_d_on     = m9(k).t_d_on;
    pts9.t_r        = m9(k).t_r;
    pts9.dv_dt_on   = m9(k).dv_dt_on;
    pts9.di_dt_on   = m9(k).di_dt_on;
    pts9.E_off      = m9(k).E_off;
    pts9.E_on       = m9(k).E_on;
    plot_diagnostics(time, VGS_all{k}, I_D_all{k}, VDS1, T1, T2, VGS_ss_m9, pts9, ...
        sprintf('MOSFET %d  —  T_1 = %d µs, 48 VDC', k, T1));
end

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
                          t_sim,  VGS1_sim,  VGS2_sim,  VGS3_sim, ...
                          t_event_us, label, colors, timewindow, include_sim)
    t_start  = t_event_us * 1e-6;
    t_end    = t_start + timewindow * 1e-9;
    mk_orig  = t_orig >= t_start & t_orig <= t_end;
    mk_new   = t_new  >= t_start & t_new  <= t_end;
    mk_sim   = t_sim  >= t_start & t_sim  <= t_end;
    t_ns_o   = (t_orig(mk_orig) - t_start) * 1e9;
    t_ns_n   = (t_new(mk_new)   - t_start) * 1e9;
    t_ns_s   = (t_sim(mk_sim)   - t_start) * 1e9;

    vgs_orig = {VGS1_orig, VGS2_orig, VGS3_orig};
    vgs_new  = {VGS1_new,  VGS2_new,  VGS3_new};
    vgs_sim  = {VGS1_sim,  VGS2_sim,  VGS3_sim};
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
        if include_sim
            plot(t_ns_s, vgs_sim{k}(mk_sim),   ':',  'Color', colors{k}, ...
                 'LineWidth', 1.5, 'DisplayName', 'Simulation');
        end
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
                       t_event_us, label, colors, timewindow, include_sim)
    t_start  = t_event_us * 1e-6;
    t_end    = t_start + timewindow * 1e-9;

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
    if include_sim
        plot(t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
        plot(t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
        plot(t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    end
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
    if include_sim
        plot(t_ns_sim, IG1_sim(mask_sim),  '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
        plot(t_ns_sim, IG2_sim(mask_sim),  '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
        plot(t_ns_sim, IG3_sim(mask_sim),  '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    end
    hold off;
    xlabel('Time (ns)');
    ylabel('I_G (A)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end

function plot_event(t, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
                    t_sim, ID1_sim, ID2_sim, ID3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                    t_event_us, label, colors, timewindow, include_sim)
    t_start  = t_event_us * 1e-6;
    t_end    = t_start + timewindow * 1e-9;

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
    if include_sim
        plot(t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
        plot(t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
        plot(t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    end
    hold off;
    xlabel('Time (ns)');
    ylabel('V_{GS} (V)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;

    nexttile;
    hold on;
    plot(t_ns,     I_D1(mask),         '-',  'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
    plot(t_ns,     I_D2(mask),         '-',  'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
    plot(t_ns,     I_D3(mask),         '-',  'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
    if include_sim
        plot(t_ns_sim, ID1_sim(mask_sim),  '--', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 simulation');
        plot(t_ns_sim, ID2_sim(mask_sim),  '--', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 simulation');
        plot(t_ns_sim, ID3_sim(mask_sim),  '--', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 simulation');
    end
    hold off;
    ylabel('Drain Current I_D (A)');
    legend('Location', 'best');
    xlim([0, timewindow]);
    grid on;
end
