%% DPT Waveform Analysis — Paralleled MOSFET Current Imbalance
% Imports CH1–CH4 from 3 DPT test folders, deduces per-MOSFET drain currents
% via rogowski coil subtraction, and plots I_D, I_G and V_GS during switch-off
% and switch-on events over a 1150 ns window.

%% Timeseries alignment strategy
% DPT pulsetrain periods are identical between empirical and simulation, but absolute timestamps are not. 
% EMPIRICAL timestamps stay fixed, since they correlate to the physically meaningful logic trigger on a 3.6V T1 rising edge
% SIMULATION timestamps are different, and simulation doesn't model logic signal or the gate driver chip,
% hence simulation data will be aligned with 1% changes in the (average) device gate voltage swing, such that gate voltage waveforms are consistently timed

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROGRAM CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% DPT Timing & Circuit Parameters
T1  = 35;   % first pulse duration (µs) - switch-off event occurs at T1
T2  = 5;    % dead time duration (µs)
T3  = 5;    % second pulse duration (µs) - switch-on event occurs at T1+T2
R_G = 5.1;  % gate resistance (Ohms)

%% LTspice: define relevant model parameters
% Passive Component Values
gate_resistor          = R_G;     % [2, R_G] Ω  — gate input resistor
snubber_capacitor      = 99e-9;    % F  — strictly: decoupling not snubber
snubber_resistor       = 2.7;      % Ω  — strictly: decoupling not snubber
DC_link_capacitor      = 132e-6;   % F
% Parasitic Inductances
gate_inductance        = 20e-9;    % [20e-9, 100e-9] H 20n  — gate input inductance
power_loop_inductance  = 3.8e-9;   % H 3.8e-9
common_mode_inductance = 1e-9;     % H 1n
DC_link_inductance     = 5.65e-9;  % H 5.65n
branch_inductance      = 1e-9;    % H 12n — bus bar inductance between parallel branches
% Parasitic Resistances
rds_on = 1.1e-3; 
DC_link_resistance     = 150e-6;   % Ω  — bus bar resistance
C_DC_ESR               = 0.48e-3;  % Ω  — DC link capacitor ESR
C_snb_ESR              = 1e-3;     % Ω  — snubber capacitor ESR
branch_resistance      = 20e-6;    % Ω 20m  — bus bar resistance between parallel branches
solder_t               = [0.05e-3, 0.2e-3, 1e-3];  % m solder thickness
% PWL voltage source for DPT gate pulses
Trise = 80e-9;
Tfall = 70e-9;

%% Program Variables
sim_timeoffset = 10e-6; % LTspice allocates 10us before T1 begins
timewindow = 1150;
colors     = {[0 0.447 0.741], [0.850 0.325 0.098], [0.466 0.674 0.188]}; % MOSFET colours: blue / orange-red / green
sim_markers = {'o', '+', '*', '.', 'x', 'square', 'diamond'}; % markers cycled per sweep run — color identifies MOSFET, marker identifies parameter value


%% Plot Control Flags
% Set to true/false to enable/disable each figure group
plt_gate_current_full    = 1;   % "Full Sample Window - Gate Currents"
plt_gate_current_events  = 1;   % "Gate Current - Switch-Off/On"
plt_vgs_compare          = 1;   % "VGS Comparison - Switch-Off/On"
plt_raw_rogowski_events  = 0;   % "Raw Rogowski & Gate Currents - Switch-Off/On"
plt_corr_rogowski_events = 0;   % "Corrected Rogowski - Switch-Off/On"
plt_raw_rogowski_full    = 0;   % "Full Sample Window - Raw Rogowski Measurements"
plt_drain_current_full   = 1;   % "Full Sample Window - Drain Currents"
plt_dc_link_current_full = 1;   % "Full Sample Window - Total DC Link Current (sum I_D1..I_D3)"
plt_drain_current_events = 1;   % "Switch-Off/On" drain current + VGS


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Section 1: Paths and Folder Definitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% sim_path = fullfile(getenv('USERPROFILE'), 'OneDrive - University of Bristol', 'grp-GRP Group 1002 - Documents', 'DPT Results', 'LTspice exports');
% sim_file = 'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_ANALYTICAL_MODEL_allsignals_prePWLtune.csv'; % OBSOLETE dataset
% sim_file = 'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_ANALYTICAL_MODEL_allsignals.csv'; % OBSOLETE dataset

ltspice_exe = fullfile(getenv('LOCALAPPDATA'), 'Programs', 'ADI', 'LTspice', 'LTspice.exe'); % system location of executable
iscad_root  = fileparts(fileparts(mfilename('fullpath')));   % up from Matlab/ to ISCAD_GIP/ % system directory containing LTspice
%asc_file    = fullfile(iscad_root, 'LTspice parasitics', ... % system directory containing .asc LTspice model files
%              'InverterUnbalancedDPT', ...
%              'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_ANALYTICAL_MODEL.asc');
asc_file    = fullfile(iscad_root, 'LTspice parasitics', ... % system directory containing .asc LTspice model files
              'InverterUnbalancedDPT', ...
              'InverterUnbalancedDPT_35-5-5DPT_48VDC_POST_EMPIRICAL_MODEL.asc');
scratch_asc = strrep(asc_file,    '.asc', '_scratch.asc'); % working copy — MATLAB patches params here before running
scratch_raw = strrep(scratch_asc, '.asc', '.raw');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EMPIRICAL DATA IMPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION DATA IMPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Section 4: Build parameter struct and expand sweep

spice_params.gate_resistor          = gate_resistor;
spice_params.snubber_capacitor      = snubber_capacitor;
spice_params.snubber_resistor       = snubber_resistor;
spice_params.DC_link_capacitor      = DC_link_capacitor;
spice_params.gate_inductance        = gate_inductance;
spice_params.power_loop_inductance  = power_loop_inductance;
spice_params.common_mode_inductance = common_mode_inductance;
spice_params.DC_link_inductance     = DC_link_inductance;
spice_params.branch_inductance      = branch_inductance;
spice_params.DC_link_resistance     = DC_link_resistance;
spice_params.C_DC_ESR               = C_DC_ESR;
spice_params.C_snb_ESR              = C_snb_ESR;
spice_params.branch_resistance      = branch_resistance;
spice_params.Trise                  = Trise;
spice_params.Tfall                  = Tfall;
spice_params.rds_on                 = rds_on;
spice_params.solder_t               = solder_t;

param_list = expand_param_sweep(spice_params);
n_runs     = numel(param_list);
fprintf('Sweep: %d simulation run(s)\n', n_runs);


%% Section 4b: Run simulations and import results

VGS_emp_mean = (VGS1 + VGS2 + VGS3) / 3;

template     = struct('params', [], 'time_sim', [], ...
    'VGS1_sim', [], 'VGS2_sim', [], 'VGS3_sim', [], ...
    'VDS1_sim', [], 'VDS2_sim', [], 'VDS3_sim', [], ...
    'ID1_sim',  [], 'ID2_sim',  [], 'ID3_sim',  [], ...
    'IG1_sim',  [], 'IG2_sim',  [], 'IG3_sim',  [], ...
    'sim_timeoffset_off', [], 'sim_timeoffset_on', []);
sweep_results = repmat(template, 1, n_runs);

for k = 1:n_runs
    fprintf('  Run %d/%d ...\n', k, n_runs);

    p = param_list{k};
    write_ltspice_scratch(asc_file, scratch_asc, p);
    system(sprintf('"%s" -b -ascii "%s"', ltspice_exe, scratch_asc));
    sim = read_ltspice_raw(scratch_raw);

    time_sim = sim('time') - sim_timeoffset;
    VGS1_sim = sim('V(v_g_low1)') - sim('V(v_s_low1)');
    VGS2_sim = sim('V(v_g_low2)') - sim('V(v_s_low2)');
    VGS3_sim = sim('V(v_g_low3)') - sim('V(v_s_low3)');
    VGS_sim_mean = (VGS1_sim + VGS2_sim + VGS3_sim) / 3;

    [offset_off, offset_on] = compute_sim_alignment( ...
        time, VGS_emp_mean, time_sim, VGS_sim_mean, T1, T2);

    sweep_results(k).params             = p;
    sweep_results(k).time_sim           = time_sim;
    sweep_results(k).VGS1_sim           = VGS1_sim;
    sweep_results(k).VGS2_sim           = VGS2_sim;
    sweep_results(k).VGS3_sim           = VGS3_sim;
    sweep_results(k).VDS1_sim           = sim('V(phase)')    - sim('V(v_s_low1)');
    sweep_results(k).VDS2_sim           = sim('V(n016)')     - sim('V(v_s_low2)');
    sweep_results(k).VDS3_sim           = sim('V(n018)')     - sim('V(v_s_low3)');
    sweep_results(k).ID1_sim            = sim('I(L_cm1[L])');
    sweep_results(k).ID2_sim            = sim('I(L_cm2[L])');
    sweep_results(k).ID3_sim            = sim('I(L_cm3[L])');
    sweep_results(k).IG1_sim            = sim('I(R_g1)');
    sweep_results(k).IG2_sim            = sim('I(R_g2)');
    sweep_results(k).IG3_sim            = sim('I(R_g3)');
    sweep_results(k).sim_timeoffset_off = offset_off;
    sweep_results(k).sim_timeoffset_on  = offset_on;
end

disp('All simulations complete.')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Section 5: Gate Current Plots

if plt_gate_current_full
    figure('Name', 'Full Sample Window - Gate Currents', 'NumberTitle', 'off');
    hold on;
    plot(time * 1e6, I_G1, 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 1');
    plot(time * 1e6, I_G2, 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 2');
    plot(time * 1e6, I_G3, 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'MOSFET 3');
    hold off;
    xlabel('Time (\mus)'); ylabel('Gate Current I_G (A)');
    title('Per-MOSFET Gate Currents - Full Sample Window');
    legend('Location', 'best'); grid on;
end

if plt_gate_current_events
    ax_ig_off = []; ax_ig_on = [];
    for k = 1:n_runs
        r       = sweep_results(k);
        lmarker = sim_markers{mod(k-1, numel(sim_markers)) + 1};
        lbl     = make_sweep_label(sweep_results, k);
        ax_ig_off = plot_ig_event(ax_ig_off, ...
            time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, ...
            r.time_sim - r.sim_timeoffset_off, ...
            r.IG1_sim, r.IG2_sim, r.IG3_sim, r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1,    'Switch-Off', colors, timewindow, lmarker, lbl);
        ax_ig_on  = plot_ig_event(ax_ig_on, ...
            time, I_G1, I_G2, I_G3, VGS1_new, VGS2_new, VGS3_new, ...
            r.time_sim - r.sim_timeoffset_on, ...
            r.IG1_sim, r.IG2_sim, r.IG3_sim, r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1+T2, 'Switch-On',  colors, timewindow, lmarker, lbl);
    end
end

if plt_vgs_compare
    ax_vgs_off = []; ax_vgs_on = [];
    for k = 1:n_runs
        r       = sweep_results(k);
        lmarker = sim_markers{mod(k-1, numel(sim_markers)) + 1};
        lbl     = make_sweep_label(sweep_results, k);
        ax_vgs_off = plot_vgs_compare(ax_vgs_off, ...
            time, VGS1, VGS2, VGS3, time, VGS1_new, VGS2_new, VGS3_new, ...
            r.time_sim - r.sim_timeoffset_off, r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1,    'Switch-Off', colors, timewindow, lmarker, lbl);
        ax_vgs_on  = plot_vgs_compare(ax_vgs_on, ...
            time, VGS1, VGS2, VGS3, time, VGS1_new, VGS2_new, VGS3_new, ...
            r.time_sim - r.sim_timeoffset_on,  r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1+T2, 'Switch-On',  colors, timewindow, lmarker, lbl);
    end
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


%% Section 7: Per-MOSFET Drain Current Deduction

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
    xlabel('Time (\mus)'); ylabel('Current (A)');
    title('Raw Rogowski Coil Measurements — Full Sample Window');
    legend('Location', 'best'); grid on;
end

if plt_drain_current_full
    figure('Name', 'Full Sample Window — Drain Currents', 'NumberTitle', 'off');
    hold on;
    plot(time * 1e6, I_D1, '-', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'M1 empirical');
    plot(time * 1e6, I_D2, '-', 'Color', colors{2}, 'LineWidth', 1.2, 'DisplayName', 'M2 empirical');
    plot(time * 1e6, I_D3, '-', 'Color', colors{3}, 'LineWidth', 1.2, 'DisplayName', 'M3 empirical');
    for k = 1:n_runs
        r       = sweep_results(k);
        lmarker = sim_markers{mod(k-1, numel(sim_markers)) + 1};
        lbl     = make_sweep_label(sweep_results, k);
        mi      = round(linspace(1, numel(r.time_sim), min(20, numel(r.time_sim))));
        plot(r.time_sim * 1e6, r.ID1_sim, '--', 'Color', colors{1}, 'Marker', lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M1 ' lbl]);
        plot(r.time_sim * 1e6, r.ID2_sim, '--', 'Color', colors{2}, 'Marker', lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M2 ' lbl]);
        plot(r.time_sim * 1e6, r.ID3_sim, '--', 'Color', colors{3}, 'Marker', lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M3 ' lbl]);
    end
    hold off;
    xlabel('Time (\mus)'); ylabel('Drain Current I_D (A)');
    title('Per-MOSFET Drain Currents — Full Sample Window');
    legend('Location', 'best'); grid on;
end

if plt_dc_link_current_full
    figure('Name', 'Full Sample Window — Total DC Link Current', 'NumberTitle', 'off');
    hold on;
    plot(time * 1e6, I_D1 + I_D2 + I_D3, '-', 'Color', colors{1}, 'LineWidth', 1.2, 'DisplayName', 'Total empirical');
    for k = 1:n_runs
        r       = sweep_results(k);
        lmarker = sim_markers{mod(k-1, numel(sim_markers)) + 1};
        lbl     = make_sweep_label(sweep_results, k);
        mi      = round(linspace(1, numel(r.time_sim), min(20, numel(r.time_sim))));
        plot(r.time_sim * 1e6, r.ID1_sim + r.ID2_sim + r.ID3_sim, '--', 'Color', colors{2}, 'Marker', lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['Total ' lbl]);
    end
    hold off;
    xlabel('Time (\mus)'); ylabel('Total DC Link Current (A)');
    title('Total DC Link Current — Full Sample Window');
    legend('Location', 'best'); grid on;
end

if plt_drain_current_events
    ax_drain_off = []; ax_drain_on = [];
    for k = 1:n_runs
        r       = sweep_results(k);
        lmarker = sim_markers{mod(k-1, numel(sim_markers)) + 1};
        lbl     = make_sweep_label(sweep_results, k);
        ax_drain_off = plot_event(ax_drain_off, ...
            time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
            r.time_sim - r.sim_timeoffset_off, ...
            r.ID1_sim, r.ID2_sim, r.ID3_sim, r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1,    'Switch-Off', colors, timewindow, lmarker, lbl);
        ax_drain_on  = plot_event(ax_drain_on, ...
            time, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
            r.time_sim - r.sim_timeoffset_on, ...
            r.ID1_sim, r.ID2_sim, r.ID3_sim, r.VGS1_sim, r.VGS2_sim, r.VGS3_sim, ...
            T1+T2, 'Switch-On',  colors, timewindow, lmarker, lbl);
    end
end

%% -------------------------------------------------------------------------
%  Local Functions
%% -------------------------------------------------------------------------

function signals = read_ltspice_raw(raw_file)
    fid = fopen(raw_file, 'r');
    n_vars   = 0;
    n_points = 0;
    var_names = {};

    while true
        line = fgetl(fid);
        if ~ischar(line), break; end
        line = strtrim(line);
        if startsWith(line, 'No. Variables:')
            n_vars = sscanf(line, 'No. Variables: %d');
        elseif startsWith(line, 'No. Points:')
            n_points = sscanf(line, 'No. Points: %d');
        elseif startsWith(line, 'Variables:')
            for k = 1:n_vars
                vline = fgetl(fid);
                if ~ischar(vline)
                    error('read_ltspice_raw: EOF in Variables block (expected %d vars, got %d)', n_vars, k-1);
                end
                vline = strtrim(vline);
                parts = regexp(vline, '\t', 'split');
                var_names{k} = strtrim(parts{2}); %#ok<AGROW>
            end
        elseif startsWith(line, 'Values:')
            break;
        end
    end

    % Values block: point index + var0 on first line, each remaining var on its own indented line
    data = zeros(n_points, n_vars);
    for i = 1:n_points
        first = fgetl(fid);                          % "<idx>\t<val0>"
        parts = regexp(first, '\t', 'split');
        data(i, 1) = str2double(parts{end});
        for j = 2:n_vars
            data(i, j) = str2double(strtrim(fgetl(fid)));
        end
    end
    fclose(fid);

    vals = cell(1, n_vars);
    for k = 1:n_vars
        vals{k} = data(:, k);
    end
    signals = containers.Map(var_names, vals);
end

function [t, v] = load_csv(filepath)
    data = readmatrix(filepath, 'NumHeaderLines', 1, 'Delimiter', ',');
    t = data(:, 1);
    v = data(:, 2);
end

function param_list = expand_param_sweep(base_params)
% Returns a cell array of scalar param structs covering the Cartesian
% product of all vector-valued fields. Scalar fields are fixed across all runs.
    names        = fieldnames(base_params);
    swept_names  = {};
    swept_vals   = {};
    for k = 1:numel(names)
        v = base_params.(names{k});
        if ~isscalar(v)
            swept_names{end+1} = names{k};      %#ok<AGROW>
            swept_vals{end+1}  = v(:)';          %#ok<AGROW>
        end
    end

    if isempty(swept_names)
        param_list = {base_params};
        return;
    end

    grids = cell(1, numel(swept_names));
    [grids{:}] = ndgrid(swept_vals{:});
    n_runs     = numel(grids{1});
    param_list = cell(1, n_runs);
    for i = 1:n_runs
        p = base_params;
        for j = 1:numel(swept_names)
            p.(swept_names{j}) = grids{j}(i);
        end
        param_list{i} = p;
    end
end

function [offset_off, offset_on] = compute_sim_alignment( ...
        time, VGS_emp_mean, time_sim, VGS_sim_mean, T1, T2)
% Aligns simulation time to empirical using the 1%-of-swing crossing on the
% mean gate voltage, computed separately for the switch-off and switch-on events.
    avg_ns = 20e-9;

    % Switch-Off: gate falls at T1
    t_off      = T1 * 1e-6;
    mask_hi_e  = time >= t_off               & time     <= t_off + avg_ns;
    mask_lo_e  = time >= t_off + 2e-6 - avg_ns & time   <= t_off + 2e-6;
    V_hi_e     = mean(VGS_emp_mean(mask_hi_e));
    V_lo_e     = mean(VGS_emp_mean(mask_lo_e));
    thresh_e   = V_hi_e - 0.01 * (V_hi_e - V_lo_e);
    mask_win   = time >= t_off & time <= t_off + 2e-6;
    [~, ix]    = min(abs(VGS_emp_mean(mask_win) - thresh_e));
    t_win      = time(mask_win);
    t_emp_off  = t_win(ix);

    mask_hi_s  = time_sim >= t_off               & time_sim <= t_off + avg_ns;
    mask_lo_s  = time_sim >= t_off + 2e-6 - avg_ns & time_sim <= t_off + 2e-6;
    V_hi_s     = mean(VGS_sim_mean(mask_hi_s));
    V_lo_s     = mean(VGS_sim_mean(mask_lo_s));
    thresh_s   = V_hi_s - 0.01 * (V_hi_s - V_lo_s);
    mask_win   = time_sim >= t_off & time_sim <= t_off + 2e-6;
    [~, ix]    = min(abs(VGS_sim_mean(mask_win) - thresh_s));
    t_win_s    = time_sim(mask_win);
    t_sim_off  = t_win_s(ix);

    offset_off = t_sim_off - t_emp_off;

    % Switch-On: gate rises at T1+T2
    t_on       = (T1 + T2) * 1e-6;
    mask_lo_e  = time >= t_on               & time     <= t_on + avg_ns;
    mask_hi_e  = time >= t_on + 2e-6 - avg_ns & time   <= t_on + 2e-6;
    V_lo_e     = mean(VGS_emp_mean(mask_lo_e));
    V_hi_e     = mean(VGS_emp_mean(mask_hi_e));
    thresh_e   = V_lo_e + 0.01 * (V_hi_e - V_lo_e);
    mask_win   = time >= t_on & time <= t_on + 2e-6;
    [~, ix]    = min(abs(VGS_emp_mean(mask_win) - thresh_e));
    t_win      = time(mask_win);
    t_emp_on   = t_win(ix);

    mask_lo_s  = time_sim >= t_on               & time_sim <= t_on + avg_ns;
    mask_hi_s  = time_sim >= t_on + 2e-6 - avg_ns & time_sim <= t_on + 2e-6;
    V_lo_s     = mean(VGS_sim_mean(mask_lo_s));
    V_hi_s     = mean(VGS_sim_mean(mask_hi_s));
    thresh_s   = V_lo_s + 0.01 * (V_hi_s - V_lo_s);
    mask_win   = time_sim >= t_on & time_sim <= t_on + 2e-6;
    [~, ix]    = min(abs(VGS_sim_mean(mask_win) - thresh_s));
    t_win_s    = time_sim(mask_win);
    t_sim_on   = t_win_s(ix);

    offset_on  = t_sim_on - t_emp_on;
end

function write_ltspice_scratch(src_asc, dst_asc, params)
% Copy src_asc → dst_asc, replacing each active .param line whose name
% matches a field in params struct with the MATLAB numeric value.
% Only replaces lines beginning with '!.param' (active directives);
% commented ';.param' lines are left untouched.
% Special field: 'rds_on' — edits Rd in a scratch copy of the .lib and
% patches the .asc .include to reference it. Stays on _L0 model.
% Rd_new = rds_on - 0.47e-3  (0.47m accounts for Rs=234u + channel ~236u)
    txt   = fileread(src_asc);
    names = fieldnames(params);
    for k = 1:numel(names)
        name = names{k};
        if strcmp(name, 'rds_on'), continue; end   % handled below
        val_str = sprintf('%.6g', params.(name));
        pat = sprintf('(!?\\.param\\s+%s\\s*=\\s*)[^\\r\\n]+', name);
        txt = regexprep(txt, pat, ['$1' val_str]);
    end
    if isfield(params, 'rds_on')
        asc_dir     = fileparts(src_asc);
        orig_lib    = fullfile(asc_dir, '..', 'OptiMOS5_80V_LTSpice.lib');
        scratch_lib = fullfile(asc_dir, '..', 'OptiMOS5_80V_LTSpice_scratch.lib');
        Rd_new      = max(params.rds_on - 0.47e-3, 1e-6);   % clamp to 1 µΩ minimum
        Rd_str      = sprintf('%.6g', Rd_new);
        lib_txt     = fileread(orig_lib);
        % Replace Rd value only within IPTG011N08NM5_L0 (0.63m is unique to this device)
        lib_txt = regexprep(lib_txt, '(Rd\s+d1\s+d2\s+)0\.63m(\s+TC=6m)', ['$1' Rd_str '$2'], 'ignorecase');
        lib_fid = fopen(scratch_lib, 'w');
        fwrite(lib_fid, lib_txt, 'char');
        fclose(lib_fid);
        % Redirect .include in scratch .asc to scratch lib
        txt = strrep(txt, 'OptiMOS5_80V_LTSpice.lib', 'OptiMOS5_80V_LTSpice_scratch.lib');
    end
    fid = fopen(dst_asc, 'w');
    fwrite(fid, txt, 'char');
    fclose(fid);
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

function ax = plot_vgs_compare(ax_in, t_orig, VGS1_orig, VGS2_orig, VGS3_orig, ...
                               t_new,  VGS1_new,  VGS2_new,  VGS3_new, ...
                               t_sim,  VGS1_sim,  VGS2_sim,  VGS3_sim, ...
                               t_event_us, label, colors, timewindow, sim_lmarker, sim_label)
    t_start = t_event_us * 1e-6;
    mk_orig = t_orig >= t_start & t_orig <= t_start + timewindow * 1e-9;
    mk_new  = t_new  >= t_start & t_new  <= t_start + timewindow * 1e-9;
    mk_sim  = t_sim  >= t_start & t_sim  <= t_start + timewindow * 1e-9;
    t_ns_o  = (t_orig(mk_orig) - t_start) * 1e9;
    t_ns_n  = (t_new(mk_new)   - t_start) * 1e9;
    t_ns_s  = (t_sim(mk_sim)   - t_start) * 1e9;
    vgs_orig = {VGS1_orig, VGS2_orig, VGS3_orig};
    vgs_new  = {VGS1_new,  VGS2_new,  VGS3_new};
    vgs_sim  = {VGS1_sim,  VGS2_sim,  VGS3_sim};
    mosfet_labels = {'MOSFET 1', 'MOSFET 2', 'MOSFET 3'};

    if isempty(ax_in)
        figure('Name', ['VGS Comparison — ' label], 'NumberTitle', 'off');
        tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, ['V_{GS} Comparison — ' label], 'FontWeight', 'bold');
        for m = 1:3
            ax(m) = nexttile; hold on; %#ok<AGROW>
            plot(t_ns_o, vgs_orig{m}(mk_orig), '--', 'Color', colors{m}, 'LineWidth', 1.2, 'DisplayName', 'Original');
            plot(t_ns_n, vgs_new{m}(mk_new),   '-',  'Color', colors{m}, 'LineWidth', 1.5, 'DisplayName', 'Re-measured');
            title(mosfet_labels{m}); ylabel('V_{GS} (V)');
            legend('Location', 'best'); xlim([0, timewindow]); grid on;
            if m == 3, xlabel('Time (ns)'); end
        end
    else
        ax = ax_in;
    end

    mi = round(linspace(1, numel(t_ns_s), min(10, numel(t_ns_s))));
    for m = 1:3
        plot(ax(m), t_ns_s, vgs_sim{m}(mk_sim), '--', 'Color', colors{m}, ...
             'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.5, 'DisplayName', ['Sim: ' sim_label]);
        legend(ax(m), 'Location', 'best');
    end
end

function ax = plot_ig_event(ax_in, t, I_G1, I_G2, I_G3, VGS1, VGS2, VGS3, ...
                            t_sim, IG1_sim, IG2_sim, IG3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                            t_event_us, label, colors, timewindow, sim_lmarker, sim_label)
    t_start  = t_event_us * 1e-6;
    mask     = t     >= t_start & t     <= t_start + timewindow * 1e-9;
    mask_sim = t_sim >= t_start & t_sim <= t_start + timewindow * 1e-9;
    t_ns     = (t(mask)         - t_start) * 1e9;
    t_ns_sim = (t_sim(mask_sim) - t_start) * 1e9;

    if isempty(ax_in)
        figure('Name', ['Gate V/I - ' label], 'NumberTitle', 'off');
        tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, ['Gate V/I - ' label], 'FontWeight', 'bold');
        ax(1) = nexttile; hold on;
        plot(t_ns, VGS1(mask), '-', 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
        plot(t_ns, VGS2(mask), '-', 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
        plot(t_ns, VGS3(mask), '-', 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
        ylabel('V_{GS} (V)'); legend('Location', 'best'); xlim([0, timewindow]); grid on;
        ax(2) = nexttile; hold on;
        plot(t_ns, I_G1(mask), '-', 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
        plot(t_ns, I_G2(mask), '-', 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
        plot(t_ns, I_G3(mask), '-', 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
        xlabel('Time (ns)'); ylabel('I_G (A)'); legend('Location', 'best'); xlim([0, timewindow]); grid on;
    else
        ax = ax_in;
    end

    mi = round(linspace(1, numel(t_ns_sim), min(10, numel(t_ns_sim))));
    plot(ax(1), t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M1 ' sim_label]);
    plot(ax(1), t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M2 ' sim_label]);
    plot(ax(1), t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M3 ' sim_label]);
    legend(ax(1), 'Location', 'best');
    plot(ax(2), t_ns_sim, IG1_sim(mask_sim), '--', 'Color', colors{1}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M1 ' sim_label]);
    plot(ax(2), t_ns_sim, IG2_sim(mask_sim), '--', 'Color', colors{2}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M2 ' sim_label]);
    plot(ax(2), t_ns_sim, IG3_sim(mask_sim), '--', 'Color', colors{3}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M3 ' sim_label]);
    legend(ax(2), 'Location', 'best');
end

function ax = plot_event(ax_in, t, I_D1, I_D2, I_D3, VGS1, VGS2, VGS3, ...
                         t_sim, ID1_sim, ID2_sim, ID3_sim, VGS1_sim, VGS2_sim, VGS3_sim, ...
                         t_event_us, label, colors, timewindow, sim_lmarker, sim_label)
    t_start  = t_event_us * 1e-6;
    mask     = t     >= t_start & t     <= t_start + timewindow * 1e-9;
    mask_sim = t_sim >= t_start & t_sim <= t_start + timewindow * 1e-9;
    t_ns     = (t(mask)         - t_start) * 1e9;
    t_ns_sim = (t_sim(mask_sim) - t_start) * 1e9;

    if isempty(ax_in)
        figure('Name', label, 'NumberTitle', 'off');
        tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, label, 'FontWeight', 'bold');
        ax(1) = nexttile; hold on;
        plot(t_ns, VGS1(mask), '-', 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
        plot(t_ns, VGS2(mask), '-', 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
        plot(t_ns, VGS3(mask), '-', 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
        ylabel('V_{GS} (V)'); legend('Location', 'best'); xlim([0, timewindow]); grid on;
        ax(2) = nexttile; hold on;
        plot(t_ns, I_D1(mask), '-', 'Color', colors{1}, 'LineWidth', 1.5, 'DisplayName', 'M1 empirical');
        plot(t_ns, I_D2(mask), '-', 'Color', colors{2}, 'LineWidth', 1.5, 'DisplayName', 'M2 empirical');
        plot(t_ns, I_D3(mask), '-', 'Color', colors{3}, 'LineWidth', 1.5, 'DisplayName', 'M3 empirical');
        xlabel('Time (ns)'); ylabel('Drain Current I_D (A)'); legend('Location', 'best'); xlim([0, timewindow]); grid on;
    else
        ax = ax_in;
    end

    mi = round(linspace(1, numel(t_ns_sim), min(10, numel(t_ns_sim))));
    plot(ax(1), t_ns_sim, VGS1_sim(mask_sim), '--', 'Color', colors{1}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M1 ' sim_label]);
    plot(ax(1), t_ns_sim, VGS2_sim(mask_sim), '--', 'Color', colors{2}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M2 ' sim_label]);
    plot(ax(1), t_ns_sim, VGS3_sim(mask_sim), '--', 'Color', colors{3}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M3 ' sim_label]);
    legend(ax(1), 'Location', 'best');
    plot(ax(2), t_ns_sim, ID1_sim(mask_sim), '--', 'Color', colors{1}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M1 ' sim_label]);
    plot(ax(2), t_ns_sim, ID2_sim(mask_sim), '--', 'Color', colors{2}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M2 ' sim_label]);
    plot(ax(2), t_ns_sim, ID3_sim(mask_sim), '--', 'Color', colors{3}, 'Marker', sim_lmarker, 'MarkerIndices', mi, 'LineWidth', 1.2, 'DisplayName', ['M3 ' sim_label]);
    legend(ax(2), 'Location', 'best');
end

function lbl = make_sweep_label(sweep_results, k)
% Returns a label showing the value of every swept parameter for run k.
% A parameter is swept if it takes more than one unique value across all runs.
% Returns 'sim' when there are no swept parameters.
    if isscalar(sweep_results)
        lbl = 'sim';
        return;
    end
    names = fieldnames(sweep_results(1).params);
    parts = {};
    for i = 1:numel(names)
        vals = arrayfun(@(r) r.params.(names{i}), sweep_results);
        if numel(unique(vals)) > 1
            parts{end+1} = sprintf('%s=%s', names{i}, format_si(sweep_results(k).params.(names{i}))); %#ok<AGROW>
        end
    end
    if isempty(parts)
        lbl = 'sim';
    else
        lbl = strjoin(parts, ', ');
    end
end

function s = format_si(v)
% Formats a scalar numeric value with an SI suffix (m, u, n, p).
    av = abs(v);
    if av == 0,      s = '0';                         return; end
    if av >= 1,      s = sprintf('%.4g',  v);          return; end
    if av >= 1e-3,   s = sprintf('%.4gm', v * 1e3);   return; end
    if av >= 1e-6,   s = sprintf('%.4gu', v * 1e6);   return; end
    if av >= 1e-9,   s = sprintf('%.4gn', v * 1e9);   return; end
    if av >= 1e-12,  s = sprintf('%.4gp', v * 1e12);  return; end
    s = sprintf('%.4g', v);
end
