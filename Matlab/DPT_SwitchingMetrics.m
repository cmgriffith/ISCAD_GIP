%% DPT Switching Metrics — T1 Parameter Sweep
% Extracts turn-off and turn-on switching characteristics from a T1 sweep at
% 48 VDC, Branch 3. Produces a diagnostic plot for T1=10 µs then tabulates
% metrics across all T1 values.
%
% Channels (Branch 3 subfolder): CH2 = VGS, CH3 = I_Dtotal, CH4 = VDS
%
% Steady-state references per test:
%   VDS_off_ss  = VDS  at T1 + 1 µs         (post turn-off settling)
%   I_ss_off    = I_D  at T1                 (pre turn-off)
%   VDS_on_ss   = VDS  at T1 + T2            (pre turn-on)
%   I_ss_on     = I_D  at T1 + T2 + 1 µs    (post turn-on settling)

%% Parameters

T2     = 5;     % dead time (µs)
VGS_ss = 15;    % assumed steady-state V_GS (V)

T1_values = [10, 20, 30, 40, 50, 70];   % µs

base_path = 'C:\Users\ossia\OneDrive - University of Bristol\grp-GRP Group 1002 - Documents\DPT Results';

folder_names = { ...
    '10-5-10DPT_48VDC_logic_VGS_IDStotal_VDS', ...
    '20-5-5DPT_48VDC_logic_VGS_IDStotal_VDS',  ...
    '30-5-5DPT_48VDC_logic_VGS_IDStotal_VDS',  ...
    '40-5-5DPT_48VDC_logic_VGS_IDStotal_VDS',  ...
    '50-5-5DPT_48VDC_logic_VGS_IDStotal_VDS',  ...
    '70-5-5DPT_48VDC_logic_VGS_IDStotal_VDS'};

%% Process Each T1

for i = 1:length(T1_values)
    T1 = T1_values(i);

    [time, VGS, I_D, VDS] = load_branch3(base_path, folder_names{i});

    % --- Steady-state reference values (nearest sample) ---
    VDS_off_ss = sample_at(time, VDS, (T1 + 1)     * 1e-6);
    I_ss_off   = sample_at(time, I_D,  T1           * 1e-6);
    VDS_on_ss  = sample_at(time, VDS, (T1 + T2)     * 1e-6);
    I_ss_on    = sample_at(time, I_D, (T1 + T2 + 1) * 1e-6);

    % --- Turn-off threshold crossings (search within 1000 ns after T1) ---
    w1 = T1       * 1e-6;
    w2 = (T1 + 1) * 1e-6;

    t_VGS_90   = find_crossing(time, VGS, 0.90 * VGS_ss,    -1, w1, w2);
    t_VDS_10up = find_crossing(time, VDS, 0.10 * VDS_off_ss, +1, w1, w2);
    t_VDS_90up = find_crossing(time, VDS, 0.90 * VDS_off_ss, +1, w1, w2);
    t_I_90dn   = find_crossing(time, I_D, 0.90 * I_ss_off,  -1, w1, w2);
    t_I_10dn   = find_crossing(time, I_D, 0.10 * I_ss_off,  -1, w1, w2);

    t_d_off   = (t_VDS_10up - t_VGS_90)              * 1e9;   % ns
    t_f       = (t_VDS_90up - t_VDS_10up)             * 1e9;   % ns
    dv_dt_off = (0.90 - 0.10) * VDS_off_ss            / ((t_VDS_90up - t_VDS_10up) * 1e9);  % V/ns
    di_dt_off = (0.10 - 0.90) * I_ss_off              / ((t_I_10dn   - t_I_90dn)   * 1e9);  % A/ns
    E_off     = energy_integral(time, VDS, I_D, t_VDS_10up, t_I_10dn);                        % J

    % --- Turn-on threshold crossings (search within 1000 ns after T1+T2) ---
    w3 = (T1 + T2)     * 1e-6;
    w4 = (T1 + T2 + 1) * 1e-6;

    t_VGS_10   = find_crossing(time, VGS, 0.10 * VGS_ss,    +1, w3, w4);
    t_VDS_90dn = find_crossing(time, VDS, 0.90 * VDS_on_ss, -1, w3, w4);
    t_VDS_10dn = find_crossing(time, VDS, 0.10 * VDS_on_ss, -1, w3, w4);
    t_I_10up   = find_crossing(time, I_D, 0.10 * I_ss_on,   +1, w3, w4);
    t_I_90up   = find_crossing(time, I_D, 0.90 * I_ss_on,   +1, w3, w4);

    t_d_on   = (t_VDS_90dn - t_VGS_10)              * 1e9;
    t_r      = (t_VDS_10dn - t_VDS_90dn)             * 1e9;
    dv_dt_on = (0.90 - 0.10) * VDS_on_ss             / ((t_VDS_10dn - t_VDS_90dn) * 1e9);
    di_dt_on = (0.90 - 0.10) * I_ss_on               / ((t_I_90up   - t_I_10up)   * 1e9);
    E_on     = energy_integral(time, VDS, I_D, t_I_10up, t_VDS_10dn);

    % --- Store ---
    results(i).T1         = T1;
    results(i).VDS_off_ss = VDS_off_ss;
    results(i).I_ss_off   = I_ss_off;
    results(i).VDS_on_ss  = VDS_on_ss;
    results(i).I_ss_on    = I_ss_on;
    results(i).t_d_off    = t_d_off;
    results(i).t_f        = t_f;
    results(i).E_off      = E_off * 1e6;    % µJ
    results(i).di_dt_off  = di_dt_off;
    results(i).dv_dt_off  = dv_dt_off;
    results(i).t_d_on     = t_d_on;
    results(i).t_r        = t_r;
    results(i).E_on       = E_on  * 1e6;
    results(i).di_dt_on   = di_dt_on;
    results(i).dv_dt_on   = dv_dt_on;

    % --- Diagnostic plot for T1 = 10 µs only ---
    if T1 == 40
        pts = struct( ...
            'VDS_off_ss', VDS_off_ss, 'I_ss_off',  I_ss_off, ...
            'VDS_on_ss',  VDS_on_ss,  'I_ss_on',   I_ss_on,  ...
            't_VGS_90',   t_VGS_90,   't_VDS_10up', t_VDS_10up, ...
            't_VDS_90up', t_VDS_90up, 't_I_90dn',   t_I_90dn, ...
            't_I_10dn',   t_I_10dn,   't_VGS_10',   t_VGS_10, ...
            't_VDS_90dn', t_VDS_90dn, 't_VDS_10dn', t_VDS_10dn, ...
            't_I_10up',   t_I_10up,   't_I_90up',   t_I_90up,  ...
            't_d_off',    t_d_off,    't_f',         t_f,       ...
            'E_off',      E_off*1e6,  'di_dt_off',   di_dt_off, ...
            'dv_dt_off',  dv_dt_off,  't_d_on',      t_d_on,    ...
            't_r',        t_r,        'E_on',        E_on*1e6,  ...
            'di_dt_on',   di_dt_on,   'dv_dt_on',    dv_dt_on);
        plot_diagnostics(time, VGS, I_D, VDS, T1, T2, VGS_ss, pts);
    end
end

%% Results Table

T = table( ...
    [results.T1]',         [results.VDS_off_ss]', [results.I_ss_off]', ...
    [results.VDS_on_ss]',  [results.I_ss_on]',    ...
    [results.t_d_off]',    [results.t_f]',         [results.E_off]',   ...
    [results.di_dt_off]',  [results.dv_dt_off]',   ...
    [results.t_d_on]',     [results.t_r]',          [results.E_on]',   ...
    [results.di_dt_on]',   [results.dv_dt_on]',     ...
    'VariableNames', { ...
        'T1_us',       'VDS_off_ss_V',  'I_ss_off_A', ...
        'VDS_on_ss_V', 'I_ss_on_A',     ...
        't_d_off_ns',  't_f_ns',        'E_off_uJ',   ...
        'di_dt_off_Ans','dv_dt_off_Vns', ...
        't_d_on_ns',   't_r_ns',        'E_on_uJ',    ...
        'di_dt_on_Ans','dv_dt_on_Vns'});

disp(T);

%% Summary Plots — Switching KPIs vs Drain Current

I_off_v  = [results.I_ss_off]';
I_on_v   = [results.I_ss_on]';
E_off_v  = [results.E_off]';          % µJ
E_on_v   = [results.E_on]';           % µJ
t_d_off_v = [results.t_d_off]';       % ns
t_f_v     = [results.t_f]';           % ns
t_d_on_v  = [results.t_d_on]';        % ns
t_r_v     = [results.t_r]';           % ns
dv_dt_off_v = [results.dv_dt_off]';   % V/ns
dv_dt_on_v  = [results.dv_dt_on]';    % V/ns
di_dt_off_v = abs([results.di_dt_off]');  % A/ns (magnitude — current falls)
di_dt_on_v  = [results.di_dt_on]';    % A/ns

c_off = [0.20 0.40 0.80];   % blue  — turn-off series
c_on  = [0.80 0.30 0.10];   % red   — turn-on series
c_tot = [0.50 0.00 0.50];   % purple — totals
lw = 1.5;  ms = 7;

figure('Name', 'DPT Parameter Sweep — KPI Summary', 'NumberTitle', 'off', ...
       'Position', [80 80 1100 860]);
tl2 = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, sprintf('Switching KPIs vs Drain Current  (48 V_{DC}, T_2 = %d µs)', T2), ...
      'FontWeight', 'bold');

% (1,1) Switching energies E_off and E_on
nexttile(1);
hold on;
plot(I_off_v, E_off_v, 'o-', 'Color', c_off, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_off);
plot(I_on_v,  E_on_v,  's-', 'Color', c_on,  'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_on);
hold off;
xlabel('I_D (A)');  ylabel('Energy (µJ)');
legend('E_{off}', 'E_{on}', 'Location', 'northwest');
title('Switching Energies');
grid on;

% (1,2) Total switching loss
nexttile(2);
plot(I_off_v, E_off_v + E_on_v, 'd-', 'Color', c_tot, 'LineWidth', lw, ...
     'MarkerSize', ms, 'MarkerFaceColor', c_tot);
xlabel('I_D (A)');  ylabel('Energy (µJ)');
title('Total Switching Loss  (E_{off} + E_{on})');
grid on;

% (2,1) Turn-off timing: t_d(off) and t_f
nexttile(3);
hold on;
plot(I_off_v, t_d_off_v, 'o-', 'Color', c_off,      'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_off);
plot(I_off_v, t_f_v,     's-', 'Color', c_off*0.55,  'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_off*0.55);
hold off;
xlabel('I_D (A)');  ylabel('Time (ns)');
legend('t_{d(off)}', 't_f', 'Location', 'best');
title('Turn-Off Timing');
grid on;

% (2,2) Turn-on timing: t_d(on) and t_r
nexttile(4);
hold on;
plot(I_on_v, t_d_on_v, 'o-', 'Color', c_on,      'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_on);
plot(I_on_v, t_r_v,    's-', 'Color', c_on*0.55,  'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_on*0.55);
hold off;
xlabel('I_D (A)');  ylabel('Time (ns)');
legend('t_{d(on)}', 't_r', 'Location', 'best');
title('Turn-On Timing');
grid on;

% (3,1) Voltage slew rates dv/dt
nexttile(5);
hold on;
plot(I_off_v, dv_dt_off_v, 'o-', 'Color', c_off, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_off);
plot(I_on_v,  dv_dt_on_v,  's-', 'Color', c_on,  'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_on);
hold off;
xlabel('I_D (A)');  ylabel('dV/dt (V/ns)');
legend('dv/dt_{off}', 'dv/dt_{on}', 'Location', 'best');
title('Voltage Slew Rates');
grid on;

% (3,2) Current slew rates di/dt
nexttile(6);
hold on;
plot(I_off_v, di_dt_off_v, 'o-', 'Color', c_off, 'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_off);
plot(I_on_v,  di_dt_on_v,  's-', 'Color', c_on,  'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', c_on);
hold off;
xlabel('I_D (A)');  ylabel('dI/dt (A/ns)');
legend('di/dt_{off}', 'di/dt_{on}', 'Location', 'best');
title('Current Slew Rates');
grid on;

%% -------------------------------------------------------------------------
%  Local Functions
%  Shared utilities (sample_at, find_crossing, energy_integral,
%  align_yyaxis_zeros) live as standalone .m files in this folder.
%% -------------------------------------------------------------------------

function [time, VGS, I_D, VDS] = load_branch3(base_path, folder)
    sub  = fullfile(base_path, folder, 'Branch 3');
    rd   = @(ch) readmatrix(fullfile(sub, sprintf('CH%d.CSV', ch)), ...
                            'NumHeaderLines', 1, 'Delimiter', ',');
    d2   = rd(2);  time = d2(:,1);  VGS = d2(:,2);
    d3   = rd(3);  I_D  = d3(:,2);
    d4   = rd(4);  VDS  = d4(:,2);
end

function plot_diagnostics(time, VGS, I_D, VDS, T1, T2, VGS_ss, pts)
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
    c_vgs = [0.15 0.55 0.15];
    c_pwr = [0.50 0.00 0.50];

    figure('Name', sprintf('Diagnostic — T_1 = %d µs', T1), ...
           'NumberTitle', 'off', 'Position', [40 40 1200 900]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Switching Metric Diagnostic  —  T_1 = %d µs', T1), ...
          'FontWeight', 'bold');

    % ---- Tile 1: Turn-Off  VDS (left y) + I_D (right y) ----
    nexttile(1);
    yyaxis left;  hold on;
    plot(ns_off(time(mk_off)), VDS(mk_off), 'Color', c_vds, 'LineWidth', 1.5);
    yline(0.10*pts.VDS_off_ss, '--', 'Color', c_vds, 'LineWidth', 0.8);
    yline(0.90*pts.VDS_off_ss, '--', 'Color', c_vds, 'LineWidth', 0.8);
    if ~isnan(pts.t_VDS_10up)
        scatter(ns_off(pts.t_VDS_10up), sample_at(time, VDS, pts.t_VDS_10up), 60, c_vds, 'filled');
    end
    if ~isnan(pts.t_VDS_90up)
        scatter(ns_off(pts.t_VDS_90up), sample_at(time, VDS, pts.t_VDS_90up), 60, c_vds, 'filled');
    end
    ylabel('V_{DS} (V)');
    yyaxis right;  hold on;
    plot(ns_off(time(mk_off)), I_D(mk_off), 'Color', c_id, 'LineWidth', 1.5);
    yline(0.90*pts.I_ss_off, '--', 'Color', c_id, 'LineWidth', 0.8);
    yline(0.10*pts.I_ss_off, '--', 'Color', c_id, 'LineWidth', 0.8);
    if ~isnan(pts.t_I_90dn)
        scatter(ns_off(pts.t_I_90dn), sample_at(time, I_D, pts.t_I_90dn), 60, c_id, 'filled');
    end
    if ~isnan(pts.t_I_10dn)
        scatter(ns_off(pts.t_I_10dn), sample_at(time, I_D, pts.t_I_10dn), 60, c_id, 'filled');
    end
    ylabel('I_D (A)');
    ax = gca;  ax.YAxis(1).Color = c_vds;  ax.YAxis(2).Color = c_id;
    xlim(xlims);  grid on;
    align_yyaxis_zeros(ax);
    title(sprintf('Turn-Off:  V_{DS} (blue, left) & I_D (red, right)\nt_d(off)=%.0f ns  t_f=%.0f ns  dv/dt=%.3f V/ns  di/dt=%.3f A/ns', ...
                  pts.t_d_off, pts.t_f, pts.dv_dt_off, pts.di_dt_off));

    % ---- Tile 2: Turn-On  VDS (left y) + I_D (right y) ----
    nexttile(2);
    yyaxis left;  hold on;
    plot(ns_on(time(mk_on)), VDS(mk_on), 'Color', c_vds, 'LineWidth', 1.5);
    yline(0.90*pts.VDS_on_ss, '--', 'Color', c_vds, 'LineWidth', 0.8);
    yline(0.10*pts.VDS_on_ss, '--', 'Color', c_vds, 'LineWidth', 0.8);
    if ~isnan(pts.t_VDS_90dn)
        scatter(ns_on(pts.t_VDS_90dn), sample_at(time, VDS, pts.t_VDS_90dn), 60, c_vds, 'filled');
    end
    if ~isnan(pts.t_VDS_10dn)
        scatter(ns_on(pts.t_VDS_10dn), sample_at(time, VDS, pts.t_VDS_10dn), 60, c_vds, 'filled');
    end
    ylabel('V_{DS} (V)');
    yyaxis right;  hold on;
    plot(ns_on(time(mk_on)), I_D(mk_on), 'Color', c_id, 'LineWidth', 1.5);
    yline(0.10*pts.I_ss_on, '--', 'Color', c_id, 'LineWidth', 0.8);
    yline(0.90*pts.I_ss_on, '--', 'Color', c_id, 'LineWidth', 0.8);
    if ~isnan(pts.t_I_10up)
        scatter(ns_on(pts.t_I_10up), sample_at(time, I_D, pts.t_I_10up), 60, c_id, 'filled');
    end
    if ~isnan(pts.t_I_90up)
        scatter(ns_on(pts.t_I_90up), sample_at(time, I_D, pts.t_I_90up), 60, c_id, 'filled');
    end
    ylabel('I_D (A)');
    ax = gca;  ax.YAxis(1).Color = c_vds;  ax.YAxis(2).Color = c_id;
    xlim(xlims);  grid on;
    align_yyaxis_zeros(ax);
    title(sprintf('Turn-On:  V_{DS} (blue, left) & I_D (red, right)\nt_d(on)=%.0f ns  t_r=%.0f ns  dv/dt=%.3f V/ns  di/dt=%.3f A/ns', ...
                  pts.t_d_on, pts.t_r, pts.dv_dt_on, pts.di_dt_on));

    % ---- Tile 3: Turn-Off VGS ----
    nexttile(3);
    hold on;
    plot(ns_off(time(mk_off)), VGS(mk_off), 'Color', c_vgs, 'LineWidth', 1.5);
    yline(0.90*VGS_ss, '--', 'Color', c_vgs, 'LineWidth', 0.8);
    if ~isnan(pts.t_VGS_90)
        scatter(ns_off(pts.t_VGS_90), sample_at(time, VGS, pts.t_VGS_90), 60, c_vgs, 'filled');
    end
    hold off;
    ylabel('V_{GS} (V)');  xlabel('Time (ns)');
    xlim(xlims);  grid on;
    title(sprintf('Turn-Off: V_{GS}\nt_d(off) = %.0f ns', pts.t_d_off));

    % ---- Tile 4: Turn-On VGS ----
    nexttile(4);
    hold on;
    plot(ns_on(time(mk_on)), VGS(mk_on), 'Color', c_vgs, 'LineWidth', 1.5);
    yline(0.10*VGS_ss, '--', 'Color', c_vgs, 'LineWidth', 0.8);
    if ~isnan(pts.t_VGS_10)
        scatter(ns_on(pts.t_VGS_10), sample_at(time, VGS, pts.t_VGS_10), 60, c_vgs, 'filled');
    end
    hold off;
    ylabel('V_{GS} (V)');  xlabel('Time (ns)');
    xlim(xlims);  grid on;
    title(sprintf('Turn-On: V_{GS}\nt_d(on) = %.0f ns', pts.t_d_on));

    % ---- Tile 5: Turn-Off Power + E_off fill ----
    nexttile(5);
    hold on;
    if ~isnan(pts.t_VDS_10up) && ~isnan(pts.t_I_10dn)
        mk_int = time >= pts.t_VDS_10up & time <= pts.t_I_10dn;
        t_int  = ns_off(time(mk_int));
        p_int  = VDS(mk_int) .* I_D(mk_int);
        fill([t_int; flipud(t_int)], [p_int; zeros(size(p_int))], ...
             [1 0.78 0.78], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    end
    plot(ns_off(time(mk_off)), VDS(mk_off) .* I_D(mk_off), 'Color', c_pwr, 'LineWidth', 1.5);
    hold off;
    ylabel('Power (W)');  xlabel('Time (ns)');
    xlim(xlims);  grid on;
    title(sprintf('Turn-Off: V_{DS}·I_D\nE_{off} = %.2f µJ', pts.E_off));

    % ---- Tile 6: Turn-On Power + E_on fill ----
    nexttile(6);
    hold on;
    if ~isnan(pts.t_I_10up) && ~isnan(pts.t_VDS_10dn)
        mk_int = time >= pts.t_I_10up & time <= pts.t_VDS_10dn;
        t_int  = ns_on(time(mk_int));
        p_int  = VDS(mk_int) .* I_D(mk_int);
        fill([t_int; flipud(t_int)], [p_int; zeros(size(p_int))], ...
             [0.78 0.78 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    end
    plot(ns_on(time(mk_on)), VDS(mk_on) .* I_D(mk_on), 'Color', c_pwr, 'LineWidth', 1.5);
    hold off;
    ylabel('Power (W)');  xlabel('Time (ns)');
    xlim(xlims);  grid on;
    title(sprintf('Turn-On: V_{DS}·I_D\nE_{on} = %.2f µJ', pts.E_on));
end

