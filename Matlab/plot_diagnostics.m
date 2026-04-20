function plot_diagnostics(time, VGS, I_D, VDS, T1, T2, VGS_ss, pts, label)
% 3×2 switching diagnostic figure.
% label (optional): string used in figure/tile title instead of default "T1 = XX µs".
    if nargin < 9 || isempty(label)
        label = sprintf('T_1 = %d µs', T1);
    end

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

    figure('Name', sprintf('Diagnostic — %s', label), ...
           'NumberTitle', 'off', 'Position', [40 40 1200 900]);
    tl = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Switching Metric Diagnostic  —  %s', label), ...
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
