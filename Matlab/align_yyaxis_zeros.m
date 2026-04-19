function align_yyaxis_zeros(ax)
% Sets ylim on both sides of a yyaxis axes so that y=0 falls at the same
% fractional position, without clipping any plotted data.
    yyaxis(ax, 'left');   yl = ylim(ax);
    yyaxis(ax, 'right');  yr = ylim(ax);
    yl = [min(yl(1), 0), max(yl(2), 0)];
    yr = [min(yr(1), 0), max(yr(2), 0)];
    fl = -yl(1) / (yl(2) - yl(1));
    fr = -yr(1) / (yr(2) - yr(1));
    f  = min(max(max(fl, fr), 0.02), 0.98);
    sl = max(-yl(1) / f, yl(2) / (1 - f));
    sr = max(-yr(1) / f, yr(2) / (1 - f));
    yyaxis(ax, 'left');   ylim(ax, [-f*sl, (1-f)*sl]);
    yyaxis(ax, 'right');  ylim(ax, [-f*sr, (1-f)*sr]);
end
