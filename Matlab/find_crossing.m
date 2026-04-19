function t_cross = find_crossing(time, signal, threshold, direction, t_start, t_end)
% Returns the time of the first threshold crossing within [t_start, t_end].
% direction: +1 = rising, -1 = falling.  Returns NaN with a warning if not found.
    mask = time >= t_start & time <= t_end;
    t_w  = time(mask);
    s_w  = signal(mask);
    d    = diff(sign(s_w - threshold));
    if direction > 0
        idx = find(d > 0, 1);
    else
        idx = find(d < 0, 1);
    end
    if isempty(idx)
        t_cross = NaN;
        warning('find_crossing: not found  thresh=%.3f  dir=%+d  window=[%.3f %.3f] µs', ...
                threshold, direction, t_start*1e6, t_end*1e6);
    else
        t_cross = t_w(idx + 1);
    end
end
