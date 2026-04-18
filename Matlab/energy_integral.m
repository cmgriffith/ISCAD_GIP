function E = energy_integral(time, VDS, I_D, t_start, t_end)
% Integrates VDS * I_D over [t_start, t_end] using the trapezoidal rule.
% Returns NaN if either bound is NaN or t_start >= t_end.
    if any(isnan([t_start, t_end])) || t_start >= t_end
        E = NaN;  return;
    end
    mask = time >= t_start & time <= t_end;
    E    = trapz(time(mask), VDS(mask) .* I_D(mask));
end
