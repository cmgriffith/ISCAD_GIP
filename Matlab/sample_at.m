function val = sample_at(time, signal, t_ref)
% Returns the signal value at the sample nearest to t_ref.
    [~, idx] = min(abs(time - t_ref));
    val = signal(idx);
end
