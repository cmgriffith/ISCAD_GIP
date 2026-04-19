L = 50e-6;
T1 = DPT(197, 48, L);
function [T1, T2, T3] = DPT(I_des, V_DC, L)
T1 = L*I_des/V_DC;
if T1 > 200e-6
    fprintf('T1 exceeds 200us, consider using a smaller load inductor\n');
end
Rds_on = 1.05e-3/3;
E = Rds_on*(0.5*T1*I_des^2);

end
fprintf('T1: %0.0f us\n', T1*10^6)