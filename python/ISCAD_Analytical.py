import numpy as np
from scipy import constants as cons
# import matplotlib.pyplot as plt
from ISCAD_parameters import *

import warnings
warnings.filterwarnings('ignore')

### STATOR RESISTANCE
Rs_DC = rho_s * l_stack / (h_slot * w_slot[0])
zeta = ( (np.pi * cons.mu_0 * f) / rho_s )**0.5 * h_slot
print("stator DC resistance = ", ("{:.2e}".format(Rs_DC)), "Ohm")
K_R = zeta * (np.sinh(2*zeta) + np.sin(2*zeta)) / (np.cosh(2*zeta) - np.cos(2*zeta)) 
Rs_AC = Rs_DC * K_R
print("stator AC resistance = ", round(Rs_AC,1), "Ohm")

### SELF INDUCTANCE

Ls_s = np.pi/2 * cons.mu_0 * (r * l_stack / ag)
print("self inductance = " , round(Ls_s,3), "H")

### LEAKAGE INDUCTANCE
# for DC current in rectangular slots
lambda_slot = (h_slot / 3*w_slot[0]) + (h_so / w_so)
Ls_sigma_DC = cons.mu_0 * l_stack * lambda_slot

K_L = 3 / (2*zeta) * (np.sinh(2*zeta) - np.sin(2*zeta)) / (np.cosh(2*zeta) - np.cos(2*zeta)) 
Ls_sigma_AC = cons.mu_0 * l_stack * lambda_slot * K_L

print("leakage inductance DC = " , round(Ls_sigma_DC,3), "H")


### MAIN INDUCTANCE
sum = 0
for k in range(1,round(Qs/2/p + 1)):
    sum = sum + np.cos( p * (k-1) * 2*np.pi / Qs)**2
Ls_m = Ls_s + sum

print("main inductance = " , round(Ls_m,3), "H")


'''
# main inductances for range of Qs
for slot in range(1,Qs+1): 
    Qs_loop = slot
    m_loop = Qs_loop/p
    sum,Ls_m = 0,0
    for k in range(1,round(m_loop/2)):
        x = np.cos( p * (k-1) * (2*np.pi / Qs_loop) ) ** 2
        sum = sum + x
    Ls_m = Ls_s * sum
    Ls_m_arr[slot-1] = Ls_m
    # print("slotcount= ", Qs_loop)
    # print(Ls_m)

# Plot bar chart
plt.bar(range(0,Qs), Ls_m_arr)
plt.xlabel("Qs (slot count)")
plt.ylabel("L_m [H] main inductance")
plt.title("Main inductance vs Q_s")
plt.axhline(
    y=Ls_s,
    color='red',
    linestyle='--',
    label="L_s self inductance"
)
plt.legend()
plt.show()
'''


### POWER ELECTRONICS

# parameters [SI]
# ----------------------USER INPUT-----------------------
V_DC = 48
f_SW = 25000
V_DC_rip_max = 2.6
Apk = 140 # [30s,330A] [5s,500A] 
# -------------------------------------------------------

C_DC_min = Apk / ( V_DC_rip_max * 2 * np.pi * f_SW )