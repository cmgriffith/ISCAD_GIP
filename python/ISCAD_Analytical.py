import numpy as np
from scipy import constants as cons
import warnings
import matplotlib.pyplot as plt
from functions import *
warnings.filterwarnings('ignore')

# load parameters from .json
file = 'ISCAD_parameters.json'
params = loadjson(file)


### STATOR RESISTANCE
Rs_DC = params.rho_s * params.l_stack / (params.h_slot * params.w_slot[0])
zeta = ( (np.pi * cons.mu_0 * params.f) / params.rho_s )**0.5 * params.h_slot
print("stator DC resistance = ", ("{:.3f}".format(Rs_DC*1e3)), "mOhm")
K_R = zeta * (np.sinh(2*zeta) + np.sin(2*zeta)) / (np.cosh(2*zeta) - np.cos(2*zeta)) 
Rs_AC = Rs_DC * K_R
print("stator AC resistance = ", ("{:.3f}".format(Rs_AC*1e3)), "mOhm")

### SELF INDUCTANCE
Ls_s = np.pi/2 * cons.mu_0 * (params.r * params.l_stack / params.ag)
print("self inductance = " , ("{:.3f}".format(Ls_s*1e3)), "mH")

### LEAKAGE INDUCTANCE
# for DC current in rectangular slots
lambda_slot = (params.h_slot / 3*params.w_slot[0]) + (params.h_so / params.w_so)
Ls_sigma_DC = cons.mu_0 * params.l_stack * lambda_slot
print("DC leakage inductance = " ,  ("{:.3f}".format(Ls_sigma_DC*1e6)), "uH")
# AC corrected to fundamental f
K_L = 3 / (2*zeta) * (np.sinh(2*zeta) - np.sin(2*zeta)) / (np.cosh(2*zeta) - np.cos(2*zeta)) 
Ls_sigma_AC = cons.mu_0 * params.l_stack * lambda_slot * K_L
print("AC leakage inductance = " ,  ("{:.3f}".format(Ls_sigma_AC*1e6)), "uH")


### MAIN INDUCTANCE
sum = 0
mutuals = np.zeros(round(params.Qs/2/params.p + 1))
for k in range(1,round(params.Qs/2/params.p + 1)):
    linkage = np.cos( params.p * (k-1) * 2*np.pi / params.Qs)**2
    mutuals[k] = linkage
    sum = sum + linkage
Ls_m = Ls_s * sum
print("main inductance = " ,  ("{:.3f}".format(Ls_m*1e3)), "mH")

print("------")
print("stator slot count: ", params.Qs)
print("pole count: ", params.p)
print("number of slots summed for main L: ", round(params.Qs/2/params.p + 1))
print("analytical self inductance sum result: ", sum)
print("analytical linkage coeffs: ", linkage)
print("mutuals array size: ", np.shape(mutuals))



sumslots = np.arange(1,round(params.Qs/2/params.p + 1))
figure1 = plt.figure()
axes1 = figure1.add_subplot(1, 1, 1)
axes1.plot(sumslots, Ls_s*mutuals[1:])
plt.xlabel("Qs Slot Number")
plt.ylabel("Mutual Inductance")
plt.show()


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

C_DC_min = Apk / ( params.V_DC_rip_max * 2 * np.pi * params.f_SW )