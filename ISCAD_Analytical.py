import numpy as np
from scipy import constants as cons
import json
import sys
# import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

# Load parameters from JSON or Python file
use_json = '--json' in sys.argv or len(sys.argv) > 1 and sys.argv[1].endswith('.json')

if use_json:
    # Load from JSON file
    json_file = 'ISCAD_parameters.json'
    if len(sys.argv) > 1 and sys.argv[1].endswith('.json'):
        json_file = sys.argv[1]
    
    with open(json_file, 'r') as f:
        params = json.load(f)
    
    # Extract parameters
    Qs = params['machine_config']['Qs']
    Qs_active = params['machine_config']['Qs_active']
    p = params['machine_config']['p']
    Qr = params['machine_config']['Qr']
    
    ag = params['machine_dimensions']['ag']
    l_stack = params['machine_dimensions']['l_stack']
    r_o = params['machine_dimensions']['r_o']
    r_so = params['machine_dimensions']['r_so']
    r = params['machine_dimensions']['r']
    r_yoke = params['machine_dimensions']['r_yoke']
    w_slot = np.array(params['machine_dimensions']['w_slot'])
    w_so = params['machine_dimensions']['w_so']
    h_slot = params['machine_dimensions']['h_slot']
    h_so = params['machine_dimensions']['h_so']
    w_bar = np.array(params['machine_dimensions']['w_bar'])
    w_bo = params['machine_dimensions']['w_bo']
    h_bar = params['machine_dimensions']['h_bar']
    
    rho_s = params['material_parameters']['rho_s']
    
    Apk = params['electrical_parameters']['Apk']
    f = params['electrical_parameters']['f']
    
    k = params['analysis_parameters']['k']
    
    # Power electronics parameters (if available)
    if 'power_electronics' in params:
        V_DC = params['power_electronics']['V_DC']
        f_SW = params['power_electronics']['f_SW']
        V_DC_rip_max = params['power_electronics']['V_DC_rip_max']
        Apk_PE = params['power_electronics']['Apk']
    
    print("Parameters loaded from JSON file")
else:
    # Load from Python file
    from ISCAD_parameters import *
    print("Parameters loaded from Python file")

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
lambda_slot = (h_slot / (3*w_slot[0])) + (h_so / w_so)
Ls_sigma_DC = cons.mu_0 * l_stack * lambda_slot

K_L = (3 / (2*zeta)) * (np.sinh(2*zeta) - np.sin(2*zeta)) / (np.cosh(2*zeta) - np.cos(2*zeta)) 
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
if not use_json or 'power_electronics' not in params:
    V_DC = 48
    f_SW = 25000
    V_DC_rip_max = 2.6
    Apk_PE = 140 # [30s,330A] [5s,500A]
else:
    Apk_PE = params['power_electronics']['Apk']
# -------------------------------------------------------

C_DC_min = Apk_PE / ( V_DC_rip_max * 2 * np.pi * f_SW )