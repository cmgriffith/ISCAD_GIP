import numpy as np
from scipy import constants as cons
import matplotlib.pyplot as plt

# parameters [SI]
# ----------------------USER INPUT-----------------------
Qs = 60 # or 60?
p = 2
r_bore = 75
l_stack = 250
ag = 0.5
b_slot = 4 # slot width
b_so = 16 # slot opening width
h_slot = 1 # slot height
h_so = 1 # slot opening height
# -------------------------------------------------------
m = Qs/p # phase count
mu_0 = cons.mu_0 # constants

### SELF INDUCTANCE

Ls_s = np.pi/2*mu_0 * (r_bore * l_stack / ag)
print("self inductance = " , round(Ls_s,3), "H")

### LEAKAGE INDUCTANCE

lambda_slot = (h_slot / 3*b_slot) + (h_so / b_so)
Ls_sigma_DC = mu_0 * l_stack * lambda_slot
print("leakage inductance DC = " , round(Ls_sigma_DC,3), "H")

# Ls_sigma_AC = mu_o * l_stack * lambda_slot * K_L


### MAIN INDUCTANCE

Ls_m_arr = np.zeros(Qs)
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

print("main inductance = " , round(Ls_m,3), "H")


'''
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