import numpy as np

# parameters [SI units]
# ----------------------USER INPUT-----------------------
# FEMM_createmodel() MACHINE CONFIG
Qs = 60 # number of stator bars
Qs_active = 60 # active slots, must be integer divisor of Qs
p = 4 # number of pole pairs
Qr = 70 # number of rotor bars

# FEMM_createmodel() MACHINE DIMENSIONS
ag = 0.5
l_stack = 250
r_o = 115 # stator outer radius [mm]
r_so = 86.5 # radius at back of slot [mm]
r = 75 # stator bore radium [mm]
r_yoke = 45 # rotor yoke radius [mm]
# slot dimensions
w_slot = np.array([3 , 3.5]) # slot width: inner, outer
w_so = 1111 # slot opening width
h_slot = 10 # slot height
h_so = 1111 # slot opening height
w_bar = np.array([2 , 2.5]) # bar width: inner , outer
w_bo = 1111 # bar opening width
h_bar = 10 # bar height



# Material parameters
rho_s = 2.8e-08 # resistivity of slot aluminium


# Electrical parameters
Apk = 120 # circuit peak amps
f = 25000 # stator currents frequency

# Analysis parameters
k = 13 # max harmonic order for DFT
# --------------------------END--------------------------

# derived parameters
m = Qs/p # phase count, also slots per pole pair
r_r = r - ag # rotor radius [mm]
parameters = [Qs, Qr, p, l_stack, h_slot, w_slot, h_bar, w_bar, ag, 
 r_o, r_so, r, r_r, r_yoke]
phaseangles = np.arange(0,2*np.pi,Qs)