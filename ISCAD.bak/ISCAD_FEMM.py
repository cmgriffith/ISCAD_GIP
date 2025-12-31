import femm
import numpy as np
from functions import *
# from functions import FEMM_createmodel, FEMM_solve, FEMM_bitmap

# parameters [SI]
# ----------------------USER INPUT-----------------------
# machine
Qs = 3 # or 60?
p = 2 # number of pole pairs
m = Qs/p # phase count, also slots per pole pair
# electrical
Apk = 200 # circuit peak amps
Aph = Apk*np.array([1,-0.5,-0.5])
circuits = np.array(['R','Y','B'])
# dimension
ag = 0.5
l_stack = 250
w_slot = 3.5 # slot width
h_slot = 10 # slot height
r_o = 115 # stator outer radius [mm]
r_so = 86.5 # radius at back of slot [mm]
r = 75 # stator bore radium [mm]
r_r = r - ag # rotor radius [mm]
r_yoke = 45 # rotor yoke radius [mm]
# -------------------------------------------------------

# create FEMM model
parameters = [Qs, Aph, circuits, l_stack, h_slot, w_slot, r_o, r_so, r, r_r, r_yoke]
FEMM_createmodel(parameters)


# solve across a rotation
bpr = 36 # bitmaps per rotation
Aph = np.zeros(Qs)
for step in range(0,bpr): # for each step
    theta = 2*np.pi/bpr*step # calc. electrical angle
    for slot in range(0,Qs): #  assign current
        Aph[slot] = Apk*np.cos( theta - slot*(2*np.pi/Qs) )
        femm.mi_addcircprop(circuits[slot], Aph[slot], 1)
    FEMM_solve(step)
    FEMM_bitmap(step,bpr)


'''
# BASH commands to render video and then GIF...
# TBC refactor the commands into a python subprocess()
ffmpeg -framerate 30 -i bitmap%d.bmp -vf "scale=1046:816" -c:v libx264 -pix_fmt yuv420p output.mp4
ffmpeg   -i output.mp4   -r 15   -vf scale=512:-1   -ss 00:00:00 -to 00:00:03   output.gif
'''
