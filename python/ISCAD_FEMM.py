import femm
import numpy as np
from functions import *
from plots import *

# Runtime Controls , Boolean
CreateModel = 0
Solve = 0
Analyze = 0
CycleRender = 0

# FILES
filename = "ISCAD_FEMM.FEM"
directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\FEMM\\"
path = directory + filename
bpr = 25 # bitmaps per rotation, if CycleRender == True

# parameters [SI units]
# ----------------------USER INPUT-----------------------

# FEMM_createmodel() MACHINE CONFIG
Qs = 60 # number of stator bars
p = 4 # number of pole pairs
Qr = 70 # number of rotor bars

# FEMM_createmodel() MACHINE DIMENSIONS
ag = 0.5
l_stack = 250
w_slot = np.array([3 , 3.5]) # slot width: inner, outer
h_slot = 10 # slot height
w_bar = np.array([2 , 2.5]) # bar width: inner , outer
h_bar = 10 # bar height
r_o = 115 # stator outer radius [mm]
r_so = 86.5 # radius at back of slot [mm]
r = 75 # stator bore radium [mm]
r_yoke = 45 # rotor yoke radius [mm]

# FEMM_currents() CIRCUITS
Qs_active = 60 # active slots, must be integer divisor of Qs
Apk = 120 # circuit peak amps

# Analysis parameters
k = 13 # max harmonic order for DFT
# --------------------------END--------------------------

# derived parameters
m = Qs/p # phase count, also slots per pole pair
r_r = r - ag # rotor radius [mm]
parameters = [Qs, Qr, p, l_stack, h_slot, w_slot, h_bar, w_bar, ag, 
 r_o, r_so, r, r_r, r_yoke]
phaseangles = np.arange(0,2*np.pi,Qs)




# create FEMM model
if CreateModel == True:
    FEMM_createmodel(parameters) # create, save, leave open in FEMM
else:
    print("Opening FEMM file: ", path)
    femm.openfemm()
    femm.opendocument(path)
    # femm.opendocument(path) # open saved
    femm.mi_probdef(0,'millimeters','planar',1e-008, l_stack)
'''
    try: 
        femm.opendocument(path)
    except: 
        print("File not found, ending script.")
        quit()
'''   


if Solve == True:
    Aph = FEMM_currents(Qs,Qs_active,p,Apk)
    FEMM_solve(0)
elif CycleRender == True:
    theta = np.arange(0,2*np.pi,1/bpr*2*np.pi)
    for step in range(0,bpr):
        FEMM_currents(Apk, theta[step])
        FEMM_solve(step)
        FEMM_bitmap(step)
else:
    print("No solve requested.")
    # quit()


if Analyze == True:
    print("Loading solution.")
    femm.mi_loadsolution()
    # View results
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag") 
    # prepare for screenshot
    femm.mo_resize(1050,800)
    femm.mo_zoomnatural()
    B_ag = FEMM_data(parameters) # get data from FEMM solution
    Bg_k = DFT(B_ag,k)
    plot_B(B_ag)
    plot_DFT(Bg_k)
    
    Aph = FEMM_currents(Qs,Qs_active,p,Apk)
    plot_Aph(Aph)


'''
# BASH commands to render video and then GIF...
# TBC refactor the commands into a python subprocess()
ffmpeg -framerate 30 -i bitmap%d.bmp -vf "scale=1046:816" -c:v libx264 -pix_fmt yuv420p output.mp4
ffmpeg   -i output.mp4   -r 15   -vf scale=512:-1   -ss 00:00:00 -to 00:00:03   output.gif
'''
