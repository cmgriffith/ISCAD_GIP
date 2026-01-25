import femm
import numpy as np
import subprocess
from functions import *
from plots import *

# load parameters from .json
file = 'ISCAD_parameters.json'
params = loadjson(file)

globals().update(params) # HORRIBLE CODE I AM SORRY, cba to rewrite


# Runtime Controls , Boolean
CreateModel = 1
Solve = 1
Analyze = 1
CycleRender = 0

# FILES
filename = "ISCAD_FEMM.FEM"
directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_FEMM.bak2\\"
path = directory + filename
bpr = 25 # bitmaps per rotation, if CycleRender == True


# derived parameters
params.m = Qs/p # phase count, also slots per pole pair
params.r_r = r - ag # rotor radius [mm]
params.phaseangles = np.arange(0,2*np.pi,Qs)

# parameters = [Qs, Qr, p, l_stack, h_slot, w_slot, h_bar, w_bar, ag, 
# r_o, r_so, r, r_r, r_yoke]



# create FEMM model
if CreateModel == True:
    FEMM_createmodel(path,params) # create, save, leave open in FEMM
else:
    print("Opening FEMM file: ", path)
    femm.openfemm()
    femm.opendocument(path)
    # femm.opendocument(path) # open saved
    femm.mi_probdef(0,'meters','planar',1e-008, l_stack)
'''
    try: 
        femm.opendocument(path)
    except: 
        print("File not found, ending script.")
        quit()
'''   


if Solve == True:
    Aph = FEMM_currents(params)
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
    B_ag = FEMM_contourplots(params) # get data from FEMM solution
    Bg_k = DFT(B_ag,k)

    Ls_s, Ls_m = FEMM_integrals(params)

    print("Numerical results...")
    print("self inductance = " ,  ("{:.3f}".format(Ls_s*1e3)), "mH")
    print("main inductance = " ,  ("{:.3f}".format(Ls_m*1e3)), "mH")

    print("Analytical results...")
    subprocess.run(["python3", "ISCAD_Analytical.py"])

    plot_B(B_ag)
    plot_DFT(Bg_k)
    Aph = FEMM_currents(params)
    plot_Aph(Aph)


'''
# BASH commands to render video and then GIF...
# TBC refactor the commands into a python subprocess()
ffmpeg -framerate 30 -i bitmap%d.bmp -vf "scale=1046:816" -c:v libx264 -pix_fmt yuv420p output.mp4
ffmpeg   -i output.mp4   -r 15   -vf scale=512:-1   -ss 00:00:00 -to 00:00:03   output.gif
'''
