import femm
import numpy as np
import os ###
import subprocess
from functions import *
from plots import *

# Runtime Controls , Boolean
CreateModel = 1
Solve = 0 # solve full machine field
Analyze = 0 # solve for ag flux, integrals
Plots = 0 # show analysis plots
CycleRender = 0 # render an animation

# FILES
filename = "ISCAD_FEMM.FEM"
directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_GIP\\python\\FEMM output\\"
# path = (directory + filename)
path = (directory, filename)
bpr = 25 # bitmaps per rotation, if CycleRender == True

# load parameters from .json
file = 'ISCAD_parameters.json'
params = loadjson(file)
def derived_params(params):
    globals().update(params) # HORRIBLE CODE I AM SORRY, cba to rewrite
    params.m = Qs/p # phase count, also slots per pole pair
    params.r_r = r - ag # rotor radius [mm]
    params.phaseangles = np.arange(0,2*np.pi,Qs)
    globals().update(params) # HORRIBLE CODE I AM SORRY, cba to rewrite
derived_params(params) # add derived values to params
globals().update(params) # update params in current python file 

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

Aph = FEMM_currents(params)
plot_Aph(Aph)

'''
# compute inductances
Ls_s, Ls_m, mutuals = FEMM_integrals(path,params)

# print results
plot_mutuals(Qs, mutuals)
print("NUMERICAL RESULTS...")
print("self inductance = " ,  ("{:.3f}".format(Ls_s*1e3)), "mH")
print("main inductance = " ,  ("{:.3f}".format(Ls_m*1e3)), "mH")
print("ANALYTICAL RESULTS...")
subprocess.run(["python3", "python/ISCAD_Analytical.py"])
'''

if Solve == True:
    Aph = FEMM_currents(params)
    FEMM_solve(0)
elif CycleRender == True:
    theta = np.arange(0,2*np.pi,1/bpr*2*np.pi)
    for step in range(0,bpr):
        FEMM_currents(Apk, theta[step])
        FEMM_solve(step)
        FEMM_bitmap(path,step)
else:
    print("No full machine solve was requested.")
    # quit()

if Analyze == True:
    if Solve != True: # ensure model has been solved
        # Aph = FEMM_currents(params)
        # FEMM_solve(path,0)
        print("")
    print("Analysing...")
    femm.mi_loadsolution()
    # View results and prepare for screenshot
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag") 
    femm.mo_resize(1050,800)
    femm.mo_zoomnatural()

    # get airgap flux from FEMM solution
    B_ag = FEMM_contourplots(params)
    Bg_k = DFT(B_ag,k)
    print("Airgap flux extracted")

    # compute inductances
    Ls_s, Ls_m, mutuals = FEMM_integrals(path,params)

    # print results
    print("Numerical results...")
    print("self inductance = " ,  ("{:.3f}".format(Ls_s*1e3)), "mH")
    print("main inductance = " ,  ("{:.3f}".format(Ls_m*1e3)), "mH")
    print("Analytical results...")
    subprocess.run(["python3", "ISCAD_Analytical.py"])
    # show plots
    if Plots == 1:
        #plot_B(B_ag)
        #plot_DFT(Bg_k)
        Aph = FEMM_currents(params)
        plot_Aph(Aph)
        plot_mutuals(mutuals)


'''
# BASH commands to render video and then GIF...
# TBC refactor the commands into a python subprocess()
ffmpeg -framerate 30 -i bitmap%d.bmp -vf "scale=1046:816" -c:v libx264 -pix_fmt yuv420p output.mp4
ffmpeg   -i output.mp4   -r 15   -vf scale=512:-1   -ss 00:00:00 -to 00:00:03   output.gif
'''
