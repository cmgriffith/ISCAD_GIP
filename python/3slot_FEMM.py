import femm
import numpy as np
import json
from functions import *

# FILES
filename = "3slot_FEMM.FEM"
directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_GIP\\python\\FEMM output\\"
path = (directory, filename)

# load parameters from .json
file = '3slot_parameters.json'
params = loadjson(file)
globals().update(params) # update params in current python file 

femm.openfemm()
femm.newdocument(0)
femm.mi_probdef(0,'meters','planar',1e-008, l_stack,10)

# Sets up materials and boundaries
femm.mi_getmaterial('Air') # Gets region materials from library - non-magnetic airgap
femm.mi_getmaterial('M-15 Steel') # Gets region materials from library - magnetic steel
femm.mi_addmaterial('Coil') # Adds a region called Coil, defaults to non-magnetic region
femm.mi_addboundprop('Zero',0,0,0,0,0,0,0,0,0) # Adds a zero potential boundary

# draw active slots
for slot in range(0,Qs):
    x = slot*spacing_slot
    femm.mi_drawline(0+x, 0, 0+x, h_slot)
    femm.mi_drawline(0+x, h_slot, w_slot+x,h_slot)
    femm.mi_drawline(w_slot+x,h_slot, w_slot+x, 0)
    femm.mi_drawline(w_slot+x, 0, 0+x, 0)

# draw back iron
femm.mi_drawline(-iron_back, -iron_back, -iron_back, h_slot + iron_top)
femm.mi_drawline(-iron_back, h_slot + iron_top, ((Qs-1)*spacing_slot + Qs*w_slot + iron_back), (h_slot + iron_top))
femm.mi_drawline((Qs-1)*spacing_slot + Qs*w_slot + iron_back, (h_slot + iron_top), ((Qs-1)*spacing_slot + Qs*w_slot + iron_back), -iron_back)
femm.mi_drawline((Qs-1)*spacing_slot + Qs*w_slot + iron_back, -iron_back, -iron_back, -iron_back)

femm.mi_zoomnatural()