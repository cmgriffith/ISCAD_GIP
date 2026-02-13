import femm
import numpy as np
import json
from functions import *
import time
start = time.time()

# Runtime Controls , Boolean
CreateModel = 1
Solve = 0
Analyse = 1
ShowFEMM = 0  # runs faster if 0

# FILES
filename = "3slot_FEMM.FEM"
# directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_GIP\\python\\FEMM output\\"
directory = "C:\\Users\\Charlie\\Desktop\\ISCAD_GIP\\python\\"
path = (directory, filename)

# load parameters from .json
file = '3slot_parameters.json'
params = loadjson(file)
globals().update(params) # update params in current python file 
derived_params(params, file) # add derived values to params
globals().update(params) # update params in current python file 

print("------")
femm.openfemm(not ShowFEMM)
femm.newdocument(0)
femm.mi_probdef(f,'meters','planar',1e-008, l_stack,30)


### MECHANICAL -----------------------------------------------------------------------------------------------

def FEMM_create3slot(path):

    # Sets up materials and boundaries
    femm.mi_getmaterial('Air') # Gets region materials from library - non-magnetic airgap
    femm.mi_getmaterial('M-15 Steel') # Gets region materials from library - magnetic steel
    femm.mi_addmaterial('Coil') # Adds a region called Coil, defaults to non-magnetic region
    femm.mi_addboundprop('Zero',0,0,0,0,0,0,0,0,0) # Adds a zero potential boundary
    femm.mi_addboundprop('Neumann',0,0,0,0,0,0,0,0,2) # Adds a zero potential boundary
    
    # set up IABC
    RHS_bound = iron_back + (Qs*gap_slot) + (Qs*w_slot) - gap_mut # x co-ord of end of slots (RHS iron?)
    LHS_bound = -iron_back - gap_slot + gap_mut
    femm.mi_makeABC(7, RHS_bound-LHS_bound, gap_slot+w_slot/2, h_slot/2, 0)
    airzone = np.array([LHS_bound - iron_back , h_slot + iron_top])
    femm.mi_addblocklabel(*airzone) # back iron set material
    femm.mi_selectlabel(*airzone)
    femm.mi_setblockprop('Air', 1, Qs+Qm+4)
    femm.mi_clearselected()    

    d = np.array([1e-4,1e-4]) # 0.1mm offset variable

    # draw back iron
    coords = np.array([(LHS_bound, -iron_back), (LHS_bound, h_slot), (RHS_bound, h_slot), (RHS_bound, -iron_back)])
    femm.mi_drawline(*coords[0],*coords[1])
    femm.mi_selectsegment(*coords[0]+d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()
    femm.mi_drawline(*coords[1],*coords[2]) # Neumann
    femm.mi_selectsegment(*coords[1]+d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+2)
    femm.mi_clearselected()
    femm.mi_drawline(*coords[2],*coords[3])
    femm.mi_selectsegment(*coords[3]+d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()
    femm.mi_drawline(*coords[3],*coords[0])
    femm.mi_selectsegment(*coords[3]-d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()
    femm.mi_addblocklabel(*coords[0]+d) # back iron set material
    femm.mi_selectlabel(*coords[0]+d)
    femm.mi_setblockprop('M-15 Steel', 1, Qs+Qm+3)
    femm.mi_clearselected()


    # draw active slots, along x from origin (0,0)
    coords = np.array([(0, 0), (0, h_slot), ( w_slot,h_slot), (w_slot,0)])
    for slot in range(0,Qs):
        shift = np.array([slot*(gap_slot+w_slot),0]) # offset
        femm.mi_drawline(*coords[0]+shift, *coords[1]+shift)
        femm.mi_drawline(*coords[1]+shift, *coords[2]+shift)
        femm.mi_drawline(*coords[2]+shift, *coords[3]+shift)
        femm.mi_drawline(*coords[3]+shift, *coords[0]+shift)
        # assign and label
        femm.mi_addblocklabel(np.mean(coords[0]+shift+coords[3]+shift), 0.5*h_slot)
        femm.mi_selectlabel(np.mean(coords[0]+shift+coords[3]+shift), 0.5*h_slot)
        femm.mi_addcircprop(str(slot),0,0)
        femm.mi_setblockprop('Coil',1,0, str(slot), 0, str(slot), 0.5)
        femm.mi_clearselected()

    # draw mutual slots
    offset = (h_slot - h_mut) / 2
    coords = np.array([(0, 0+offset), (0, h_slot-offset), ( w_mut,h_slot-offset), (w_mut,0+offset)]) # mutual slot coords from origin
    coords = coords[:] + (-gap_slot,0) + (gap_mut,0) # coord correction to start drawing from one space left of origin
    # print("gap_slot: ",gap_slot)
    # print("mut slot coords: ",coords)
    for gap in range(0,Qs+1):
        slotshift = np.array([gap*(w_slot+gap_slot),0]) # transloaction for drawing in each slot's gap
        for mut in range(0,Qm):
            mutshift = np.array([mut*(gap_mut+w_mut),0]) # translocation for drawing each mutual across a slot gap
            femm.mi_drawline(*(coords[0]+mutshift+slotshift), *(coords[1]+mutshift+slotshift))
            femm.mi_drawline(*(coords[1]+mutshift+slotshift), *(coords[2]+mutshift+slotshift))
            femm.mi_drawline(*(coords[2]+mutshift+slotshift), *(coords[3]+mutshift+slotshift))
            femm.mi_drawline(*(coords[3]+mutshift+slotshift), *(coords[0]+mutshift+slotshift))
            # femm.mi_addblocklabel(mutshift[0]+slotshift[0]+np.average(coords[0]+coords[3]), np.average(coords[0,1]+coords[2,1]))
            # femm.mi_selectlabel(mutshift[0]+slotshift[0]+np.average(coords[0]+coords[3]), np.average(coords[0,1]+coords[2,1]))
            femm.mi_addblocklabel(*(mutshift+slotshift+np.mean([coords[0],coords[2]], axis=0)))
            femm.mi_selectlabel(*(mutshift+slotshift+np.mean([coords[0],coords[2]], axis=0)))
            femm.mi_addcircprop(str(Qs+gap*Qm+mut),0,0)
            femm.mi_setblockprop('Coil',1,0, str(Qs+gap*Qm+mut),0, str(Qs+gap*Qm+mut), 0.5)
            femm.mi_clearselected()
            # print("gap:",gap," mut:",mut," total:",(Qs+gap*Qm+mut))
            # print("Qm: ", (Qs+gap*Qm+mut), "slotshift: ",slotshift, " gapshift: ", mutshift)


    # draw top iron
    coords = np.array([(LHS_bound, h_slot), (LHS_bound, (h_slot + iron_top)), (RHS_bound, (h_slot + iron_top)), (RHS_bound, h_slot)])

    femm.mi_drawline(*coords[0],*coords[1])
    femm.mi_selectsegment(*coords[1]-d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()
    femm.mi_drawline(*coords[1],*coords[2])
    femm.mi_selectsegment(*coords[1]+d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()
    femm.mi_drawline(*coords[2],*coords[3])
    femm.mi_selectsegment(*coords[3]+d)
    femm.mi_setsegmentprop('Neumann',1,1,0,Qs+Qm+1)
    femm.mi_clearselected()

    femm.mi_addblocklabel(*coords[0]+d) # back iron set material
    femm.mi_selectlabel(*coords[0]+d)
    femm.mi_setblockprop('M-15 Steel', 1, Qs+Qm+3)
    femm.mi_clearselected()

    ### ELECTRICAL -------------------------------------------------------------------------------------------
    Aph = np.zeros(Qs)
    phaseangle = np.arange(0,2*np.pi,1/Qs*2*np.pi)
    for slot in range(0,Qs):
        Aph[slot] = Apk*np.sin(phaseangle[slot] + 0.5*np.pi )
        femm.mi_setcurrent(str(slot),Aph[slot])

    # print("slot currents: ", Aph)
    # print(type(Aph[0]))


    femm.mi_saveas("".join(path))
    end = time.time()
    print("Model Created in ", ("{:.2f}".format((end-start))), "s")
    print("Model saved to: ","".join(path))
    femm.mi_zoomnatural()
    print("------")


### SOLVE -------------------------------------------------------------------------------------------------

if CreateModel == 1:
    print("Creating FEMM file...")
    FEMM_create3slot(path)
else:
    print("Opening FEMM file: ", path)
    femm.openfemm()
    femm.opendocument("".join(path))

if Solve == 1:
    print("Solving FEMM model...")
    start = time.time()
    femm.mi_analyze(0)
    femm.mi_loadsolution()
    femm.mo_shownames()
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag")
    femm.mo_zoomnatural()
    end = time.time()
    print("Model Solved in ", ("{:.2f}".format((end-start))), "s")
    print("------")


if Analyse == 1:
    print("Analysing FEMM model...")
    start = time.time()
    Aph = np.zeros(Qs)
    active_slot = 0
    Aph[active_slot] = Apk
    # Aph = [Apk,0,0]
    for slot in range(0,Qs):
        femm.mi_setcurrent(str(slot),Aph[slot])
    
    femm.mi_analyze(0)
    femm.mi_loadsolution()
    femm.mo_zoomnatural()
    femm.mo_showdensityplot(1,0,3,0,"jmag") 

    # SELF INDUCTANCE
    femm.mo_seteditmode('area')
    femm.mo_groupselectblock(active_slot)
    AJ = femm.mo_blockintegral(0)
    femm.mo_clearblock()
    current,*_ = femm.mo_getcircuitproperties(str(active_slot))
    Ls_s = AJ / (current**2)

    # MUTUAL INDUCTANCE
    mutuals = [];
    for slot in range(0,Qs+Qm_total):
        if slot == active_slot:
            pass
        else:
            femm.mo_seteditmode('area')
            femm.mo_groupselectblock(str(slot))
            A = femm.mo_blockintegral(1) # A integral
            a = femm.mo_blockintegral(5) # cross section area
            femm.mo_clearblock()
            mutuals.append((n/(np.float64(Apk)*a)) * A)
    Ls_m = Ls_s + np.sum(np.absolute(mutuals))

    end = time.time()
    print("Model Analysed in ", ("{:.2f}".format((end-start))), "s")
    print("------")

    print("Numerical results...")
    print("self inductance = " ,  ("{:.3f}".format(Ls_s*1e3)), "mH")
    print("main inductance = " ,  ("{:.3f}".format(Ls_m*1e3)), "mH")
    print("------")


(input("Press ENTER to close FEMM"))
