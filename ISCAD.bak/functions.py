import femm
import numpy as np

def FEMM_solve(step):
    # Save, analyze
    femm.mi_analyze(0)
    femm.mi_loadsolution()
    # View results
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag") 
    # prepare for screenshot
    femm.mo_resize(1050,800)
    femm.mo_zoomnatural()
    print("solved: ",step)

def FEMM_bitmap(step,bpr):
    filename = "bitmap_" + str(step)
    directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\FEMM\\"
    path = directory + filename
    femm.mo_savebitmap(path)
    print("saved bitmap to: ",path)

def FEMM_createmodel(parameters):

    Qs, Aph, circuits, l_stack, h_slot, w_slot, r_o, r_so, r, r_r, r_yoke = parameters
    # Initialises and defines model units, type, accuracy and length
    femm.openfemm()
    femm.newdocument(0)
    femm.mi_probdef(0,'millimeters','planar',1e-008, l_stack)

    # Function creates a FEMM model for a stator with 3 independent massive slots
    # non-periodic winding, whole stator modelled

    # Sets up materials and boundaries
    femm.mi_getmaterial('Air') # Gets region materials from library - non-magnetic airgap
    femm.mi_getmaterial('M-15 Steel') # Gets region materials from library - magnetic steel
    femm.mi_addmaterial('Coil') # Adds a region called Coil, defaults to non-magnetic region
    femm.mi_addboundprop('Zero',0,0,0,0,0,0,0,0,0) # Adds a zero potential boundary

    # [yoke, rotor, bore, outer]
    # nodes on x axis, positive and negative
    femm.mi_addnode(r_yoke,0)
    femm.mi_addnode(-r_yoke,0)
    femm.mi_addnode(r_r,0)
    femm.mi_addnode(-r_r,0)
    femm.mi_addnode(r,0)
    femm.mi_addnode(-r,0)
    femm.mi_addnode(r_o,0)
    femm.mi_addnode(-r_o,0)
    # circumference construction
    femm.mi_addarc(r_yoke,0,-r_yoke,0,180,1)
    femm.mi_addarc(-r_yoke,0,r_yoke,0,180,1)
    femm.mi_addarc(r_r,0,-r_r,0,180,1)
    femm.mi_addarc(-r_r,0,r_r,0,180,1)
    femm.mi_addarc(r,0,-r,0,180,1)
    femm.mi_addarc(-r,0,r,0,180,1)
    femm.mi_addarc(r_o,0,-r_o,0,180,1)
    femm.mi_addarc(-r_o,0,r_o,0,180,1)

    # set zero potential boundary at stator outer and yoke inner 
        # mi_setarcsegmentprop(maxsegdeg, ’propname’, hide, group)
    femm.mi_selectarcsegment(r_o,0)
    femm.mi_setarcsegmentprop(1,'Zero',0,0)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(-r_o,-1)
    femm.mi_setarcsegmentprop(1,'Zero',0,0)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(r_yoke,0)
    femm.mi_setarcsegmentprop(1,'Zero',0,0)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(-r_yoke,-1)
    femm.mi_setarcsegmentprop(1,'Zero',0,0)
    femm.mi_clearselected()

    # Assigns statorsteel region
    femm.mi_addblocklabel(0,(r_o+r_so)/2) # Adds stator region block label
    femm.mi_selectlabel(0,(r_o+r_so)/2) # Selects label
    femm.mi_setblockprop('M-15 Steel', 1, 0, 0, 0, 0, 0) # Sets region to magnetic steel
    femm.mi_clearselected()
    # Assigns airgap region
    femm.mi_addblocklabel(0, (r+r_r)/2) # Adds airgap region block label
    femm.mi_selectlabel(0, (r+r_r)/2) # Selects label
    femm.mi_setblockprop('Air', 1, 0, 0, 0, 0, 0) # sets selected regions to air
    femm.mi_clearselected()
    # Assigns rotor steel region
    femm.mi_addblocklabel(0,(r_r+r_yoke)/2) # Adds stator region block label
    femm.mi_selectlabel(0,(r_r+r_yoke)/2) # Selects label
    femm.mi_setblockprop('M-15 Steel', 1, 0, 0, 0, 0, 0) # Sets region to magnetic steel
    femm.mi_clearselected()
    # Assigns yoke zero boundary / no mesh
    femm.mi_addblocklabel((r_yoke)/2,0) # Adds airgap region block label
    femm.mi_selectlabel((r_yoke)/2,0) # Selects label
    femm.mi_setblockprop('<No Mesh>', 1, 0, 0, 0, 0, 0) # sets selected regions to air
    femm.mi_clearselected()
    

    if r_so != r + h_slot:
        print("slot dimension mismatch! = ", (r_so-(r+h_slot)))
    else:
        print ("slot dimension good")

    # draw slots

    # polar coordinates, positive y axis as angle reference
    slot_angle = np.arange(0,2*np.pi,2*np.pi/Qs) # Slot separation angles array in radians =2.094..
    Th_si = np.atan((w_slot/2)/r) # inner half pitch angle
    Th_so = np.atan((w_slot/2)/r_so) # outer half pitch angle
    # print("slot half pitches",Th_si,Th_so)

    # polar coordinates arrays for slot nodes
    R_slots = np.array([r,r_so,r_so,r])
    Th_slots = np.array([-Th_si,-Th_so,Th_so,Th_si])
    # R_center = np.array()

    # cartesian coordinates arrays for slot nodes
    x = np.multiply( R_slots , np.sin(Th_slots + slot_angle[...,None]) )
    y = np.multiply( R_slots , np.cos(Th_slots + slot_angle[...,None]) )    
    # row = slot(3), column = coord(4)
    print("x",x)
    print("y",y)

    # x_center = 
    # y_center


    for s in range(0,Qs): # for each slot
        for n in range(0,len(R_slots)): # for each coordinate in slot
            femm.mi_addnode(x[s,n],y[s,n]) # add node at coordinate
    for s in range(0,Qs): # for each slot
        for n in range(0,len(R_slots)-1):
            # add segments. stator bore forms inner arc segment
            femm.mi_addsegment( x[s,n] , y[s,n] , x[s,n+1] , y[s,n+1] )

    # for n in range (0,Qs):
    #     femm.mi_addcircprop(circuits[n], Aph[n], 1)

    # set coil block properties, 
    # start at 12 o'clock slot, slots into groups 1,2,3 clockwise
    for s in range(0,Qs):
    #  for n in range(0,len(R_slots))
        x_centre = np.mean(x[s,:])
        y_centre = np.mean(y[s,:])
        # add block labels, assign to circuits
        femm.mi_addblocklabel(x_centre,y_centre)
        femm.mi_selectlabel(x_centre,y_centre)
        femm.mi_setblockprop('Coil', 1, 0, circuits[s], 0, (s+1), 0.5)
        femm.mi_clearselected()

    femm.mi_zoomnatural() # Displays the model outline to fit window

    # Saves model - note the directory file structure with '\\' replacing the normal '\'
    filename = "ISCAD_FEMM.FEM"
    directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\FEMM\\"
    path = directory + filename
    femm.mi_saveas(path)
