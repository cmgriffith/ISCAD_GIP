import femm
import numpy as np
import json


class DotDict(dict):
    """Dictionary subclass that allows dot notation access to nested keys"""
# to flatten hierarchy, include
# START
    def __init__(self, data=None):
        super().__init__()
        if data:
            flat = self._flatten_dict(data)
            for k, v in flat.items():
                self[k] = v

    def _flatten_dict(self, d):
        """Recursively flatten a nested dictionary"""
        items = {}
        for key, value in d.items():
            if isinstance(value, dict):
                items.update(self._flatten_dict(value))
            elif isinstance(value, list): ###
                items[key] = np.array(value) ###
            else:
                items[key] = value
        return items
# END        
    def __getattr__(self, key):
        try:
            value = self[key]
            # Recursively convert nested dicts to DotDict
            if isinstance(value, dict):
                return DotDict(value)
            if isinstance(value, list):    ###      # safety net
                value = np.array(value)###
                self[key] = value          ###      # cache conversion
            return value
        except KeyError:
            raise AttributeError(f"'DotDict' object has no attribute '{key}'")
    
    def __setattr__(self, key, value):
        if isinstance(value, list): ###
            value = np.array(value) ###
        self[key] = value
    
    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'DotDict' object has no attribute '{key}'")


def loadjson(file):
    # Read the JSON file
    with open(file, 'r') as f:
        data = json.load(f)
    # Convert to DotDict for dot notation access
    data = DotDict(data)
    # Now you can access with dot notation
    # print(data.machine.ag)  # Output: 0.0005
    # print(type(data.machine.ag))  # Output: <class 'float'>
    return data


def FEMM_solve(step):
    # Save, analyze
    print("solving model...")
    femm.mi_analyze(0)
    femm.mi_loadsolution()
    # View results
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag") 
    # prepare for screenshot
    femm.mo_resize(1050,800)
    femm.mo_zoomnatural()
    print("solved: ",step)


def FEMM_bitmap(step):
    filename = "bitmap_" + str(step)
    directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\FEMM\\"
    path = directory + filename
    femm.mo_savebitmap(path)
    print("saved bitmap to: ",path)


def FEMM_currents(params,theta=0):
    # create Qs evenly spaced angles within 2pi range
    print("Assigning circuit currents...")
    phaseangle = np.arange(0,2*np.pi,1/Qs*2*np.pi) # TBC redundant calls...
    Aph = np.zeros(Qs) # TBC redundant calls...
    factor = Qs / Qs_active
    is_divisible = lambda dividend, divisor : dividend % divisor == 0
    pp = p / 2 # pole pairs
    for slot in range(0,Qs):
        if is_divisible(slot+1 , factor) == True:
            # print(slot, " is active slot")
            print(theta, pp)
            Aph[slot] = Apk*np.sin( (phaseangle[slot] - theta) * pp ) # round for fp precision fix
            femm.mi_setcurrent(str(slot),Aph[slot])
        else:
            # print(slot, " is inactive slot")
            Aph[slot] = 0
            femm.mi_setcurrent(str(slot),Aph[slot])
        # print(Aph)
    return Aph


def FEMM_contourplots(params): # get data from FEMM solution

    # draw contour inside airgap, for MMF
    r_ag = r - ag/2
    femm.mo_seteditmode('contour')
    femm.mo_addcontour(0, r_ag)
    femm.mo_addcontour(0, -r_ag)
    femm.mo_bendcontour(180,1)
    femm.mo_addcontour(0, -r_ag)
    femm.mo_addcontour(0, r_ag)
    femm.mo_bendcontour(180,1)
    # save the data
    femm.mo_makeplot(2,361*4,'B_ag.txt',1) # Writes normal flux density (Bn) to text file
    data = np.loadtxt("B_ag.txt", dtype=float) # load data into python variable
    B_ag = data[:,1]
    return B_ag
    # print(B_ag)

def FEMM_integrals(params): # get data from FEMM solution
    ### block integrals for inductances
    # self inductance
    print("SELF INDUCTANCE")
    Ls_s = np.zeros(Qs)
    sum = 0
    for slot in range(0,Qs):
        femm.mo_seteditmode('area')
        femm.mo_groupselectblock(slot)
        current,*_ = femm.mo_getcircuitproperties(str(slot))
        if current != 0:
            Ls_s[slot] = femm.mo_blockintegral(0) / (current**2) #p15 self inductance
            sum += Ls_s[slot]
            print("currrent ", current, "Ls_s ", Ls_s[slot], " sum ", sum)
        femm.mo_clearblock()
    Ls_s = sum / np.count_nonzero(Ls_s)
    # mutual inductances
    print("MUTUAL INDUCTANCE")
    femm.mi_saveas("C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_FEMM.bak2\\backup.FEM")
    current = np.zeros(Qs)
    mutuals = np.zeros(Qs)
    for slot in range(0,Qs): # stash circuit currents
        current[slot],*_ = femm.mo_getcircuitproperties(str(slot))
    for slot in range(1,Qs): # for each slot except 0th
        for off in range(0,Qs):
            femm.mi_setcurrent(str(off),0) # turn off all slots
        if current[slot] != 0:
            femm.mi_setcurrent(str(slot),current[slot]) # turn on the next slot
            print("solving with slot ", slot, " only, with amps ", current[slot])

            FEMM_solve(0)
            femm.mi_loadsolution()
            femm.mo_showdensityplot(1,0,3,0,"jmag") 
            femm.mo_seteditmode('area')
            femm.mo_groupselectblock(slot)
            A = femm.mo_blockintegral(1) # A integral
            a = femm.mo_blockintegral(5) # cross section area      

            mutuals[slot] = 0.5/(current[slot]*a)*A
            femm.mo_clearblock()
        else:
            print("skipping slot ", slot, " due to 0 current")

    Ls_m = Ls_s + np.sum(np.absolute(mutuals[1:])) # first element is reference slot
    print("mutuals = ", mutuals)

    # revert to backup and show solution
    femm.opendocument("C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\ISCAD_FEMM.bak2\\backup.FEM")
    femm.mi_probdef(0,'meters','planar',1e-008, l_stack)
    femm.mi_loadsolution()
    femm.mo_hidepoints()
    femm.mo_showdensityplot(1,0,1.8,0,"bmag") 
    femm.mo_resize(1050,800)
    femm.mo_zoomnatural()

    return Ls_s, Ls_m


def DFT(x, N_hmax):
    # digital fourier transform
    # Undertakes a Digital Fourier Transform for harmonic order 1 to N_hmax
    # Uses complex form Xk = 2/N sum[ x(n).exp(-i2(pi)nk/N) ]
    # which yields harmonic amplitudes
    # Uniform spaced samples over range 0 to 360 degrees elect.for fundamental
    N = len(x) # 361
    theta = np.arange(N)*2*np.pi/N # 361 length array, elements 0 to 360
    Xk = np.zeros(N_hmax,dtype=complex) # (1,13) array of 13 zeros
    for k in range(0,N_hmax): # for each element, 0th to 12th index
        for n in range(0,N): # for each element, 0th to 360th index
            # 1st to 13th elements of Xk
            # python index 0 to 12
            Xk[k] = Xk[k] + x[n]*np.exp(-1j*theta[n]*(k+1)) # K+1 !!!!
        Xk[k] = 2*Xk[k]/N
    return Xk

# Bg_k = DFT(B_ag,13) # Harmonic decomposition up to 13th harmonic,
# np.savetxt('Bg_k.txt', np.abs(Bg_k), fmt='%.3f') #


# create model with input parameters, assign block properties inc. coil, save to file
# slots named clockwise by index number converted to string
def FEMM_createmodel(path,params):
    globals().update(params) # HORRIBLE CODE I AM SORRY, cba to rewrite

    # Initialises and defines model units, type, accuracy and length
    print("Opening FEMM...")
    femm.openfemm()
    femm.newdocument(0)
    femm.mi_probdef(0,'meters','planar',1e-008, l_stack,10)

    # Function creates a FEMM model for a stator with 3 independent massive slots
    # non-periodic winding, whole stator modelled

    # Sets up materials and boundaries
    femm.mi_getmaterial('Air') # Gets region materials from library - non-magnetic airgap
    femm.mi_getmaterial('M-15 Steel') # Gets region materials from library - magnetic steel
    femm.mi_addmaterial('Coil') # Adds a region called Coil, defaults to non-magnetic region
    femm.mi_addboundprop('Zero',0,0,0,0,0,0,0,0,0) # Adds a zero potential boundary

    # DRAW CIRCUMFERENCES [yoke, rotor, bore, outer]
    print("Creating circumferences...")
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
    femm.mi_setarcsegmentprop(1,'Zero',0,Qs+Qr+3)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(-r_o,-1)
    femm.mi_setarcsegmentprop(1,'Zero',0,Qs+Qr+3)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(r_yoke,0)
    femm.mi_setarcsegmentprop(1,'Zero',0,Qs+Qr+3)
    femm.mi_clearselected()
    femm.mi_selectarcsegment(-r_yoke,-1)
    femm.mi_setarcsegmentprop(1,'Zero',0,Qs+Qr+3)
    femm.mi_clearselected()
    # Assigns statorsteel region
    femm.mi_addblocklabel(0,(r_o+r_so)/2) # Adds stator region block label
    femm.mi_selectlabel(0,(r_o+r_so)/2) # Selects label
    femm.mi_setblockprop('M-15 Steel', 1, 0, 0, 0, Qs+Qr, 0) # Sets region to magnetic steel
    femm.mi_clearselected()
    # Assigns airgap region
    femm.mi_addblocklabel(0, (r+r_r)/2) # Adds airgap region block label
    femm.mi_selectlabel(0, (r+r_r)/2) # Selects label
    femm.mi_setblockprop('Air', 1, 0, 0, 0, Qs+Qr+1, 0) # sets selected regions to air
    femm.mi_clearselected()
    # Assigns rotor steel region
    femm.mi_addblocklabel(0,(r_r+r_yoke)/2) # Adds stator region block label
    femm.mi_selectlabel(0,(r_r+r_yoke)/2) # Selects label
    femm.mi_setblockprop('M-15 Steel', 1, 0, 0, 0, Qs+Qr, 0) # Sets region to magnetic steel
    femm.mi_clearselected()
    # Assigns yoke zero boundary / no mesh
    femm.mi_addblocklabel((r_yoke)/2,0) # Adds airgap region block label
    femm.mi_selectlabel((r_yoke)/2,0) # Selects label
    femm.mi_setblockprop('<No Mesh>', 1, 0, 0, 0, Qs+Qr+2, 0) # sets selected regions to air
    femm.mi_clearselected()
    

    if r_so != r + h_slot:
        print("slot dimension mismatch! = ", (r_so-(r+h_slot)))
    else:
        print ("slot dimension good")

    # DRAW SLOTS
    print("Drawing slots...")
    R = w_slot / 2
    th = np.atan( ( R[1]-R[0] ) / ( r_so - R[1] - r - R[0]) ) * 180 / np.pi
    # polar coordinates, positive y axis as angle reference
    slot_angle = np.arange(0,2*np.pi,2*np.pi/Qs) # Slot separation angles in radians
    Th_si = np.atan( R[0] / (r+R[0]) ) # inner half pitch angle
    Th_so = np.atan( R[1] / (r_so-R[1]) ) # outer half pitch angle
    # polar coordinates arrays for slot nodes
    R_slot = np.array([ r+R[0] , r_so-R[1] , r_so-R[1] , r+R[0] ])
    Th_slot = np.array([-Th_si,-Th_so,Th_so,Th_si])
    # cartesian coordinates arrays for slot nodes
    x_s = np.multiply( R_slot , np.sin(Th_slot + slot_angle[...,None]) )
    y_s = np.multiply( R_slot , np.cos(Th_slot + slot_angle[...,None]) )    
    for slot in range(0,Qs): # for each slot Qs
        for n in range(0,3,2): # draw 2 radial lines of slot
            femm.mi_drawline( x_s[slot,n] , y_s[slot,n] , x_s[slot,n+1] , y_s[slot,n+1] )
        femm.mi_addarc( x_s[slot,2] , y_s[slot,2] , x_s[slot,1] , y_s[slot,1] , 180 , 1 ) # outer slot arc
        femm.mi_addarc( x_s[slot,0] , y_s[slot,0] , x_s[slot,3] , y_s[slot,3] , 180-(2*th), 1 ) # inner slot arc
        # get coil region centre coords
        x_centre = np.mean(x_s[slot,:])
        y_centre = np.mean(y_s[slot,:])
        # add block labels, assign to circuits
        femm.mi_addblocklabel(x_centre,y_centre)
        femm.mi_selectlabel(x_centre,y_centre)
         # ’blockname’, automesh, meshsize, ’incircuit’, magdir, group, turns
        femm.mi_addcircprop(str(slot), 0, 0) # 0 = parallel connected flag.
        femm.mi_setblockprop('Coil', 1, 0, str(slot), 0, str(slot), 0.5) # assign to circuit str(slot)
        femm.mi_clearselected()

    # DRAW ROTOR BARS
    print("Drawing rotor bars...")
    R = w_bar / 2
    # polar coordinates, positive y axis as angle reference
    bar_angle = np.arange(0,2*np.pi,2*np.pi/Qr) # bar separation angles in radians
    Th_ri = np.atan( R[0] / (r_r - h_bar + R[0]) ) # inner half pitch angle
    Th_ro = np.atan( R[1] / (r_r - R[1]) ) # outer half pitch angle
    # polar coordinates arrays for bar nodes
    R_bar = np.array([ r_r-R[1] , r_r-h_bar+R[0] , r_r-h_bar+R[0] , r_r-R[1] ])
    Th_bar = np.array([-Th_ro,-Th_ri,Th_ri,Th_ro])
    # cartesian coordinates arrays for bar nodes
    x_r = np.multiply( R_bar , np.sin(Th_bar + bar_angle[...,None]) )
    y_r = np.multiply( R_bar , np.cos(Th_bar + bar_angle[...,None]) )    
    for bar in range(0,Qr): # for each bar Qr
        for n in range(0,3,2): # draw 2 radial lines of bar
            femm.mi_drawline( x_r[bar,n] , y_r[bar,n] , x_r[bar,n+1] , y_r[bar,n+1] )
        femm.mi_addarc( x_r[bar,1] , y_r[bar,1] , x_r[bar,2] , y_r[bar,2] , 180 , 1 ) # outer slot arc
        femm.mi_addarc( x_r[bar,3] , y_r[bar,3] , x_r[bar,0] , y_r[bar,0] , 180, 1 ) # inner slot arc
        # get coil region centre coords
        x_centre = np.mean(x_r[bar,:])
        y_centre = np.mean(y_r[bar,:])
        # add block labels, assign to circuits
        femm.mi_addblocklabel(x_centre,y_centre)
        femm.mi_selectlabel(x_centre,y_centre)
         # ’blockname’, automesh, meshsize, ’incircuit’, magdir, group, turns
        # femm.mi_addcircprop(str(bar), 0, 0) # 0 = parallel connected flag.
        femm.mi_setblockprop('Air', 1, 0, 0, 0, str(Qs+bar), 0) # assign to circuit str(slot)
        femm.mi_clearselected()

    femm.mi_zoomnatural() # Displays the model outline to fit window
    print("Done")

    # Saves model - note the directory file structure with '\\' replacing the normal '\'
    # filename = "ISCAD_FEMM.FEM"
    # directory = "C:\\users\\chxps15\\My Documents\\UoBMechElec\\GIP\\FEMM\\"
    # path = directory + filename
    femm.mi_saveas(path)

    print("Steel in group ", Qs+Qr)
    print("Air in group ", Qs+Qr+1)
    print("NoMesh in group ", Qs+Qr+2)
    print("Boundaries in group ", Qs+Qr+3)

    print("Model saved to: ",path)
