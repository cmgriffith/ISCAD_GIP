import femm
import numpy as np

# parameters [SI]
# ----------------------USER INPUT-----------------------
Qs = 3 # or 60?
p = 2 # number of pole pairs
m = Qs/p # phase count, also slots per pole pair
Apk = 100 # circuit peak amps

r_bore = 75
r_yoke = 45
r_slot_outer = 86.5 # radius at back of slot [mm]
r_stat_outer = 115 # Stator outer radius [mm]
ag = 0.5
l_stack = 250

w_slot = 3.5 # slot width
h_slot = 10 # slot height


# mu_0 = cons.mu_0 # constants
# -------------------------------------------------------



# used in script
D_o = 2*r_stat_outer
r_o = D_o/2
D = 2*r_bore
r = D/2
D_so = 2*r_slot_outer
r_so = D_so/2
D_r = D - 2*ag
r_r = D_r/2


# FEMM START
femm.openfemm()
femm.newdocument(0)
femm.mi_probdef(0,'millimeters','planar',1e-008, l_stack)

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
femm.mi_selectarcsegment(r_yoke,0)
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
# Assigns yoke air region
femm.mi_addblocklabel(0, (r_yoke)/2) # Adds airgap region block label
femm.mi_selectlabel(0, (r_yoke)/2) # Selects label
femm.mi_setblockprop('Air', 1, 0, 0, 0, 0, 0) # sets selected regions to air
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

circuits = np.array(['R','Y','B'])
Aph = Apk*np.array([1,-0.5,-0.5])
for n in range (0,Qs):
    femm.mi_addcircprop(circuits[n], Aph[n], 1)

# set coil block properties, 
# start at 12 o'clock slot, slots into groups 1,2,3 clockwise
for s in range(0,Qs):
   #  for n in range(0,len(R_slots))
    x_centre = np.mean(x[s,:])
    y_centre = np.mean(y[s,:])
    # add block labels
    femm.mi_addblocklabel(x_centre,y_centre)
    femm.mi_selectlabel(x_centre,y_centre)
    # femm.mi_selectgroup(s+1)
    femm.mi_setblockprop('Coil', 1, 0, circuits[s], 0, (s+1), 0.5)
    print(circuits[s])
    femm.mi_clearselected()


femm.mi_zoomnatural() # Displays the model outline to fit window





'''
Qs = 10
p = 2

print(round(2.500000001))

for k in range(1,round(Qs/(2*p))):
    print(k)

mu_0 = 4*np.pi* 10 ** (-7)


x,y = 2,3

print(x)
print(y)

'''