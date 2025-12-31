'''
% Script to create an AC machine stator stator segment FEMM model
%
% Phil Mellor 2022
%
% A two-pole segment for a multi-pole stator is created.
% The stator dimensions are variables.
% The number of pole pairs (2 or greater) and slots slots per pole pair
% are also input variables
% The model sets up the boundary conditions and region labels,
% however the winding phase pattern (circuits) are not assigned to the coil region
'''
# Input parameters
Stator_outer = 160 # Stator outer diameter Do [mm]
Slot_back_dia = 136 # Diameter at back of slot Ds [mm]
Stator_bore = 100 # Stator bore/active diameter D [mm]
Active_length = 120 # Length of stator L and rotor core packs [mm]
Air_gap = 3 # Stator to rotor air gap clearance [mm]
wtt_pu = 0.667 # Dimensionless slot width, as a proportion of slot pitch [0<wtt_pu<1]
p_pairs = 3 # Number of pole pairs
q_pp = 18 # Number of slots per pole pair, evenly distributed around stator circumference




openfemm()
newdocument(0)
mi_probdef(0,'millimeters','planar',1e-008, Active_length) ; # % Defines model units, type, accuracy and length

# 
def FEMM_createmodel(p_pairs, q_pp, Stator_outer, Slot_back_dia, Stator_bore, Air_gap, wtt_pu):

# Sets up three phase winding 'circuits' and (1, -0.5, -0.5)*I_phase current pattern
I_phase = 80 ; # Phase current magnitude
mi_addcircprop('R',I_phase, 1)
mi_addcircprop('Y', -0.5*I_phase, 1)
mi_addcircprop('B', -0.5*I_phase, 1)
# Saves model - note the directory file structure with '\\' replacing the normal '\'
filename = strcat('Winding_with_',num2str(q_pp/6),'_slots_per_pole_per_phase.FEM')
directory ='C:\Users\\eephm\\Documents\\' 
mi_saveas([directory filename])

function FEMM_dbl_layer_stator(p ,q_pp, Do, Ds, D, Air_gap, wtt_pu)

# Function creates a FEMM model for a stator with 3 independent massive slots
# non-periodic winding, whole stator modelled
#
# p = pole number
# q_pp = slots per pole pair
# Do = stator outer diameter, units as set in FEMM problem definition
# Ds = diameter of back of slot, units as set in FEMM problem definition
# D = stator bore diameter, units as set in FEMM problem definition
# Air_gap = air gap length, units as set in FEMM problem definition
# wtt_pu = tooth tip width as a proportion of the slot pitch at the stator bore


# Sets up materials and boundaries
mi_getmaterial('Air') ; % Gets region materials from library - non-magnetic airgap
mi_getmaterial('M-15 Steel') ; % Gets region materials from library - magnetic steel
mi_addmaterial('Coil') ; % Adds a region called Coil, defaults to non-magnetic region
mi_addboundprop('Zero',0,0,0,0,0,0,0,0,0) ; % Adds a zero potential boundary
mi_addboundprop('SA',0,0,0,0,0,0,0,0,4) ; % Adds a periodic boundary
mi_addboundprop('SB',0,0,0,0,0,0,0,0,4) ; % Adds a second periodic boundary

# Creates outer regions of model, model bounded to a single pole pair
# Outer region co-ordinates defined as an array [xsb, ysb]
pp_angle = 360/p ; % Pole angle in degrees
Dr = D - 2*Air_gap ; % Rotor outer diameter
x_sb = [-Dr, -D, -Do, Do, D, Dr]/2*sin(pp_angle*pi/360) ;
y_sb = [Dr, D, Do, Do, D, Dr]/2*cos(pp_angle*pi/360) ;
mi_addnode(x_sb, y_sb); % Adds nodes
mi_addsegment(x_sb(1), y_sb(1), x_sb(2), y_sb(2)) ; % LH air gap side
mi_addsegment(x_sb(2), y_sb(2), x_sb(3), y_sb(3)) ; % LH stator side
mi_addsegment(x_sb(4), y_sb(4), x_sb(5), y_sb(5)) ; % RH air gap side
mi_addsegment(x_sb(5), y_sb(5), x_sb(6), y_sb(6)) ; % RH stator side
mi_addarc(x_sb(6), y_sb(6), x_sb(1), y_sb(1), pp_angle, 1) ; % Rotor bore
mi_addarc(x_sb(5), y_sb(5), x_sb(2), y_sb(2), pp_angle, 1) ; % Stator bore
mi_addarc(x_sb(4), y_sb(4), x_sb(3), y_sb(3), pp_angle, 1) ; % Stator back
# Sets periodic boundary conditions on sides and zero potential at stator back
mi_selectsegment((x_sb(2)+x_sb(3))/2, (y_sb(2)+y_sb(3))/2) ;
mi_selectsegment((x_sb(4)+x_sb(5))/2, (y_sb(4)+y_sb(5))/2) ;
mi_setsegmentprop('SA',0,0,0,0) ;
mi_clearselected() ;
mi_selectsegment((x_sb(1)+x_sb(2))/2, (y_sb(1)+y_sb(2))/2) ;
mi_selectsegment((x_sb(5)+x_sb(6))/2, (y_sb(5)+y_sb(6))/2) ;
mi_setsegmentprop('SB',0,0,0,0) ;
mi_clearselected() ;
mi_selectarcsegment(0, Do/2) ;
mi_setarcsegmentprop(1,'Zero',0,0) ;
mi_clearselected() ;
# Assigns steel region
mi_addblocklabel(0, (Do+Ds)/4) ; % Adds stator region block label
mi_selectlabel(0, (Do+Ds)/4)
# Selects label
mi_setblockprop('M-15 Steel', 1, 0, 0, 0, 0, 0) ; % Sets region to magnetic steel
mi_clearselected() ;% Assigns airgap region
mi_addblocklabel(0, (Dr+D)/4) ; % Adds stator region block label
mi_selectlabel(0, (Dr+D)/4)
# Selects label
mi_setblockprop('Air', 1, 0, 0, 0, 0, 0) ; % sets selected regions to air
mi_clearselected() ;
mi_zoomnatural() ; % Displays the model outline to fit window

# Creates a polar co-ordinate array for the slot in polar co-ordinates (R_sl, Th_sl)
# Co-ordinates are ordered clockwise around the slot starting at the stator bore
# Profile is for an open slot with a parallel sided tooth
# As the course develops more complex slot profiles can be built

Sl_pitch = 2*pi/(p*q_pp) # Slot pitch angle in radians
Dms = sqrt((D^2 + Ds^2)/2) ; # Mid slot diameter, approx. equal coil side areas
Th_slopen = (1 - wtt_pu)*Sl_pitch/2 ; # Half slot opening angle at stator bore
Th_slmid = Sl_pitch/2 - atan(D/Dms*tan(wtt_pu*Sl_pitch/2) ); # Slot middle half pitch
Th_slback = Sl_pitch/2 - atan(D/Ds*tan(wtt_pu*Sl_pitch/2) ); # Slot back half pitch
N_slco = 6 ; # Number of co-ordinates used to define the slot profile
R_sl = [D, Dms, Ds, Ds, Dms, D]/2 
Th_sl = [-Th_slopen, -Th_slmid, -Th_slback, Th_slback,Th_slmid, Th_slopen]
# Adds the slots as nodes and the slot side surfaces as line segments
# and allocates region properties. Starts at the LHS of the model.
sc_angle = (2*(1:q_pp)-q_pp-1)*(pp_angle/q_pp*pi/360) ; # Slot centre angle

for n = 1:q_pp
x_sl = R_sl.*sin(Th_sl + sc_angle(n))
y_sl = R_sl.*cos(Th_sl + sc_angle(n))
mi_addnode(x_sl, y_sl) ; % Adds slot nodes

for m = 2:N_slco
mi_addsegment(x_sl(m-1), y_sl(m-1), x_sl(m), y_sl(m)) ; % Adds line segment
end
mi_addsegment(x_sl(2), y_sl(2), x_sl(5), y_sl(5)) ; % Splits slot into 2 regions

# Centroid co-ordinates for front of slot region
x_ctf = (x_sl(1)+x_sl(2)+x_sl(5)+x_sl(6))/4 ;
y_ctf = (y_sl(1)+y_sl(2)+y_sl(5)+y_sl(6))/4 ;
# Centroid co-ordinates for back of slot region
x_ctb = (x_sl(2)+x_sl(3)+x_sl(4)+x_sl(5))/4 ;
y_ctb = (y_sl(2)+y_sl(3)+y_sl(4)+y_sl(5))/4 ;
mi_addblocklabel(x_ctf,y_ctf) ; % Adds block label for front coil region
mi_selectlabel(x_ctf,y_ctf) ; % Selects label and assigns material properties
mi_setblockprop('Coil', 1, 0, 0, 0, 2*n-1, 0) ; % Assigns to group number 2*n-1
mi_clearselected();
mi_addblocklabel(x_ctb,y_ctb) ; % Adds block label for back coil region
mi_selectlabel(x_ctb,y_ctb) ; % Selects label and assigns material properties
mi_setblockprop('Coil', 1, 0, 0, 0, 2*n, 0) ; % Assigns to group number 2*n
mi_clearselected();
end
mi_zoomnatural() ; % Sizes model display in FEMM to window
end