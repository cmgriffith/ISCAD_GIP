%LTspice_parameter_getter_for_pulses






I_target=140; %Target current to reach/ stabilise at in the RL bar
R_bar=0.853e-3; %Resistance in RL bar
V_in = 48; %Input voltage



Gate_resistance=0.79; % SUBJECT TO CHANGE MOSFET DEPENDENT



R_on = 0.56e-3; %On resistance of mosfet SUBJECT TO CHANGE (USING TYPICAL FOR NOW)
R_high=R_on; %Resistance of high mosfet
R_low = R_on; %Resistance of low mosfet

R_inductor =0; %Inductor resistance think its 0 as it is series RL have a think about it as it will have a value if we dont model as a bar but model as a resistance and inductor instead (in practical work)   

R_total = R_on + R_bar;

display(R_total);


D= (I_target*R_total)/V_in;

display(D);
display(D*100);



Frequency = 1e3; %Frequency
T_per = 1/Frequency; % Time period
T_on = D*T_per; %Time on
T_rise = 1e-9;  %Rise time
T_fall = 1e-9;  %Fall time
T_dead =50e-9;  %Dead time



display(T_on);



%NEED TO CONSIDER HOW GATE RESISTORS ARE AFFECTING CURRENT IN THE RL actually this shouldnt effect that i think it may depend on how fast rising is etc but it should be ok, BUT MAY DO AN IDEAL FIRST AND THEN COME BACK LATER AND MEASURE THE EFFECTS