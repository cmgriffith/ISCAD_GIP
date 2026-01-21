%{
% read parameters into matlab as a struct()
file = "ISCAD_parameters.json";
str = fileread(file);
A = jsondecode(str);
A.machine.Qs
%}

% read parameters into matlab as class properties
params = Analytical.LoadFromJSON('Parameters.json');