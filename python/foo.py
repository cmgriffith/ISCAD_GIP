import numpy as np
from scipy import constants as cons

rho_s = 2.8e-08 # resistivity of slot aluminium
f = 25000

x = (np.pi * cons.mu_0 * f)

a = (np.pi * cons.mu_0 * f) / rho_s

b = ( (np.pi * cons.mu_0 * f) / rho_s )**0.5

print(x,a,b)