import numpy as np
from functions import function


'''
circuits = np.array(['R','Y','B'])

for n in range(1,3):
    print(n)

print("liug",len(circuits))
'''
Qs = 3

A = 1
B = 2
C = 3

varlist = [A,B,C]

function(varlist)

Aph = np.zeros(Qs)
print(Aph)