import numpy as np

x = np.array([(1,2),(3,4),(5,6)])


y = np.array([2,0])
yy = np.array([1,0])

z = (x[0]+y+yy)
# print((x[0]+y+yy))
# print(*(x[0]+y+yy))

a = x[0]
# print(x[0] + x[1])
# print(np.average(x[0] + x[1]))
print(np.mean([x[0], x[1]]))
print(np.average([x[0], x[1]], axis=0))
print(np.mean([x[0], x[1]], axis=0))
print(1 + np.mean([x[0], x[1]], axis=0))
b = np.average(x[0]+x[1])

# file = 'ISCAD_parameters.json'

# print(type('ISCAD_parameters.json'))
# print(type(file))

# def check(file, foo=0):
#     if file == ('ISCAD_parameters.json'):
#         print("ture")
#     else:
#         print("hee")

# check(file)

# from foo2 import *