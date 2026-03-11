import numpy as np
import csv
import json
from functions import *
from datetime import datetime


# file = 'ISCAD_parameters.json'

# print(type('ISCAD_parameters.json'))
# print(type(file))

# def check(file, foo=0):
#     if file == ('ISCAD_parameters.json'):
#         print("ture")
#     else:
#         print("hee")

# check(file

# from foo2 import *



Qs = 3
Qm_total = 0

Ls_s = 100e-6
arr = [];
print(arr)


arr.append("0")
arr.append(3)



print("self inductances = ", end=" ")
for active_slot in range(0,Qs):
    print("{:.3f}".format(Ls_s*1e3), end=" "), 
print("mH")

# # load parameters from .json
# file = '3slot_parameters.json'
# params = loadjson(file)
# globals().update(params) # update params in current python file 
# derived_params(params, file) # add derived values to params
# globals().update(params) # update params in current python file 


# timestamp = datetime.now().strftime('%H%M%S_%d%m%Y')
# filename = f'python/testcsv_{timestamp}.csv'

# # # with headers
# # with open(filename, 'w', newline='') as csvfile:
# #     fieldnames = ['Qs', 'w_slot']
# #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
# #     writer.writeheader()
# #     # Write a single row with only the specified fields
# #     writer.writerow({key: params[key] for key in fieldnames})


# # all keys
# with open('python/testcsv.csv', 'w', newline='') as csvfile:
#     fieldnames = list(params.keys())  # Get all keys from DotDict
#     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
#     writer.writeheader()
#     writer.writerow(dict(params))  # Convert DotDict to regular dict


# # # without headers
# # with open('python/testcsv.csv', 'w', newline='') as csvfile:
# #     fieldnames = ['parameter', 'value']
# #     writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
# #     # Write each key-value pair as a row (no header)
# #     for key in ['Qs', 'w_slot']:
# #         writer.writerow({'parameter': key, 'value': params[key]})