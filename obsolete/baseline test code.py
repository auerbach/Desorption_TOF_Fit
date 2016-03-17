# -*- coding: utf-8 -*-

import numpy as np

sig = np.arange(30.,100.)
t = np.arange(3.0,10.0, .1)

r1 = np.logical_and(t>=3.5, t<=4.0)
r2 = np.logical_and(t>=9.0, t<=9.5)
r3 = np.logical_or(r1, r2)

#print('r1 = ', r1)
#print('r2 = ', r2)
# print('r3 = ', r3)
print('sig[r3]= ', sig[r3])

range =(3.5, 4.0, 9.0, 9.5) 
n_segments = int(len(range) / 2)
avg_segments = np.array(range).reshape(n_segments, 2)

first = True
for r in avg_segments:
    c = np.logical_and(t >= r[0], t <=r [1])
    if first:
        cc = c
        first = False
    else: 
        cc = np.logical_or(cc, c)
        
print('sig[cc]= ', sig[cc])
