#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 22:20:36 2021

@author: hargy
"""

import numpy as np
import matplotlib.pyplot as plt

### Random Number Correlation Test ###

N = 10000
rands = []
#rands1 = []
for i in range(0,N):
    a = np.random.rand()
    #b = np.random.rand()
    rands.append(a)
    #rands1.append(b)

mean_x = np.mean(rands)

### Variance
sigma = 0
#sigma1 = 0
for num in range(0,N):
    a = (rands[num] - mean_x)**2
    #b = (rands1[num] - mean_x1)**2
    sigma+=a
    #sigma1+=b

var = sigma/(N-1)
#var1 = sigma1/(N-1)

### Correlation

k_vals = np.linspace(1,N-1,int(N/10))
corr_list = []

#print(k_vals)

for k in k_vals:
    sum_corr = 0
    ki = int(round(k))
    for j in range(0,N):
        c = rands[j%N]*rands[(j+ki)%N]
        sum_corr+=c
    corr_list.append(sum_corr/N)
    
mean_corr = np.mean(corr_list)


print('Mean: '+str(mean_x))
print('Variance: '+str(var))
print('Mean Correlation: '+str(mean_corr))

#print('Correlation: '+ str(sum_corr))


