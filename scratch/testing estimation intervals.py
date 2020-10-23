# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 19:22:02 2020

@author: eugen
"""
import numpy as np

A = np.array([[0.2,0.3,0.5],\
              [0.3,0.5,0.2],\
              [0.4,0.1,0.5],\
              [0.9,0.1,0.0],\
              [0.3,0.3,0.4],\
              [0.1,0.2,0.7]])

theta = np.array((0.1, 0.2, 0.3))
pi = np.array((0.001, 0.001, 0.1, 0.001, 0.001, 0.001))
sens = 0.99
spec = 0.99
realproby = (1-pi) * np.array((A @ theta)) + pi #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
numsamples = (100 * np.ones(A.shape[0])).astype('int')
ydata =np.random.binomial(numsamples,realprobz)

import Falsification_Sim_EstimationMethods as simEst

out = simEst.Est_BernoulliProjection(A,ydata,numsamples,sens,spec)
out['intProj']
out['endProj']
out['covar']
out['90upper_int']
out['90lower_int']
out['95upper_int']
out['95lower_int']
out['99upper_int']
out['99lower_int']

out2 = simEst.Est_LinearProjection(A,ydata,numsamples,sens,spec)
out2['intProj']
out2['90upper_int']
out2['90lower_int']
out2['95upper_int']
out2['95lower_int']
out2['99upper_int']
out2['99lower_int']
out2['endProj']
out2['90upper_end']
out2['90lower_end']
out2['95upper_end']
out2['95lower_end']
out2['99upper_end']
out2['99lower_end']

intervalStorage = []
for i in range(10000):
    numsamples = (100 * np.ones(A.shape[0])).astype('int')
    ydata = np.random.binomial(numsamples,realprobz)
    out2 = simEst.Est_LinearProjection(A,ydata,numsamples,sens,spec)
    vec = []
    for j in range(A.shape[1]):
        if theta[j] < out2['90upper_int'][j] and theta[j] > out2['90lower_int'][j]:
            vec.append(1)
        else:
            vec.append(0)
    intervalStorage.append(vec)
    
np.mean([intervalStorage[i][0] for i in range(len(intervalStorage))])
np.mean([intervalStorage[i][1] for i in range(len(intervalStorage))])
np.mean([intervalStorage[i][2] for i in range(len(intervalStorage))])
    
    
    
    
    
    