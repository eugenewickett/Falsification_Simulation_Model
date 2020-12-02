# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 20:40:09 2020

@author: eugen
"""
import numpy as np
import Falsification_Sim_Modules as simMod
import Falsification_Sim_EstimationMethods as simEst
#CREATING MAP CONFIDENCE INTERVALS
#UNTRACKED
Sens,Spec,wt = 0.95,0.95,0.1
pImp = np.array((0.001, 0.2, 0.1))
pOut = np.array((0.001, 0.001, 0.2, 0.001, 0.1, 0.001))
Q = np.array([[0.5,0.2,0.3],
              [0.4,0.3,0.3],
              [0.1,0.6,0.3],
              [0.2,0.25,0.55],
              [0.3,0.4,0.3],
              [0.3,0.05,0.65]])
realproby = (1-pOut) * np.array((Q @ pImp)) + pOut #optimal testing
realprobz = realproby * Sens + (1-realproby) * (1-Spec) #real testing
N = (1000 * np.ones(Q.shape[0])).astype('int')
Y = np.random.binomial(N,realprobz)
(n,m) = Q.shape
beta0 = -4.5*np.ones(n+m)


estimDict = simEst.Est_UntrackedMLE(Q,Y,N,Sens,Spec,wt)

hess = estimDict['hess']

print(estimDict['intProj'])
print(estimDict['endProj'])
print(np.linalg.det(hess))

#TRACKED
Sens,Spec,wt = 0.95,0.95,0.1
pImp = np.array((0.001, 0.2, 0.1))
pOut = np.array((0.001, 0.001, 0.2, 0.001, 0.1, 0.001))
n = len(pOut)
m = len(pImp)
trackP = np.array([pImp]*n)+np.array([1-pImp]*n)*\
            np.array([pOut]*m).transpose()
trackPtildes = Sens*trackP + (1-Spec)*(1-trackP)
numSamps = 100
N = np.zeros(shape=(n,m))+numSamps
Y = np.random.binomial(numSamps,trackPtildes)

beta0 = -4.5*np.ones(n+m)

estimDict = simEst.Est_TrackedMLE(N,Y,Sens,Spec)

hess = estimDict['hess']

print(estimDict['intProj'])
print(estimDict['endProj'])
print(np.linalg.det(hess))
for k in range(n+m):
    print(hess[k,k])













