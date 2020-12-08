# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 20:40:09 2020

@author: eugen
"""
import numpy as np
import Falsification_Sim_Modules as simMod
import Falsification_Sim_EstimationMethods as simEst
import scipy.stats as spstat
#CREATING MAP CONFIDENCE INTERVALS
#UNTRACKED
Sens,Spec,wt = 0.95,0.95,0.1
pImp = np.array((0.01, 0.2, 0.1))
pOut = np.array((0.01, 0.01, 0.2, 0.01, 0.1, 0.01))
Q = np.array([[0.5,0.2,0.3],
              [0.4,0.3,0.3],
              [0.1,0.6,0.3],
              [0.2,0.25,0.55],
              [0.15,0.55,0.3],
              [0.3,0.05,0.65]])
realproby = (1-pOut) * np.array((Q @ pImp)) + pOut #optimal testing
realprobz = realproby * Sens + (1-realproby) * (1-Spec) #real testing
(n,m) = Q.shape
beta0 = -4.5*np.ones(n+m)
N = (1000 * np.ones(n)).astype('int')
np.random.seed(90) # Seeds 6, 8, 16, 18, 27 give negative Hessian diagonals
Y = np.random.binomial(N,realprobz)
estimDict = simEst.Est_UntrackedMLE(Q,Y,N,Sens,Spec,wt)
print(np.diag(estimDict['hess']))

hess = estimDict['hess']
#np.diag(hess)
#np.sqrt(np.diag(np.linalg.inv(hess)))

print(estimDict['90upper_int'])
print(estimDict['intProj'])
print(estimDict['90lower_int'])

print(estimDict['90upper_end'])
print(estimDict['endProj'])
print(estimDict['90lower_end'])

'''
print(estimDict['95upper_int'])
print(estimDict['intProj'])
print(estimDict['95lower_int'])

print(estimDict['99upper_int'])
print(estimDict['intProj'])
print(estimDict['99lower_int'])
'''



#TRACKED
Sens,Spec,wt = 0.95,0.95,0.1
pImp = np.array((0.01, 0.2, 0.1))
pOut = np.array((0.01, 0.01, 0.2, 0.01, 0.1, 0.01))
n = len(pOut)
m = len(pImp)
trackP = np.array([pImp]*n)+np.array([1-pImp]*n)*\
            np.array([pOut]*m).transpose()
trackPtildes = Sens*trackP + (1-Spec)*(1-trackP)
beta0 = -4.5*np.ones(n+m)
numSamps = 300
N = np.zeros(shape=(n,m))+numSamps
np.random.seed(16) # Seeds 0,15,16,17,19,20,24,27 give negative Hessian diagonals
Y = np.random.binomial(numSamps,trackPtildes)
estimDict = simEst.Est_TrackedMLE(N,Y,Sens,Spec)
print(np.diag(estimDict['hess']))

hess = estimDict['hess']

print(estimDict['90upper_int'])
print(estimDict['intProj'])
print(estimDict['90lower_int'])

print(estimDict['90upper_end'])
print(estimDict['endProj'])
print(estimDict['90lower_end'])

'''
print(estimDict['95upper_int'])
print(estimDict['intProj'])
print(estimDict['95lower_int'])

print(estimDict['99upper_int'])
print(estimDict['intProj'])
print(estimDict['99lower_int'])
'''



















