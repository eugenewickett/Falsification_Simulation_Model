# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 15:48:44 2020

@author: eugen
"""

import Falsification_Sim_Modules as simModules
import scipy.optimize as spo
import numpy as np

Sens,Spec,wt = 0.95,0.95,0.1
pImp = np.array((0.001, 0.2, 0.1))
pOut = np.array((0.001, 0.001, 0.2, 0.001, 0.1, 0.001))
n = len(pOut)
m = len(pImp)
trackPtildes = np.zeros(shape=(n,m))
for i in range(n):
    for j in range(m):
        currP = pOut[i]+pImp[j]-pOut[i]*pImp[j]
        trackPtildes[i,j]=Sens*currP+(1-Spec)*(1-currP)

N = np.zeros(shape=(n,m))
Y = np.zeros(shape=(n,m))
numSamps = 100
for i in range(n):
    for j in range(m):
        N[i,j] += numSamps
        Y[i,j] += np.random.binomial(numSamps,trackPtildes[i,j])

beta0 = -3*np.ones(n+m)


import Falsification_Sim_Modules as simModules
L0 = simModules.TRACKED_NegLikeFunc(beta0,N,Y,Sens,Spec,wt)

dL0 = simModules.TRACKED_NegLikeFunc_Jac(beta0,N,Y,Sens,Spec,wt)

for k in range(m+n):
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    
    L1 = simModules.TRACKED_NegLikeFunc(beta1,N,Y,Sens, Spec,wt)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])

bds = spo.Bounds(np.zeros(n+m)-8, np.zeros(m+n)+8)
opVal = spo.minimize(simModules.TRACKED_NegLikeFunc,
                             beta0,
                             args=(N,Y,Sens,Spec,wt),
                             method='L-BFGS-B',
                             options={'disp': False},
                             bounds=bds)
simModules.invlogit(opVal.x)        

pVec,numMat,posMat,sens,spec,RglrWt = opVal.x,N,Y,Sens,Spec,wt




n,m = numMat.shape
th=pVec[:m]
py = pVec[m:]
pMat = np.zeros(shape=(n,m))
betaInitial = -6*np.ones(m+n)
for i in range(n):
    for j in range(m):
        pMat[i,j] = simModules.invlogit(th[j])+(1-simModules.invlogit(th[j]))\
                    *simModules.invlogit(py[i])
pMatTilda = np.zeros(shape=(n,m))
for i in range(n):
    for j in range(m):
        pMatTilda[i,j] = sens*pMat[i,j] + (1-spec)*(1-pMat[i,j])

#Grab importers partials first, then outlets
partialsVec = []

for impInd in range(m):
    currImpPartial = np.sum([posMat[a,impInd]*(1-simModules.invlogit(py[a]))*(sens+spec-1)/pMatTilda[a,impInd]
                            - (numMat[a,impInd]-posMat[a,impInd])*(1-simModules.invlogit(py[a]))*(sens+spec-1)/(1-pMatTilda[a,impInd])                                
                             for a in range(n)])
    partialsVec.append(currImpPartial)
for outInd in range(n):
    if py[outInd] > betaInitial[m+outInd-1]:
        c = 1
    elif py[outInd] < betaInitial[m+outInd-1]:
        c = -1
    else:
        c = 0
    currOutPartial = np.sum([posMat[outInd,b]*(1-simModules.invlogit(th[b]))*(sens+spec-1)/pMatTilda[outInd,b]
                            - (numMat[outInd,b]-posMat[outInd,b])*(1-simModules.invlogit(th[b]))*(sens+spec-1)/(1-pMatTilda[outInd,b])
                            - RglrWt*c                               
                             for b in range(m)])
    partialsVec.append(currOutPartial)

[partialsVec[i]*-1 for i in range(len(partialsVec))]











