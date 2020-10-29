# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:43:07 2020

@author: eugen
"""

import numpy as np
import random

theta = [0.1,0.2,0.3]
pi = [0.001,0.2,0.05,0.4,0.001,0.001]

m=len(theta)
n=len(pi)

sens=0.95
spec=0.95

numSamp=10000

newdata=[]
for i in range(numSamp):
    currImp = random.randint(0,m-1)
    currEnd = random.randint(0,n-1)
    currP = theta[currImp] + (1-theta[currImp])*pi[currEnd]
    currPtild = sens*currP + (1-spec)*(1-currP)
    sampResult = np.random.binomial(1,currPtild)
    newdata.append([currImp,currEnd,sampResult])

newdata.sort()

# convert data into matrix
N = np.zeros(shape=(n,m))
Y = np.zeros(shape=(n,m))
for samp in newdata:
    j,i,res = samp[0], samp[1], samp[2]
    N[i,j] += 1
    Y[i,j] += res

def newLikeFunc(pVec,numMat, posMat,sens,spec):
    # pVec should be [importers, outlets]
    n,m = numMat.shape
    th=pVec[:m]
    py = pVec[m:]
    pMat = np.zeros(shape=(n,m))
    for i in range(n):
        for j in range(m):
            pMat[i,j] = th[j]+(1-th[j])*py[i]
    pMatTilda = np.zeros(shape=(n,m))
    for i in range(n):
        for j in range(m):
            pMatTilda[i,j] = sens*pMat[i,j] + (1-spec)*(1-pMat[i,j])
    
    L = np.sum(np.multiply(Y,np.log(pMatTilda))+np.multiply(np.subtract(numMat,posMat),\
               np.log(1-pMatTilda)))
    
    return L

#test
pVec1 = [0.05, 0.23, 0.2, 0.0, 0.3, 0.0, 0.6, 0.0, 0.0]
pVec2 = [0.07, 0.23, 0.32, 0.0, 0.25, 0.0, 0.5, 0.0, 0.0]
pVec3 = [0.09, 0.21, 0.29, 0.0, 0.23, 0.0, 0.48, 0.0, 0.0]
pVec4 = [0.1,0.2,0.3, 0.001,0.2,0.05,0.4,0.001,0.001]

newLikeFunc(pVec1,N,Y,sens,spec)
newLikeFunc(pVec2,N,Y,sens,spec)
newLikeFunc(pVec3,N,Y,sens,spec)
newLikeFunc(pVec4,N,Y,sens,spec)

import Falsification_Sim_EstimationMethods as simEst

th,pi = simEst.Est_SampleMLE_Optimizer(newdata,3,6,sens,spec)





