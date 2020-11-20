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

#TRACKED
L0 = simModules.TRACKED_NegLogLikeFunc(beta0,N,Y,Sens,Spec,wt)
dL0 = np.array(simModules.TRACKED_NegLogLikeFunc_Jac(beta0,N,Y,Sens,Spec,wt))

for k in range(m+n):
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    L1 = simModules.TRACKED_NegLogLikeFunc(beta1,N,Y,Sens, Spec,wt)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])

bds = spo.Bounds(np.zeros(n+m)-8, np.zeros(m+n)+8)
opVal = spo.minimize(simModules.TRACKED_NegLogLikeFunc,
                             beta0,
                             args=(N,Y,Sens,Spec,wt),
                             method='L-BFGS-B',jac=simModules.TRACKED_NegLogLikeFunc_Jac,
                             options={'disp': False},
                             bounds=bds)

#jac=simModules.TRACKED_NegLikeFunc_Jac,
simModules.invlogit(opVal.x)        

pVec,numMat,posMat,sens,spec,RglrWt = opVal.x,N,Y,Sens,Spec,wt

#UNTRACKED
Q = np.array([[0.5,0.2,0.3],
              [0.4,0.3,0.3],
              [0.1,0.6,0.3],
              [0.2,0.2,0.6],
              [0.3,0.4,0.3],
              [0.3,0.2,0.5],])
realproby = (1-pOut) * np.array((Q @ pImp)) + pOut #optimal testing
realprobz = realproby * Sens + (1-realproby) * (1-Spec) #real testing
N = (1000 * np.ones(Q.shape[0])).astype('int')
Y = np.random.binomial(N,realprobz)



L0 = simModules.UNTRACKED_NegLogLikeFunc(beta0,N,Y,Sens,Spec,Q,wt)
dL0 = simModules.UNTRACKED_NegLogLikeFunc_Jac(beta0,N,Y,Sens,Spec,Q,wt)
for k in range(m+n):
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    L1 = simModules.UNTRACKED_NegLogLikeFunc(beta1,N,Y,Sens, Spec,Q,wt)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])

bds = spo.Bounds(beta0-8, beta0+8)
opval = spo.minimize(simModules.UNTRACKED_NegLogLikeFunc, beta0,
                     args=(N,Y,Sens,Spec,Q,wt),
                     method='L-BFGS-B',
                     jac=simModules.UNTRACKED_NegLogLikeFunc_Jac,
                     options={'disp': False},
                     bounds=bds)
print(simModules.invlogit(opval.x))
print(pImp)
print(pOut)

#testing UNTRACKED posterior samples
lklhdEst_M, lklhdEst_Madapt, lklhdEst_delta = 500, 5000, 0.4 
postSamps = simModules.GeneratePostSamps_UNTRACKED(N,Y,Q,Sens,Spec,wt,\
                                                  lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta)

import matplotlib.pyplot as plt
np.mean(simModules.invlogit(postSamps[:,2]))
fig = plt.figure()
ax = fig.add_axes([0,0,2,1])
ax.set_xlabel('Intermediate Node',fontsize=16)
ax.set_ylabel('Est. model parameter distribution',fontsize=16)
for i in range(m):
    plt.hist(simModules.invlogit(postSamps[:,i]))

fig = plt.figure()
ax = fig.add_axes([0,0,2,1])
ax.set_xlabel('End Node',fontsize=16)
ax.set_ylabel('Est. model parameter distribution',fontsize=16)
for i in range(n):
    plt.hist(simModules.invlogit(postSamps[:,m+i]))

meanSampVec = []
for i in range(m):
    meanSampVec.append(np.mean(simModules.invlogit(postSamps[:,i])))
for i in range(n):
    meanSampVec.append(np.mean(simModules.invlogit(postSamps[:,i+m])))
meanSampVec = [round(meanSampVec[i],3) for i in range(len(meanSampVec))]



grad1=simModules.mynegloglik_grad(opval.x, N, Y, sens, spec, Q)
grad1_norm = grad1/np.linalg.norm(grad1)
simModules.UNTRACKED_NegLogLikeFunc_Jac(opval.x,N,Y,sens,spec,Q,0)

#testing TRACKED posterior samples
N = np.zeros(shape=(n,m))
Y = np.zeros(shape=(n,m))
numSamps = 100
for i in range(n):
    for j in range(m):
        N[i,j] += numSamps
        Y[i,j] += np.random.binomial(numSamps,trackPtildes[i,j])

lklhdEst_M, lklhdEst_Madapt, lklhdEst_delta = 500, 5000, 0.4 
postSamps = simModules.GeneratePostSamps_TRACKED(N,Y,Sens,Spec,wt,\
                                                  lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta)


fig = plt.figure()
ax = fig.add_axes([0,0,2,1])
ax.set_xlabel('Intermediate Node',fontsize=16)
ax.set_ylabel('Est. model parameter distribution',fontsize=16)
for i in range(m):
    plt.hist(simModules.invlogit(postSamps[:,i]))

fig = plt.figure()
ax = fig.add_axes([0,0,2,1])
ax.set_xlabel('End Node',fontsize=16)
ax.set_ylabel('Est. model parameter distribution',fontsize=16)
for i in range(n):
    plt.hist(simModules.invlogit(postSamps[:,m+i]))

meanSampVec = []
for i in range(m):
    meanSampVec.append(np.mean(simModules.invlogit(postSamps[:,i])))
for i in range(n):
    meanSampVec.append(np.mean(simModules.invlogit(postSamps[:,i+m])))
meanSampVec = [round(meanSampVec[i],3) for i in range(len(meanSampVec))]





    