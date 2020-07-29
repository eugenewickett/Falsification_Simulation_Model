# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:09:23 2020

@author: eugen
"""
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import numpy.random as npr

def PlumleeEstimates(ydata, numsamples, A, sens, spec):
        beta0 = -5 * np.ones(A.shape[1]+A.shape[0])
        
        def invlogit(beta):
            return np.exp(beta)/(np.exp(beta)+1)
        def logit(p):
            return np.log(p/(1-p)) 
        def mynegloglik(beta, ydata, numsamples, A, sens, spec):
            betaI = beta[0:A.shape[1]]
            betaJ = beta[A.shape[1]:]
            probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
            probsz = probs*sens + (1-probs) * (1-spec)
            return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
                + 20 * 1/4*np.sum(np.abs((betaJ - beta[A.shape[1]:]))) #have to regularize to prevent problems
        
        bounds = spo.Bounds(beta0-2, beta0+8)
        opval = spo.minimize(mynegloglik, beta0,
                             args=(ydata, numsamples, A, sens, spec),
                             method='L-BFGS-B',
                             options={'disp': False},
                             bounds=bounds)
        return invlogit(opval.x)[0:A.shape[1]], invlogit(opval.x)[A.shape[1]:]

#test case
A = np.matrix([[0.5, 0.25, 0.25], 
    [0.001, 0.499, 0.5],
   [0.75, 0.001, 0.249],
    [0.499, 0.001, 0.5],
    [0.001, 0.249, 0.75],
    [0.001, 0.001, 0.998]])
pI = np.array((0.001, 0.2, 0.001))
pJ = np.array((0.001, 0.001, 0.2, 0.001, 0.001, 0.001))
sens = 0.99
spec = 0.99
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
numsamples = (100 * np.ones(A.shape[0])).astype('int')
ydata =np.random.binomial(numsamples,realprobz)


# METROPOLIS-HASTINGS
p = 0.10 # initialize w presumed SF rate
lapScale = 1 # initialize w 1; LOOK TO MODIFY
jumpDistStDev = 0.5 # standard deviation for the jumping distribution; LOOK TO MODIFY
x0 = npr.laplace(logit(p),lapScale,A.shape[0]+A.shape[1]) 
xinv = invlogit(x0)
tol = 0.01

ptArr = []
ptArr.append(x0)
numGens = 10000 # how many points to generate?
numRej = 0 # counter for rejections
likelihoodJumpArr = []
for i in range(numGens):
    xCurr = ptArr[-1] #Last submitted point
    xProp = xCurr + npr.normal(loc=0.0,scale=jumpDistStDev,size=len(xCurr)) #Proposed point
    aRatio = mynegloglik(xProp, ydata, numsamples, A, sens, spec) \
         - mynegloglik(xCurr, ydata, numsamples, A, sens, spec)
    if npr.uniform(size=1) < np.exp(aRatio)-tol:
        ptArr.append(xProp)
        xCurr = np.copy(xProp)
        likelihoodJumpArr.append(np.exp(aRatio))
    else:
        ptArr.append(np.copy(xCurr))
        numRej += 1
  
#plot of resulting distribution

testVec = []
for rw in ptArr:
    testVec.append(rw[1])
plt.plot(invlogit(testVec))
plt.hist(invlogit(testVec))
ptArr[9999]

plt.plot(likelihoodJumpArr[1000:])
plt.hist(likelihoodJumpArr[1000:])




















