# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:49:53 2020

@author: eugen
"""
import numpy as np
import scipy.optimize as spo

def PlumleeEstimates(ydata, numsamples, A, sens, spec, rglrWt = 0.1):
    ydata = np.array(ydata)
    numsamples = np.array(numsamples)
    beta0 = -6 * np.ones(A.shape[1]+A.shape[0])
    def invlogit_INTERIOR(beta):
        return np.exp(beta)/(np.exp(beta)+1)
    def mynegloglik_INTERIOR(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit_INTERIOR(betaJ)) * np.matmul(A,invlogit_INTERIOR(betaI)) + invlogit_INTERIOR(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
            + rglrWt*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bounds = spo.Bounds(beta0-8, beta0+8)
    opval = spo.minimize(mynegloglik_INTERIOR, beta0+1,
                         args=(ydata, numsamples, A, sens, spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bounds)
    return invlogit_INTERIOR(opval.x)[0:A.shape[1]], invlogit_INTERIOR(opval.x)[A.shape[1]:]



Amat = np.array([[0.6,0.2,0.2],
                 [0.0,0.0,1.0],
                 [0.2,0.7,0.1],
                 [0.2,0.5,0.3],
                 [0.9,0.1,0.0],
                 [0.1,0.5,0.4]])
ydat = [20,10,40,30,10,50]
nsamp = [100,100,100,100,100,100]




iProj, eProj = PlumleeEstimates(ydat,nsamp,Amat,0.99,0.99)
print(iProj)
print(eProj)









