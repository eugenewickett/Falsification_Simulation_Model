# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:15:52 2020

@author: FloYd
"""
import numpy
import scipy.optimize as spo


def PlumleeEstimates(ydata, numsamples, A, sens, spec):
    beta0 = -5 * numpy.ones(A.shape[1]+A.shape[0])
    
    def invlogit(beta):
        return numpy.exp(beta)/(numpy.exp(beta)+1)
    def logit(p):
        return numpy.log(p/(1-p)) 
    def mynegloglik(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit(betaJ)) * numpy.array(A @ invlogit(betaI)) + invlogit(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -numpy.sum(ydata * numpy.log(probsz) + (numsamples-ydata) * numpy.log(1-probsz)) \
            + 4 * 1/2*numpy.sum((betaJ - beta0[A.shape[1]:]) ** 2) #have to regularize to prevent problems
    
    bounds = spo.Bounds(beta0-5, beta0+5)
    opval = spo.minimize(mynegloglik, -3 * numpy.ones(A.shape[1]+A.shape[0]),
                         args=(ydata, numsamples, A, sens, spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bounds)
    return invlogit(opval.x)[0:A.shape[1]], invlogit(opval.x)[A.shape[1]:]


#test case
A = numpy.matrix([[0.5, 0.25, 0.25], 
    [0.001, 0.499, 0.5],
    [0.75, 0.001, 0.249],
    [0.499, 0.001, 0.5],
    [0.001, 0.249, 0.75],
    [0.001, 0.001, 0.998]])
pI = numpy.array((0.001, 0.2, 0.001))
pJ = numpy.array((0.001, 0.001, 0.2, 0.001, 0.001, 0.001))
sens = 0.99
spec = 0.99
realproby = (1-pJ) * numpy.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
numsamples = (1000 * numpy.ones(A.shape[0])).astype('int')
ydata =numpy.random.binomial(numsamples,realprobz)

importerhat, outlethat = PlumleeEstimates(ydata, numsamples, A, sens, spec)
print(importerhat)
print(outlethat)