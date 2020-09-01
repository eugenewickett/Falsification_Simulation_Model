# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:10:27 2020

@author: eugen
"""
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import os


def PlumleeEstimates(ydata, nsamp, A, sens, spec):
        beta0 = -5 * np.ones(A.shape[1]+A.shape[0])
        
        def invlogit(beta):
            return np.exp(beta)/(np.exp(beta)+1)
        def logit(p):
            return np.log(p/(1-p)) 
        def mynegloglik(beta, ydata, nsamp, A, sens, spec):
            betaI = beta[0:A.shape[1]]
            betaJ = beta[A.shape[1]:]
            probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
            probsz = probs*sens + (1-probs) * (1-spec)
            return -np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz)) \
                + 4 * 1/4*np.sum(np.abs((betaJ - beta[A.shape[1]:]))) #have to regularize to prevent problems
        
        bounds = spo.Bounds(beta0-2, beta0+8)
        opval = spo.minimize(mynegloglik, beta0,
                             args=(ydata, nsamp, A, sens, spec),
                             method='L-BFGS-B',
                             options={'disp': False},
                             bounds=bounds)
        return invlogit(opval.x)[0:A.shape[1]], invlogit(opval.x)[A.shape[1]:]

def invlogit(beta):
    return np.array(np.exp(beta)/(np.exp(beta)+1))

def invlogit_grad(beta):
    return np.diag(np.exp(beta)/((np.exp(beta)+1) ** 2))



def myloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = (1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz))

#def myloglik_grad(beta, ydata, nsamp, A, sens, spec): #EUGENE
#    betaI = beta[0:A.shape[1]]
#    betaJ = beta[A.shape[1]:]
#    probs = (1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ)
#    
#    probsz = probs*sens + (1-probs) * (1-spec)
#    gradVec = ydata/probsz - (nsamp-probsz)/(1-probsz)
#    return gradVec

def myloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
    
    probs_dirJ = -(invlogit_grad(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0]))
    probs_dirI =   (ilJ + ((1- ilJ).T * np.array(A @ invlogit_grad(betaI))).T)
    
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  np.sum((ydata / probsz  - (nsamp-ydata) / (1-probsz)) *  (probz_dirI),1) 
    negloglikJ =  np.sum((ydata / probsz  - (nsamp-ydata) / (1-probsz)) *  (probz_dirJ),1) 

    return np.concatenate((negloglikI,negloglikJ))

def mynegloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
    
    probs_dirJ = -(invlogit_grad(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0]))
    probs_dirI =   (ilJ + ((1- ilJ).T * np.array(A @ invlogit_grad(betaI))).T)
    
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  -np.sum((ydata / probsz  - (nsamp-ydata) / (1-probsz)) *  (probz_dirI),1) 
    negloglikJ =  -np.sum((ydata / probsz  - (nsamp-ydata) / (1-probsz)) *  (probz_dirJ),1) 

    return np.concatenate((negloglikI,negloglikJ))

#def mynegloglik(beta, ydata, nsamp, A, sens, spec):
#    betaI = beta[0:A.shape[1]]
#    betaJ = beta[A.shape[1]:]
#    probs = (1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ)
#    probsz = probs*sens + (1-probs) * (1-spec)
#    return -np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz))

def mynegloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]   
    probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return -np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz))




# TEST DATA
A = np.matrix([[0.5, 0.25, 0.25], 
    [0.001, 0.499, 0.5],
   [0.75, 0.001, 0.249],
    [0.499, 0.001, 0.5],
    [0.001, 0.249, 0.75],
    [0.001, 0.001, 0.998]])
beta0 = -5 * np.ones(A.shape[1]+A.shape[0])
pI = np.array((0.001, 0.2, 0.001))
pJ = np.array((0.001, 0.001, 0.2, 0.001, 0.001, 0.001))
sens = 0.99
spec = 0.99
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
nsamp = (100 * np.ones(A.shape[0])).astype('int')
ydata =np.squeeze(np.random.binomial(nsamp,realprobz))

# Check our likelihood functions are working ok
myloglik(beta0, ydata, nsamp, A, sens, spec)
myloglik_grad(beta0, ydata, nsamp, A, sens, spec)

L0 = mynegloglik(beta0,ydata, nsamp, A, sens, spec)

def mylogprior(beta, ydata, nsamp, A, sens, spec):
    #betaJ = beta[A.shape[1]:]
    return -0.25*np.sum(np.abs(beta + 3))

def mylogprior_grad(beta, ydata, nsamp, A, sens, spec):
    #betaI = beta[0:A.shape[1]]
    #betaJ = beta[A.shape[1]:]
    return -0.25*np.squeeze(1*(beta >= -3) - 1*(beta <= -3))

L0 = mylogprior(beta0,ydata, nsamp, A, sens, spec)
dL0 = mylogprior_grad(beta0,ydata, nsamp, A, sens, spec)
for k in range(0,beta0.shape[0]): # check that the gradient values make sense
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    
    L1 = mylogprior(beta1,ydata, nsamp, A, sens, spec)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])

def mylogpost(beta, ydata, nsamp, A, sens, spec):
    return mylogprior(beta, ydata, nsamp, A, sens, spec)+myloglik(beta, ydata, nsamp, A, sens, spec)

def mylogpost_grad(beta, ydata, nsamp, A, sens, spec):
    return mylogprior_grad(beta, ydata, nsamp, A, sens, spec)+myloglik_grad(beta, ydata, nsamp, A, sens, spec)


L0 = mylogpost(beta0,ydata, nsamp, A, sens, spec)
dL0 = mylogpost_grad(beta0,ydata, nsamp, A, sens, spec)
for k in range(0,beta0.shape[0]):
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    
    L1 = mylogpost(beta1,ydata, nsamp, A, sens, spec)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])
    
    
bounds = spo.Bounds(beta0-2, beta0+8)
opval = spo.minimize(mynegloglik, beta0,
                     args=(ydata, nsamp, A, sens, spec),
                     method='L-BFGS-B',
                     options={'disp': False},
                     bounds=bounds)
opval.x




# NUTS SAMPLER
dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\packages\\NUTS-master'
from dirStr import nuts
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, dirStr)

import nuts

#sampler = NUTSSampler(A.shape[0]+A.shape[1], mynegloglik)
#samples = sampler.run_mcmc( theta0, M, Madapt, delta )

D = beta0.shape[0]
M = 100
Madapt = 10
delta = 0.2

def exampletargetfornuts(beta):
    """
    Example of a target distribution that could be sampled from using NUTS.
    (Although of course you could sample from it more efficiently)
    Doesn't include the normalizing constant.
    """

    return mylogpost(beta,ydata, nsamp, A, sens, spec), mylogpost_grad(beta,ydata, nsamp, A, sens, spec)

samples, lnprob, epsilon = nuts6(exampletargetfornuts, M, Madapt, beta0, delta)
plt.hist(invlogit(samples[:,1]))
    



