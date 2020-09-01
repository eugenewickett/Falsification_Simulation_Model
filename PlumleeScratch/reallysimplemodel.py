# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:15:52 2020

@author: FloYd
"""
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import numpy.random as npr
import os

dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\PlumleeScratch'
os.chdir(dirStr)

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
            + 4 * 1/2*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bounds = spo.Bounds(beta0-2, beta0+8)
    opval = spo.minimize(mynegloglik, beta0+5,
                         args=(ydata, nsamp, A, sens, spec),
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
sens = 0.90
spec = 0.90
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
nsamp = (1000 * np.ones(A.shape[0])).astype('int')
ydata =np.random.binomial(nsamp,realprobz)

importerhat, outlethat = PlumleeEstimates(ydata, nsamp, A, sens, spec)
#print(importerhat)
#print(outlethat)

### CHANGE INITIAL BETA

beta0 = -2 * np.ones(A.shape[1]+A.shape[0])
beta0[0] = -1
beta0[A.shape[1]+2] = -6

def invlogit(beta):
    return np.array(np.exp(beta)/(np.exp(beta)+1))


def invlogit_grad(beta):
    return np.exp(beta)/((np.exp(beta)+1) ** 2)


def logit(p):
    return np.log(p/(1-p))




def mynegloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    
    probs = (1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return -np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz))

def myloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    
    probs = (1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return np.sum(ydata * np.log(probsz) + (nsamp-ydata) * np.log(1-probsz))

L0 = myloglik(beta0,ydata, nsamp, A, sens, spec)

def myloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = ((1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ))
    probs_dirJ = -invlogit_grad(betaJ) * np.squeeze(np.array(A @ invlogit(betaI))) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0])).T
    ilK = np.array([invlogit_grad(betaI),]*(betaJ.shape[0]))
    probs_dirI =  np.multiply(np.multiply((1- ilJ ) , A) , ilK).T
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  np.array(probz_dirI @ (-(ydata / probsz  - (nsamp-ydata) / (1-probsz))).T)
    negloglikJ =  np.array((probz_dirJ * (-(ydata / probsz  - (nsamp-ydata) / (1-probsz)))).T)
    return -np.squeeze(np.concatenate([negloglikI,negloglikJ]))

dL0 = myloglik_grad(beta0,ydata, nsamp, A, sens, spec)
for k in range(0,beta0.shape[0]):
    beta1 = 1*beta0[:]
    beta1[k] = beta1[k] + 10**(-5)
    
    L1 = myloglik(beta1,ydata, nsamp, A, sens, spec)
    print((L1-L0) * (10 **(5)))
    print(dL0[k])



def mylogprior(beta, ydata, nsamp, A, sens, spec):
    betaJ = beta[A.shape[1]:]
    return -0.25*np.sum(np.abs(beta + 3))

def mylogprior_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    return -0.25*np.squeeze(1*(beta >= -3) - 1*(beta <= -3))

L0 = mylogprior(beta0,ydata, nsamp, A, sens, spec)
dL0 = mylogprior_grad(beta0,ydata, nsamp, A, sens, spec)
for k in range(0,beta0.shape[0]):
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




# METROPOLIS-HASTINGS
p = 0.05 # initialize w presumed SF rate
lapScale = 1 # initialize w 1; LOOK TO MODIFY
jumpDistStDev = 0.01 # standard deviation for the jumping distribution; LOOK TO MODIFY
x0 = npr.laplace(logit(p),lapScale,A.shape[0]+A.shape[1]) 
xinv = np.exp(x0)/(1+np.exp(x0))

ptArr = []
ptArr.append(x0)
numGens = 10000 # how many points to generate?
numRej = 0 # counter for rejections
for i in range(numGens):
    xCurr = ptArr[-1] #Last submitted point
    xProp = xCurr + npr.normal(loc=0.0,scale=jumpDistStDev,size=len(xCurr)) #Proposed point
    aRatio = mynegloglik(xProp, ydata, nsamp, A, sens, spec) \
         - mynegloglik(xCurr, ydata, nsamp, A, sens, spec)
    if npr.uniform(size=1) < np.exp(aRatio):
        ptArr.append(xProp)
        xCurr = np.copy(xProp)
    else:
        ptArr.append(np.copy(xCurr))
        numRej += 1
  
#plot of resulting distribution
output = np.array(ptArr)
output.shape
plt.hist(invlogit(output[1,:]))


#If using anaconda, call: conda install -c conda-forge
from nuts import nuts6


D = beta0.shape[0]
M = 1000
Madapt = 1000
delta = 0.2

def exampletargetfornuts(beta):
    """
    Example of a target distribution that could be sampled from using NUTS.
    (Although of course you could sample from it more efficiently)
    Doesn't include the normalizing constant.
    """

    return mylogpost(beta,ydata, nsamp, A, sens, spec), mylogpost_grad(beta,ydata, nsamp, A, sens, spec)

samples, lnprob, epsilon = nuts6(exampletargetfornuts, M, Madapt, opval.x, delta)
plt.hist(invlogit(samples[:,1]))

