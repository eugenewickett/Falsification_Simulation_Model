# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:48:57 2020

@author: eugen

This file contains the methods used for estimating SF prevalence from end node
testing results. The inputs to these methdos are:
    1) A:   The estimated transition matrix between m intermediate nodes and n
            end nodes, with n rows and m columns
    2) PosData: A vector of length n containing the respective positive samples
            found.
    3) NumSamples: A vector of length n containing the number of samples (not
            including stockouts) collected at each end node
    4) Sens: Diagnostic sensitivity
    5) Spec: Diagnostic specificity
    6) RglrWt=0.1: Regularization weight (only used for MLE optimization)
    7) M=500, Madapt=5000, delta=0.4: Parameters only used for NUTS sampling
"""

import numpy as np
import scipy.optimize as spo
#import scipy.special as sps
import Falsification_Sim_Modules as simModules

def GetUsableSampleVectors(A,PosData,NumSamples):
    '''
    Takes in vectors of sample amounts, sample positives, and a transition
    matrix A, and returns the same items but suitable for manipulation. Also
    returns a list of two lists containing the [rows],[cols] of removed indices.    
    '''
    n = len(NumSamples)
    m = len(A[0])
    # Grab the zeros lists first
    zeroInds = [[],[]]
    zeroInds[0] = [i for i in range(n) if (NumSamples[i]==0)]
    zeroInds[1] = [i for i in range(m) if (np.sum(A[:,i])==0)]
    
    #Adjust the vectors, doing NumSamples last
    idx = np.argwhere(np.all(A[..., :] == 0, axis=0))
    adjA = np.delete(A, idx, axis=1)
    adjA = np.delete(adjA,zeroInds[0],0)
    adjPosData = [PosData[i] for i in range(n) if (NumSamples[i] > 0)]
    adjNumSamples = [NumSamples[i] for i in range(n) if (NumSamples[i] > 0)]
    
    return adjA, adjPosData, adjNumSamples, zeroInds
    
########################### SF RATE ESTIMATORS ###########################
def Est_LinearProjection(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,\
                         Madapt=5000,delta=0.4): 
    '''
    Linear Projection Estimate: Uses the (estimated) transition matrix, A, 
    and the (estimated) percentage SF at each end node, X, calculaed as PosData
    / NumSamples
    '''
    # Grab 'usable' data
    adjA, adjPosData, adjNumSamples, zeroInds = GetUsableSampleVectors(A,PosData\
                                                                       ,NumSamples)
    
    X = np.array([adjPosData[i]/adjNumSamples[i] for i in range(len(adjNumSamples))])
    intProj = np.dot(np.linalg.inv(np.dot(adjA.T,adjA)),np.dot(adjA.T,X))
    endProj = np.subtract(X,np.dot(adjA,intProj))
    #Insert 'nan' where we didn't have any samples
    for i in range(len(zeroInds[0])):
        endProj = np.insert(endProj,zeroInds[0][i],np.nan)
    for i in range(len(zeroInds[1])):
        intProj = np.insert(intProj,zeroInds[1][i],np.nan)
    return np.ndarray.tolist(intProj.T), np.ndarray.tolist(endProj.T)

def Est_BernoulliProjection(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,\
                            Madapt=5000,delta=0.4):
    '''
    MLE of a Bernoulli variable, using iteratively reweighted least squares;
    see Wikipedia page for notation
    '''
    # Grab 'usable' data
    big_m = A.shape[1]
    adjA, adjPosData, adjNumSamples, zeroInds = GetUsableSampleVectors(A,PosData\
                                                                       ,NumSamples)
    
    A = np.array(adjA)
    X = np.array([adjPosData[i]/adjNumSamples[i] for i in range(len(adjNumSamples))])
    currGap = 10000
    tol = 1e-2
    n = A.shape[0] # Number of end nodes
    m = A.shape[1] # Number of intermediate nodes
    X = np.reshape(X,(n,1))
    w_k = np.zeros([m,1])
    while currGap > tol:
        mu_k = []
        for i in range(n):
            mu_k.append(float(1/(1+np.exp(-1*((np.dot(w_k.T,A[i])))))))
        Sdiag = []
        for i in range(n):
            Sdiag.append(mu_k[i]*(1-mu_k[i]))            
        mu_k = np.reshape(mu_k,(n,1))
        S_k = np.diag(Sdiag)
        w_k1 = np.dot(np.linalg.inv(np.dot(A.T,np.dot(S_k,A))),np.dot(A.T, np.subtract(np.add(np.dot(np.dot(S_k,A),w_k),X),mu_k)))
        if np.linalg.norm(w_k-w_k1) > currGap+tol:
            print('BERNOULLI ALGORITHM COULD NOT CONVERGE')
            intProj = np.zeros([big_m,1])
            endProj = np.zeros([len(NumSamples),1])
            return np.ndarray.tolist(np.squeeze(intProj.T)), np.ndarray.tolist(np.squeeze(endProj.T))
        else:
            currGap = np.linalg.norm(w_k-w_k1)
        w_k = np.copy(w_k1)
    # Now our importer SF rates are calculated; figure out variance + Wald statistics
    covarMat_Bern = np.linalg.inv(np.dot(A.T,np.dot(S_k,A)))
    w_Var = np.diag(covarMat_Bern)
    wald_stats = []
    for j in range(m):
        wald_stats.append(float((w_k[j]**2)/w_Var[j]))
    
    # Convert to intermediate and end node estimates
    intProj = np.ndarray.tolist(simModules.invlogit(w_k.T.tolist()[0]))
    errs_Bern = np.subtract(X,mu_k)
    endProj = errs_Bern.T.tolist()[0]
    #Insert 'nan' where we didn't have any samples
    for i in range(len(zeroInds[0])):
        endProj = np.insert(endProj,zeroInds[0][i],np.nan)
    for i in range(len(zeroInds[1])):
        intProj = np.insert(intProj,zeroInds[1][i],np.nan)
    # Could also return the covariance matrix and Wald statistics if needed
    return intProj, endProj

def Est_MLE_Optimizer(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,\
                      Madapt=5000,delta=0.4):
    '''
    Uses the L-BFGS-B method of the SciPy Optimizer to maximize the
    log-likelihood of different SF rates for a given set of testing data and 
    diagnostic capabilities
    '''
   
    PosData = np.array(PosData)
    NumSamples = np.array(NumSamples)
    n = A.shape[0]
    m = A.shape[1]
    beta0 = -6 * np.ones(m+n) + np.random.uniform(-1,1,m+n)
    #Need to figure out why this likelihood function wasn't working when called
    # via the Modules...
    def mynegloglik_INTERIOR(beta,A,PosData,NumSamples,Sens,Spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-simModules.invlogit(betaJ)) * np.matmul(A,simModules.invlogit(betaI)) + simModules.invlogit(betaJ)
        probsz = probs*Sens + (1-probs) * (1-Spec)
        return -np.sum(PosData * np.log(probsz) + (NumSamples-PosData) * np.log(1-probsz)) \
            + RglrWt*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bds = spo.Bounds(beta0-8, beta0+8)
    opval = spo.minimize(mynegloglik_INTERIOR, beta0,
                         args=(A,PosData,NumSamples,Sens,Spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bds)
    '''
    #Insert 'nan' where we didn't have any samples
    for i in range(len(zeroInds[0])):
        endProj = np.insert(endProj,zeroInds[0][i],np.nan)
    for i in range(len(zeroInds[1])):
        intProj = np.insert(intProj,zeroInds[1][i],np.nan)
    '''
    return simModules.invlogit(opval.x)[0:A.shape[1]].tolist(), simModules.invlogit(opval.x)[A.shape[1]:].tolist()

def Est_NUTS(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,Madapt=5000,delta=0.4):
    '''
    Returns the mean estimate of M NUTS samples, using the Madapt and delta
    parameters and given testing data
    '''
    def postForNUTS(beta):
        """
        Example of a target distribution that could be sampled from using NUTS.
        (Although of course you could sample from it more efficiently)
        Doesn't include the normalizing constant.
        """
        return  simModules.mylogpost(beta,PosData,NumSamples,A,Sens,Spec),\
                simModules.mylogpost_grad(beta,PosData,NumSamples,A,Sens,Spec)

    beta0 = -2 * np.ones(A.shape[1] + A.shape[0])
    samples, _, _ = simModules.nuts6(postForNUTS,M,Madapt,beta0,delta)
    intMeans = [np.mean(simModules.invlogit(samples[:,i])) for i in range(A.shape[1])]
    endMeans = [np.mean(simModules.invlogit(samples[:,A.shape[1]+i])) for i in range(A.shape[0])]
    return intMeans, endMeans


########################### END SF RATE ESTIMATORS ###########################