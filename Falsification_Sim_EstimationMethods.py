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
import scipy.stats as spstat
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
    # Initialize output dictionary
    outDict = {}
    # Grab 'usable' data
    adjA, adjPosData, adjNumSamples, zeroInds = GetUsableSampleVectors(A,PosData\
                                                                       ,NumSamples)

    X = np.array([adjPosData[i]/adjNumSamples[i] for i in range(len(adjNumSamples))])
    AtA_inv = np.linalg.inv(np.dot(adjA.T,adjA)) # Store so we only calculate once
    intProj = np.dot(AtA_inv,np.dot(adjA.T,X))
    endProj = np.subtract(X,np.dot(adjA,intProj))
    # Generate variance of intermediate projections
    H = np.dot(np.dot(adjA,AtA_inv),adjA.T)
    X_fitted = np.dot(H,X)
    resids = np.subtract(X,X_fitted)
    sampVar = np.dot(resids.T,resids)/(adjA.shape[0]-adjA.shape[1])
    covarInt = sampVar*AtA_inv
    covarInt_diag = np.diag(covarInt)
    varEnds = [sampVar*(1-H[i][i]) for i in range(len(X))]
    t90 = spstat.t.ppf(0.95,adjA.shape[0]-adjA.shape[1])
    t95 = spstat.t.ppf(0.975,adjA.shape[0]-adjA.shape[1])
    t99 = spstat.t.ppf(0.995,adjA.shape[0]-adjA.shape[1])
    int90upper = [intProj[i]+t90*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    int90lower = [intProj[i]-t90*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    int95upper = [intProj[i]+t95*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    int95lower = [intProj[i]-t95*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    int99upper = [intProj[i]+t99*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    int99lower = [intProj[i]-t99*np.sqrt(covarInt_diag[i]) for i in range(len(intProj))]
    end90upper = [endProj[i]+t90*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    end90lower = [endProj[i]-t90*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    end95upper = [endProj[i]+t95*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    end95lower = [endProj[i]-t95*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    end99upper = [endProj[i]+t99*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    end99lower = [endProj[i]-t99*np.sqrt(varEnds[i]) for i in range(len(intProj))]
    #Insert 'nan' where we didn't have any samples
    for i in range(len(zeroInds[0])):
        endProj = np.insert(endProj,zeroInds[0][i],np.nan)
        end90upper = np.insert(end90upper,zeroInds[0][i],np.nan)
        end90lower = np.insert(end90upper,zeroInds[0][i],np.nan)
        end95upper = np.insert(end90upper,zeroInds[0][i],np.nan)
        end95lower = np.insert(end90upper,zeroInds[0][i],np.nan)
        end99upper = np.insert(end90upper,zeroInds[0][i],np.nan)
        end99lower = np.insert(end90upper,zeroInds[0][i],np.nan)
    for i in range(len(zeroInds[1])):
        intProj = np.insert(intProj,zeroInds[1][i],np.nan)
        int90upper = np.insert(int90upper,zeroInds[1][i],np.nan)
        int90lower = np.insert(int90lower,zeroInds[1][i],np.nan)
        int95upper = np.insert(int95upper,zeroInds[1][i],np.nan)
        int95lower = np.insert(int95lower,zeroInds[1][i],np.nan)
        int99upper = np.insert(int99upper,zeroInds[1][i],np.nan)
        int99lower = np.insert(int99lower,zeroInds[1][i],np.nan)
    
    outDict['intProj'] = np.ndarray.tolist(intProj.T)
    outDict['endProj'] = np.ndarray.tolist(endProj.T)
    outDict['covarInt'] = covarInt_diag
    outDict['varEnd'] = varEnds
    outDict['90upper_int'] = int90upper
    outDict['90lower_int'] = int90lower
    outDict['95upper_int'] = int95upper
    outDict['95lower_int'] = int95lower
    outDict['99upper_int'] = int99upper
    outDict['99lower_int'] = int99lower
    outDict['90upper_end'] = end90upper
    outDict['90lower_end'] = end90lower
    outDict['95upper_end'] = end95upper
    outDict['95lower_end'] = end95lower
    outDict['99upper_end'] = end99upper
    outDict['99lower_end'] = end99lower
    return outDict

def Est_BernoulliProjection(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,\
                            Madapt=5000,delta=0.4):
    '''
    MLE of a Bernoulli variable, using iteratively reweighted least squares;
    see Wikipedia page for notation
    '''
    # Initialize output dictionary
    outDict = {}
    
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
    outDict['intProj'] = intProj
    outDict['endProj'] = endProj
    outDict['covar'] = w_Var
    outDict['waldStats'] = wald_stats
    # Confidence intervals: 90%, 95%, 99%
    z90 = spstat.norm.ppf(0.95)
    z95 = spstat.norm.ppf(0.975)
    z99 = spstat.norm.ppf(0.995)
    outDict['90upper_int'] = [simModules.invlogit(w_k[i]+z90*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['90lower_int'] = [simModules.invlogit(w_k[i]-z90*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['95upper_int'] = [simModules.invlogit(w_k[i]+z95*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['95lower_int'] = [simModules.invlogit(w_k[i]-z95*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['99upper_int'] = [simModules.invlogit(w_k[i]+z99*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['99lower_int'] = [simModules.invlogit(w_k[i]-z99*np.sqrt(w_Var[i]))[0] for i in range(m)]
    return outDict

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

def Est_SampleMLE_Optimizer(sampleWiseData,numImp,numOut,Sens,Spec,\
                            RglrWt=0.1,M=500,Madapt=5000,delta=0.4):
    '''
    Forms MLE sample-wise - DOES NOT use A, but instead a list of samples,
    which should be entered [importer index, outlet index, result in {0,1}].
    Then uses the L-BFGS-B method of the SciPy Optimizer to maximize the
    log-likelihood of different SF rates for a given set of testing data and 
    diagnostic capabilities
    '''
    N = np.zeros(shape=(numOut,numImp))
    Y = np.zeros(shape=(numOut,numImp))
    for samp in sampleWiseData:
        j,i,res = samp[0], samp[1], samp[2]
        N[i,j] += 1
        Y[i,j] += res
    def newNegLikeFunc(pVec,numMat, posMat,sens,spec):
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
        # Regularize?
        L = np.sum(np.multiply(Y,np.log(pMatTilda))+np.multiply(np.subtract(numMat,posMat),\
                   np.log(1-pMatTilda))) - RglrWt*np.sum(np.abs(py))
        return L*-1
    p_0 = 0.01 * np.ones(numImp+numOut) + np.random.uniform(-0.001,0.001,numImp+numOut)
    bds = spo.Bounds(p_0-p_0, np.ones(numImp+numOut))
    opVal = spo.minimize(newNegLikeFunc, p_0,
                         args=(N,Y,Sens,Spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bds)
    
    return opVal.x[0:numImp].tolist(), opVal.x[numImp:].tolist()
    

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
    intMeans = [simModules.invlogit(np.mean(samples[:,i])) for i in range(A.shape[1])]
    endMeans = [simModules.invlogit(np.mean(samples[:,A.shape[1]+i])) for i in range(A.shape[0])]
    return intMeans, endMeans


########################### END SF RATE ESTIMATORS ###########################