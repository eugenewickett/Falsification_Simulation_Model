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
import scipy.special as sps

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
    end90upper = [endProj[i]+t90*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    end90lower = [endProj[i]-t90*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    end95upper = [endProj[i]+t95*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    end95lower = [endProj[i]-t95*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    end99upper = [endProj[i]+t99*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    end99lower = [endProj[i]-t99*np.sqrt(varEnds[i]) for i in range(len(endProj))]
    #Insert 'nan' where we didn't have any samples
    for i in range(len(zeroInds[0])):
        endProj = np.insert(endProj,zeroInds[0][i],np.nan)
        end90upper = np.insert(end90upper,zeroInds[0][i],np.nan)
        end90lower = np.insert(end90lower,zeroInds[0][i],np.nan)
        end95upper = np.insert(end95upper,zeroInds[0][i],np.nan)
        end95lower = np.insert(end95lower,zeroInds[0][i],np.nan)
        end99upper = np.insert(end99upper,zeroInds[0][i],np.nan)
        end99lower = np.insert(end99lower,zeroInds[0][i],np.nan)
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
            #print('BERNOULLI ALGORITHM COULD NOT CONVERGE')
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
    intProj = np.ndarray.tolist(sps.expit(w_k.T.tolist()[0]))
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
    outDict['90upper_int'] = [sps.expit(w_k[i]+z90*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['90lower_int'] = [sps.expit(w_k[i]-z90*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['95upper_int'] = [sps.expit(w_k[i]+z95*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['95lower_int'] = [sps.expit(w_k[i]-z95*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['99upper_int'] = [sps.expit(w_k[i]+z99*np.sqrt(w_Var[i]))[0] for i in range(m)]
    outDict['99lower_int'] = [sps.expit(w_k[i]-z99*np.sqrt(w_Var[i]))[0] for i in range(m)]
    return outDict

def Est_UntrackedMLE(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,\
                      Madapt=5000,delta=0.4,beta0_List=[]):
    '''
    Uses the L-BFGS-B method of the SciPy Optimizer to maximize the
    log-likelihood of different SF rates for a given set of UNTRACKED testing
    data in addition to diagnostic capabilities
    '''
    outDict = {} # Output dictionary
    PosData = np.array(PosData)
    NumSamples = np.array(NumSamples)
    numOut = A.shape[0]
    numImp = A.shape[1]
    if beta0_List == []: # We do not have any initial points to test; generate a generic initial point
        beta0_List.append(-6 * np.ones(numImp+numOut) + np.random.uniform(-1,1,numImp+numOut))
    
    #Loop through each possible initial point and store the optimal solution likelihood values
    likelihoodsList = []
    solsList = []
    bds = spo.Bounds(np.zeros(numImp+numOut)-8, np.zeros(numImp+numOut)+8)
    for curr_beta0 in beta0_List:
        opVal = spo.minimize(UNTRACKED_NegLogLikeFunc,
                             curr_beta0,
                             args=(NumSamples,PosData,Sens,Spec,A,RglrWt),
                             method='L-BFGS-B',
                             jac=UNTRACKED_NegLogLikeFunc_Jac,
                             options={'disp': False},
                             bounds=bds)
        likelihoodsList.append(opVal.fun)
        solsList.append(opVal.x)
    best_x = solsList[np.argmin(likelihoodsList)]
    
    #Generate confidence intervals
    #First we need to generate the information matrix
    #Expected positives vector at the outlets
    pi_hat = sps.expit(best_x[numImp:])
    theta_hat = sps.expit(best_x[:numImp])
    #y_Expec = (1-Spec) + (Sens+Spec-1) *(pi_hat + (1-pi_hat)*(A @ theta_hat))
    #Insert it into our hessian
    hess = UNTRACKED_NegLogLikeFunc_Hess(best_x,NumSamples,PosData,\
                                                       Sens,Spec,A)
    
    hess_invs = [i if i >= 0 else np.nan for i in 1/np.diag(hess)]
    z90 = spstat.norm.ppf(0.95)
    z95 = spstat.norm.ppf(0.975)
    z99 = spstat.norm.ppf(0.995)
    imp_Interval90 = z90*np.sqrt(hess_invs[:numImp])
    imp_Interval95 = z95*np.sqrt(hess_invs[:numImp])
    imp_Interval99 = z99*np.sqrt(hess_invs[:numImp])
    out_Interval90 = z90*np.sqrt(hess_invs[numImp:])
    out_Interval95 = z95*np.sqrt(hess_invs[numImp:])
    out_Interval99 = z99*np.sqrt(hess_invs[numImp:])
    
    outDict['90upper_int'] = sps.expit(best_x[:numImp] + imp_Interval90)
    outDict['90lower_int'] = sps.expit(best_x[:numImp] - imp_Interval90)
    outDict['95upper_int'] = sps.expit(best_x[:numImp] + imp_Interval95)
    outDict['95lower_int'] = sps.expit(best_x[:numImp] - imp_Interval95)
    outDict['99upper_int'] = sps.expit(best_x[:numImp] + imp_Interval99)
    outDict['99lower_int'] = sps.expit(best_x[:numImp] - imp_Interval99)
    outDict['90upper_end'] = sps.expit(best_x[numImp:] + out_Interval90)
    outDict['90lower_end'] = sps.expit(best_x[numImp:] - out_Interval90)
    outDict['95upper_end'] = sps.expit(best_x[numImp:] + out_Interval95)
    outDict['95lower_end'] = sps.expit(best_x[numImp:] - out_Interval95)
    outDict['99upper_end'] = sps.expit(best_x[numImp:] + out_Interval99)
    outDict['99lower_end'] = sps.expit(best_x[numImp:] - out_Interval99)
    
    outDict['intProj'] = theta_hat
    outDict['endProj'] = pi_hat
    outDict['hess'] = hess  
    
    return outDict


def Est_TrackedMLE(N,Y,Sens,Spec,RglrWt=0.1,M=500,Madapt=5000,delta=0.4,beta0_List=[]):
    '''
    Forms MLE sample-wise - DOES NOT use A, but instead matrices N and Y,
    which record the positives and number of tests for each (outlet,importer) combination.
    Then uses the L-BFGS-B method of the SciPy Optimizer to maximize the
    log-likelihood of different SF rates for a given set of testing data and 
    diagnostic capabilities
    '''
    outDict = {}
    (numOut,numImp) = N.shape  
    if beta0_List == []: # We do not have any initial points to test; generate a generic initial point
        beta0_List.append(-6 * np.ones(numImp+numOut) + np.random.uniform(-1,1,numImp+numOut))
    
    #Loop through each possible initial point and store the optimal solution likelihood values
    likelihoodsList = []
    solsList = []
    bds = spo.Bounds(np.zeros(numImp+numOut)-8, np.zeros(numImp+numOut)+8)
    for curr_beta0 in beta0_List:
        opVal = spo.minimize(TRACKED_NegLogLikeFunc,
                             curr_beta0,
                             args=(N,Y,Sens,Spec,RglrWt),
                             method='L-BFGS-B',
                             jac=TRACKED_NegLogLikeFunc_Jac,
                             options={'disp': False},
                             bounds=bds)
        likelihoodsList.append(opVal.fun)
        solsList.append(opVal.x)
    best_x = solsList[np.argmin(likelihoodsList)]
    
    #Generate confidence intervals
    #First we need to generate the information matrix
    #Expected positives vector at the outlets
    pi_hat = sps.expit(best_x[numImp:])
    theta_hat = sps.expit(best_x[:numImp])
    #y_Expec = (1-Spec) + (Sens+Spec-1) *(np.array([theta_hat]*numOut)+np.array([1-theta_hat]*numOut)*np.array([pi_hat]*numImp).transpose())
    #Insert it into our hessian
    hess = TRACKED_NegLogLikeFunc_Hess(best_x,N,Y,Sens,Spec)
    
    hess_invs = [i if i >= 0 else np.nan for i in 1/np.diag(hess)]
    z90 = spstat.norm.ppf(0.95)
    z95 = spstat.norm.ppf(0.975)
    z99 = spstat.norm.ppf(0.995)
    imp_Interval90 = z90*np.sqrt(hess_invs[:numImp])
    imp_Interval95 = z95*np.sqrt(hess_invs[:numImp])
    imp_Interval99 = z99*np.sqrt(hess_invs[:numImp])
    out_Interval90 = z90*np.sqrt(hess_invs[numImp:])
    out_Interval95 = z95*np.sqrt(hess_invs[numImp:])
    out_Interval99 = z99*np.sqrt(hess_invs[numImp:])
    
    outDict['90upper_int'] = sps.expit(best_x[:numImp] + imp_Interval90)
    outDict['90lower_int'] = sps.expit(best_x[:numImp] - imp_Interval90)
    outDict['95upper_int'] = sps.expit(best_x[:numImp] + imp_Interval95)
    outDict['95lower_int'] = sps.expit(best_x[:numImp] - imp_Interval95)
    outDict['99upper_int'] = sps.expit(best_x[:numImp] + imp_Interval99)
    outDict['99lower_int'] = sps.expit(best_x[:numImp] - imp_Interval99)
    outDict['90upper_end'] = sps.expit(best_x[numImp:] + out_Interval90)
    outDict['90lower_end'] = sps.expit(best_x[numImp:] - out_Interval90)
    outDict['95upper_end'] = sps.expit(best_x[numImp:] + out_Interval95)
    outDict['95lower_end'] = sps.expit(best_x[numImp:] - out_Interval95)
    outDict['99upper_end'] = sps.expit(best_x[numImp:] + out_Interval99)
    outDict['99lower_end'] = sps.expit(best_x[numImp:] - out_Interval99)
    
    outDict['intProj'] = theta_hat
    outDict['endProj'] = pi_hat
    outDict['hess'] = hess
    
    return outDict
    

def Est_PostSamps_Untracked(A,PosData,NumSamples,Sens,Spec,RglrWt=0.1,M=500,Madapt=5000,delta=0.4):
    '''
    Returns the mean estimate of M NUTS samples, using the Madapt and delta
    parameters and given testing data
    '''
    samples = GeneratePostSamps_UNTRACKED(NumSamples,PosData,A,Sens,Spec,RglrWt,M,Madapt,delta)
    intMeans = [sps.expit(np.mean(samples[:,i])) for i in range(A.shape[1])]
    endMeans = [sps.expit(np.mean(samples[:,A.shape[1]+i])) for i in range(A.shape[0])]
    return intMeans, endMeans

def Est_PostSamps_Tracked(Nmat,Ymat,Sens,Spec,RglrWt=0.1,M=500,Madapt=5000,delta=0.4):
    '''
    Returns the mean estimate of M NUTS samples, using the Madapt and delta
    parameters and given testing data
    '''
    samples = GeneratePostSamps_TRACKED(Nmat,Ymat,Sens,Spec,RglrWt,M,Madapt,delta)
    intMeans = [sps.expit(np.mean(samples[:,i])) for i in range(Nmat.shape[1])]
    endMeans = [sps.expit(np.mean(samples[:,Nmat.shape[1]+i])) for i in range(Nmat.shape[0])]
    return intMeans, endMeans
########################### END SF RATE ESTIMATORS ###########################
    
def UNTRACKED_NegLogLikeFunc(betaVec,numVec,posVec,sens,spec,transMat,RglrWt):
    # betaVec should be [importers, outlets]
    n,m = transMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    betaInitial = -6*np.ones(m+n)
    pVec = sps.expit(py)+(1-sps.expit(py))*np.matmul(transMat,sps.expit(th))
    pVecTilde = sens*pVec + (1-spec)*(1-pVec)
    
    L = np.sum(np.multiply(posVec,np.log(pVecTilde))+np.multiply(np.subtract(numVec,posVec),\
               np.log(1-pVecTilde))) - RglrWt*np.sum(np.abs(py-betaInitial[m:]))
    return L*-1

def UNTRACKED_NegLogLikeFunc_Jac(betaVec,numVec,posVec,sens,spec,transMat,RglrWt):
    # betaVec should be [importers, outlets]
    n,m = transMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    betaInitial = -6*np.ones(m+n)
    pVec = sps.expit(py)+(1-sps.expit(py))*np.matmul(transMat,sps.expit(th))
    pVecTilde = sens*pVec + (1-spec)*(1-pVec)
    
    #Grab importers partials first, then outlets
    impPartials = np.sum(posVec[:,None]*transMat*(sps.expit(th)-sps.expit(th)**2)*(sens+spec-1)*\
                     np.array([(1-sps.expit(py))]*m).transpose()/pVecTilde[:,None]\
                     - (numVec-posVec)[:,None]*transMat*(sps.expit(th)-sps.expit(th)**2)*(sens+spec-1)*\
                     np.array([(1-sps.expit(py))]*m).transpose()/(1-pVecTilde)[:,None]\
                     ,axis=0)
    outletPartials = posVec*(1-np.matmul(transMat,sps.expit(th)))*(sps.expit(py)-sps.expit(py)**2)*\
                        (sens+spec-1)/pVecTilde - (numVec-posVec)*(sps.expit(py)-sps.expit(py)**2)*\
                        (sens+spec-1)*(1-np.matmul(transMat,sps.expit(th)))/(1-pVecTilde)\
                        - RglrWt*np.squeeze(1*(py >= betaInitial[m:]) - 1*(py <= betaInitial[m:]))
    
    retVal = np.concatenate((impPartials,outletPartials))*-1
    
    return retVal

def UNTRACKED_NegLogLikeFunc_Hess(betaVec,numVec,posVec,sens,spec,transMat):
    # betaVec should be [importers, outlets]
    n,m = transMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    
    zVec = sps.expit(py)+(1-sps.expit(py))*np.matmul(transMat,sps.expit(th))
    zVecTilde = sens*zVec+(1-spec)*(1-zVec)
    sumVec = np.matmul(transMat,sps.expit(th))
    
    #initialize a Hessian matrix
    hess = np.zeros((n+m,n+m))
    # get off-diagonal entries first; importer-outlet entries
    for triRow in range(n):
        for triCol in range(m):
            outBeta,impBeta = py[triRow],th[triCol]
            outP,impP = sps.expit(outBeta),sps.expit(impBeta)
            s,r=sens,spec
            c1 = transMat[triRow,triCol]*(s+r-1)*(sps.expit(impBeta)-sps.expit(impBeta)**2)
            yDat,nSam = posVec[triRow],numVec[triRow]
            elem = c1*(1-outP)*(yDat*( (s+r-1)*(-sumVec[triRow]*(outP**2-outP) - outP + outP**2) )\
                    /( s*(sumVec[triRow]*(1 - outP) + outP) +\
                   (1-r)*(-sumVec[triRow]*(1 - outP) + 1 - outP) )**2 -\
                    (nSam - yDat)*((-r + 1-s)*(-sumVec[triRow]*(-outP + outP**2)-outP+outP**2))\
                     /(-s*(sumVec[triRow]*(1 - outP) + outP) - (1-r)*(-sumVec[triRow]*(1 - outP) +\
                   1 - outP) + 1)**2) +\
                    c1*(yDat/(s*(sumVec[triRow]*(1 - outP) + outP) + (-r + 1)*(-sumVec[triRow]*(1 - outP) +\
                   1 - outP)) - (nSam - yDat)/( -s*(sumVec[triRow]*(1 - outP) +\
                   outP) - (1-r)*(-sumVec[triRow]*(1 - outP) + 1 - outP) + 1))*( outP**2 - outP)
            hess[m+triRow,triCol] = elem
            hess[triCol,m+triRow] = elem
    # get off-diagonals for importer-importer entries
    for triCol in range(m-1):
        for triCol2 in range(triCol+1,m):
            elem = 0
            for i in range(n):
                nextPart = (sens+spec-1)*transMat[i,triCol]*(1-sps.expit(py[i]))*(sps.expit(th[triCol])-sps.expit(th[triCol])**2)*\
                (-posVec[i]*(sens+spec-1)*(1-sps.expit(py[i]))*transMat[i,triCol2]*(sps.expit(th[triCol2]) - sps.expit(th[triCol2])**2)            /\
                 (zVecTilde[i]**2)
                - (numVec[i]-posVec[i])*(sens+spec-1)*(1-sps.expit(py[i]))*transMat[i,triCol2]*(sps.expit(th[triCol2]) - sps.expit(th[triCol2])**2) /\
                ((1-zVecTilde[i])**2) )
                
                elem += nextPart
            hess[triCol,triCol2] = elem
            hess[triCol2,triCol] = elem
    # importer diagonals next
    impPartials = np.zeros(m)
    for imp in range(m):
        currPartial = 0
        for outlet in range(n):
            outBeta,impBeta = py[outlet],th[imp]
            outP,impP = sps.expit(outBeta),sps.expit(impBeta)
            s,r=sens,spec                      
            c1 = transMat[outlet,imp]*(s+r-1)*(1-outP)            
            c3 = (1-outP)*transMat[outlet,imp]
            yDat,nSam = posVec[outlet],numVec[outlet]
            currElem = c1*(yDat/(zVecTilde[outlet]) - (nSam - yDat)/(1-zVecTilde[outlet]))\
                       *(impP - 3*(impP**2) + 2*(impP**3)) +\
                       c1*(impP - impP**2)*(yDat*((s+r-1)*c3*(\
                       (impP**2)-impP) )/(zVecTilde[outlet])**2 -\
                       (nSam - yDat)*((s+r-1)*(c3*impP - c3*(impP**2)))/\
                       (1-zVecTilde[outlet])**2)
            currPartial += currElem
        impPartials[imp] = currPartial
    
    # outlet diagonals next
    outletPartials = np.zeros(n)
    for outlet in range(n):
        outBeta = py[outlet]
        outP = sps.expit(outBeta)
        s,r=sens,spec
        c1 = transMat[outlet] @ sps.expit(th)
        c2 = (r + s - 1)
        yDat,nSam = posVec[outlet],numVec[outlet]
        currPartial = (1-c1)*(yDat/(zVecTilde[outlet]) -\
                    (nSam - yDat)/(1-zVecTilde[outlet]))*c2*(outP -\
                    3*(outP**2) + 2*(outP**3)) + \
                      (1-c1)*(outP - outP**2 )*(yDat*(-c2*(c1*(-outP + outP**2 )+ outP -outP**2 ) )/\
                    (zVecTilde[outlet])**2 - (nSam - yDat)*(c2*(c1*(-outP + outP**2) +\
                     outP - outP**2 ))/( -s*(c1*(1 - outP) +\
                     outP) - (1-r)*(1-c1*(1 - outP)  - outP) + 1 )**2)*c2
        outletPartials[outlet] = currPartial
    
    diags = np.diag(np.concatenate((impPartials,outletPartials)))
    
    hess = (hess + diags)*-1 + 0.2 #* np.eye(hess.shape[0]) # include prior term   
    return hess

def UNTRACKED_LogPrior(beta,numVec,posVec,sens,spec,transMat):
    '''
    Prior is a Laplace distribution with parameters: mu=-3, scale=4
    '''
    #-0.25*np.sum(np.abs(beta + 3)) - 0.001 * np.sum((beta + 3) ** 2)
    return - 0.1 * np.sum((beta + 3) ** 2)


def UNTRACKED_LogPrior_Grad(beta, nsamp, ydata, sens, spec, A):
    '''
    Prior is a Laplace distribution with parameters: mu=-3, scale=4
    '''
    #-0.25*np.squeeze(1*(beta >= -3) - 1*(beta <= -3)) - 0.002 * np.sum(beta + 3)
    return - 0.2 * np.sum(beta + 3)


def UNTRACKED_LogPost(betaVec,numVec,posVec,sens,spec,transMat):
    return UNTRACKED_LogPrior(betaVec,numVec,posVec,sens,spec,transMat)\
           -UNTRACKED_NegLogLikeFunc(betaVec,numVec,posVec,sens,spec,transMat,0)

def UNTRACKED_LogPost_Grad(beta, nsamp, ydata, sens, spec, A):
    return UNTRACKED_LogPrior_Grad(beta, nsamp, ydata, sens, spec, A)\
           -UNTRACKED_NegLogLikeFunc_Jac(beta,nsamp,ydata,sens,spec,A,0)

def GeneratePostSamps_UNTRACKED(numSamples,posData,A,sens,spec,regWt,M,Madapt,delta):
    def UNTRACKEDtargetForNUTS(beta):
        return UNTRACKED_LogPost(beta,numSamples,posData,sens,spec,A),\
               UNTRACKED_LogPost_Grad(beta,numSamples,posData,sens,spec,A)

    beta0 = -2 * np.ones(A.shape[1] + A.shape[0])
    samples, lnprob, epsilon = nuts6(UNTRACKEDtargetForNUTS,M,Madapt,beta0,delta)
    
    return samples
###### END UNTRACKED FUNCTIONS ######
    
###### BEGIN UNTRACKED FUNCTIONS ######
def TRACKED_NegLogLikeFunc(betaVec,numMat,posMat,sens,spec,RglrWt):
    # betaVec should be [importers, outlets]
    n,m = numMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    betaInitial = -6*np.ones(m+n)
    pMat = np.array([sps.expit(th)]*n)+np.array([(1-sps.expit(th))]*n)*\
            np.array([sps.expit(py)]*m).transpose()
    pMatTilde = sens*pMat+(1-spec)*(1-pMat)
    
    L = np.sum(np.multiply(posMat,np.log(pMatTilde))+np.multiply(np.subtract(numMat,posMat),\
               np.log(1-pMatTilde))) - RglrWt*np.sum(np.abs(py-betaInitial[m:]))
    return L*-1

def TRACKED_NegLogLikeFunc_Jac(betaVec,numMat,posMat,sens,spec,RglrWt):
    # betaVec should be [importers, outlets]
    n,m = numMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    betaInitial = -6*np.ones(m+n)
    pMat = np.array([sps.expit(th)]*n)+np.array([(1-sps.expit(th))]*n)*\
            np.array([sps.expit(py)]*m).transpose()
    pMatTilde = sens*pMat+(1-spec)*(1-pMat)
    
    #Grab importers partials first, then outlets
    impPartials = np.sum(posMat*(sps.expit(th)-sps.expit(th)**2)*(sens+spec-1)*\
                     np.array([(1-sps.expit(py))]*m).transpose()/pMatTilde\
                     - (numMat-posMat)*(sps.expit(th)-sps.expit(th)**2)*(sens+spec-1)*\
                     np.array([(1-sps.expit(py))]*m).transpose()/(1-pMatTilde)\
                     ,axis=0)
    outletPartials = np.sum((sens+spec-1)*(posMat*(sps.expit(py)-sps.expit(py)**2)[:,None]*\
                     np.array([(1-sps.expit(th))]*n)/pMatTilde\
                     - (numMat-posMat)*(sps.expit(py)-sps.expit(py)**2)[:,None]*\
                     np.array([(1-sps.expit(th))]*n)/(1-pMatTilde))\
                     ,axis=1) - RglrWt*np.squeeze(1*(py >= betaInitial[m:]) - 1*(py <= betaInitial[m:]))
       
    retVal = np.concatenate((impPartials,outletPartials))*-1
    
    return retVal

def TRACKED_NegLogLikeFunc_Hess(betaVec,numMat,posMat,sens,spec):
    # betaVec should be [importers, outlets]
    n,m = numMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    
    zMat = np.array([sps.expit(th)]*n)+np.array([(1-sps.expit(th))]*n)*\
            np.array([sps.expit(py)]*m).transpose()
    zMatTilde = sens*zMat+(1-spec)*(1-zMat)
    
    hess = np.zeros((n+m,n+m))
    # get off-diagonal entries first
    for triRow in range(n):
        for triCol in range(m):
            outBeta,impBeta = py[triRow],th[triCol]
            outP,impP = sps.expit(outBeta),sps.expit(impBeta)
            s,r=sens,spec
            z = outP + impP - outP*impP
            zTilde = zMatTilde[triRow,triCol]
            yDat,nSam = posMat[triRow,triCol],numMat[triRow,triCol]
            elem = (1-impP)*(outP - outP**2)*(yDat*((1-r-s)*(impP-impP**2)*(1-outP))/\
                    zTilde**2-(nSam-yDat)*((s+r-1)*(impP-impP**2-outP*impP+outP*\
                    (impP**2)))/(1-zTilde)**2)*\
                    (r+s-1) + (yDat/zTilde - (nSam - yDat)/(1-zTilde ))\
                    *(outP - outP**2)*(impP**2 -impP)*(r + s - 1)
            hess[m+triRow,triCol] = elem
            hess[triCol,m+triRow] = elem
    
    # importer diagonals next
    impPartials = np.zeros(m)
    for imp in range(m):
        currPartial = 0
        for outlet in range(n):
            outBeta,impBeta = py[outlet],th[imp]
            outP,impP = sps.expit(outBeta),sps.expit(impBeta)
            s,r=sens,spec
            z = outP + impP - outP*impP
            zTilde = s*z + (1-r)*(1-z)
            yDat,nSam = posMat[outlet,imp],numMat[outlet,imp]
            currElem = (1-outP)*(s+r-1)*(yDat/zTilde-(nSam-yDat)/(1-zTilde))*\
                        (impP - 3*(impP)**2 + 2*(impP)**3)+\
                        (((1-outP)*(impP-impP**2)*(s+r-1))**2)*\
                        (-yDat/zTilde**2-(nSam-yDat)/(1-zTilde)**2)
            currPartial += currElem
        impPartials[imp] = currPartial
    
    # outlet diagonals next
    outletPartials = np.zeros(n)
    for outlet in range(n):
        currPartial = 0
        for imp in range(m):
            outBeta,impBeta = py[outlet],th[imp]
            outP,impP = sps.expit(outBeta),sps.expit(impBeta)
            s,r=sens,spec
            z = outP + impP - outP*impP
            zTilde = s*z + (1-r)*(1-z)
            yDat,nSam = posMat[outlet,imp],numMat[outlet,imp]
            currElem = (1 - impP)*(yDat/zTilde-(nSam-yDat)/(1-zTilde))*\
                        (r+s-1)*(outP - 3*(outP**2) + 2*(outP**3)) +\
                        (1-impP)*(outP - outP**2 )*(s+r-1)*\
                        (yDat*((1-r-s)*(outP-outP**2)*(1-impP) )/(zTilde**2) -\
                        (nSam-yDat)*((s+r-1)*(outP-outP**2)*(1-impP))/(1-zTilde)**2)
            currPartial += currElem
        outletPartials[outlet] = currPartial
    
    diags = np.diag(np.concatenate((impPartials,outletPartials)))
    
    hess = (hess + diags)*-1 + 0.2 #* np.eye(hess.shape[0]) # include prior term
    return hess

def TRACKED_LogPrior(beta, numVec, posVec, sens, spec):
    '''
    Prior is a Laplace distribution with parameters: mu=-3, scale=4
    '''
    #-0.25*np.sum(np.abs(beta + 3)) - 0.001 * np.sum((beta + 3) ** 2)
    return - 0.1 * np.sum((beta + 3) ** 2)

def TRACKED_LogPrior_Grad(beta, nsamp, ydata, sens, spec):
    '''
    Prior is a Laplace distribution with parameters: mu=-3, scale=4
    '''
    #-0.25*np.squeeze(1*(beta >= -3) - 1*(beta <= -3)) - 0.002 * np.sum(beta + 3)
    return -0.2 * np.sum(beta + 3)

def TRACKED_LogPost(beta,N,Y,sens,spec):
    return TRACKED_LogPrior(beta,N,Y,sens,spec)\
           -TRACKED_NegLogLikeFunc(beta,N,Y,sens,spec,0)

def TRACKED_LogPost_Grad(beta, N, Y, sens, spec):
    return TRACKED_LogPrior_Grad(beta, N, Y, sens, spec)\
           -TRACKED_NegLogLikeFunc_Jac(beta,N,Y,sens,spec,0)
def GeneratePostSamps_TRACKED(N,Y,sens,spec,regWt,M,Madapt,delta):
    def TRACKEDtargetForNUTS(beta):
        return TRACKED_LogPost(beta,N,Y,sens,spec),\
               TRACKED_LogPost_Grad(beta,N,Y,sens,spec)

    beta0 = -2 * np.ones(N.shape[1] + N.shape[0])
    samples, lnprob, epsilon = nuts6(TRACKEDtargetForNUTS,M,Madapt,beta0,delta)
    
    return samples
###### END TRACKED FUNCTIONS ######







#### Necessary NUTS module ####
"""
This package implements the No-U-Turn Sampler (NUTS) algorithm 6 from the NUTS
paper (Hoffman & Gelman, 2011).

Content
-------

The package mainly contains:
  nuts6                     return samples using the NUTS
  test_nuts6                example usage of this package

and subroutines of nuts6:
  build_tree                the main recursion in NUTS
  find_reasonable_epsilon   Heuristic for choosing an initial value of epsilon
  leapfrog                  Perfom a leapfrog jump in the Hamiltonian space
  stop_criterion            Compute the stop condition in the main loop


A few words about NUTS
----------------------

Hamiltonian Monte Carlo or Hybrid Monte Carlo (HMC) is a Markov chain Monte
Carlo (MCMC) algorithm that avoids the random walk behavior and sensitivity to
correlated parameters, biggest weakness of many MCMC methods. Instead, it takes
a series of steps informed by first-order gradient information.

This feature allows it to converge much more quickly to high-dimensional target
distributions compared to simpler methods such as Metropolis, Gibbs sampling
(and derivatives).

However, HMC's performance is highly sensitive to two user-specified
parameters: a step size, and a desired number of steps.  In particular, if the
number of steps is too small then the algorithm will just exhibit random walk
behavior, whereas if it is too large it will waste computations.

Hoffman & Gelman introduced NUTS or the No-U-Turn Sampler, an extension to HMC
that eliminates the need to set a number of steps.  NUTS uses a recursive
algorithm to find likely candidate points that automatically stops when it
starts to double back and retrace its steps.  Empirically, NUTS perform at
least as effciently as and sometimes more effciently than a well tuned standard
HMC method, without requiring user intervention or costly tuning runs.

Moreover, Hoffman & Gelman derived a method for adapting the step size
parameter on the fly based on primal-dual averaging.  NUTS can thus be used
with no hand-tuning at all.

In practice, the implementation still requires a number of steps, a burning
period and a stepsize. However, the stepsize will be optimized during the
burning period, and the final values of all the user-defined values will be
revised by the algorithm.

reference: arXiv:1111.4246
"The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte
Carlo", Matthew D. Hoffman & Andrew Gelman
"""

from numpy import log, exp, sqrt

def leapfrog(theta, r, grad, epsilon, f):
    """ Perfom a leapfrog jump in the Hamiltonian space
    INPUTS
    ------
    theta: ndarray[float, ndim=1]
        initial parameter position

    r: ndarray[float, ndim=1]
        initial momentum

    grad: float
        initial gradient value

    epsilon: float
        step size

    f: callable
        it should return the log probability and gradient evaluated at theta
        logp, grad = f(theta)

    OUTPUTS
    -------
    thetaprime: ndarray[float, ndim=1]
        new parameter position
    rprime: ndarray[float, ndim=1]
        new momentum
    gradprime: float
        new gradient
    logpprime: float
        new lnp
    """
    # make half step in r
    rprime = r + 0.5 * epsilon * grad
    # make new step in theta
    thetaprime = theta + epsilon * rprime
    #compute new gradient
    logpprime, gradprime = f(thetaprime)
    # make half step in r again
    rprime = rprime + 0.5 * epsilon * gradprime
    return thetaprime, rprime, gradprime, logpprime



def find_reasonable_epsilon(theta0, grad0, logp0, f, epsilonLB = 0.005, epsilonUB = 0.5):
    """ Heuristic for choosing an initial value of epsilon """
    epsilon = (1)
    r0 = np.random.normal(0., 1., len(theta0))

    # Figure out what direction we should be moving epsilon.
    _, rprime, gradprime, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
    # brutal! This trick make sure the step is not huge leading to infinite
    # values of the likelihood. This could also help to make sure theta stays
    # within the prior domain (if any)
    k = 1.
    while np.isinf(logpprime) or np.isinf(gradprime).any():
        k *= 0.5
        _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon * k, f)

    epsilon = np.minimum(np.maximum(0.5 * k * epsilon, 2.*epsilonLB),epsilonUB/(2.))
    # acceptprob = np.exp(logpprime - logp0 - 0.5 * (np.dot(rprime, rprime.T) - np.dot(r0, r0.T)))
    # a = 2. * float((acceptprob > 0.5)) - 1.
    logacceptprob = logpprime-logp0-0.5*(np.dot(rprime, rprime)-np.dot(r0,r0))
    a = 1. if logacceptprob > np.log(0.5) else -1.
    # Keep moving epsilon in that direction until acceptprob crosses 0.5.
    # while ( (acceptprob ** a) > (2. ** (-a))):
    while a * logacceptprob > -a * np.log(2):
        epsilon = epsilon * (1.5 ** a)
        if epsilon < epsilonLB or epsilon > epsilonUB:
            break
        _, rprime, _, logpprime = leapfrog(theta0, r0, grad0, epsilon, f)
        # acceptprob = np.exp(logpprime - logp0 - 0.5 * ( np.dot(rprime, rprime.T) - np.dot(r0, r0.T)))
        logacceptprob = logpprime-logp0-0.5*(np.dot(rprime, rprime)-np.dot(r0,r0))

    #print("find_reasonable_epsilon=", epsilon) EOW commented out

    return epsilon


def stop_criterion(thetaminus, thetaplus, rminus, rplus):
    """ Compute the stop condition in the main loop
    dot(dtheta, rminus) >= 0 & dot(dtheta, rplus >= 0)

    INPUTS
    ------
    thetaminus, thetaplus: ndarray[float, ndim=1]
        under and above position
    rminus, rplus: ndarray[float, ndim=1]
        under and above momentum

    OUTPUTS
    -------
    criterion: bool
        return if the condition is valid
    """
    dtheta = thetaplus - thetaminus
    return (np.dot(dtheta, rminus.T) >= 0) & (np.dot(dtheta, rplus.T) >= 0)


def build_tree(theta, r, grad, logu, v, j, epsilon, f, joint0):
    """The main recursion."""
    if (j == 0):
        # Base case: Take a single leapfrog step in the direction v.
        thetaprime, rprime, gradprime, logpprime = leapfrog(theta, r, grad, v * epsilon, f)
        joint = logpprime - 0.5 * np.dot(rprime, rprime.T)
        # Is the new point in the slice?
        nprime = int(logu < joint)
        # Is the simulation wildly inaccurate?
        sprime = int((logu - 1000.) < joint)
        # Set the return values---minus=plus for all things here, since the
        # "tree" is of depth 0.
        thetaminus = thetaprime[:]
        thetaplus = thetaprime[:]
        rminus = rprime[:]
        rplus = rprime[:]
        gradminus = gradprime[:]
        gradplus = gradprime[:]
        # Compute the acceptance probability.
        alphaprime = min(1., np.exp(joint - joint0))
        #alphaprime = min(1., np.exp(logpprime - 0.5 * np.dot(rprime, rprime.T) - joint0))
        nalphaprime = 1
    else:
        # Recursion: Implicitly build the height j-1 left and right subtrees.
        thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime = build_tree(theta, r, grad, logu, v, j - 1, epsilon, f, joint0)
        # No need to keep going if the stopping criteria were met in the first subtree.
        if (sprime == 1):
            if (v == -1):
                thetaminus, rminus, gradminus, _, _, _, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2 = build_tree(thetaminus, rminus, gradminus, logu, v, j - 1, epsilon, f, joint0)
            else:
                _, _, _, thetaplus, rplus, gradplus, thetaprime2, gradprime2, logpprime2, nprime2, sprime2, alphaprime2, nalphaprime2 = build_tree(thetaplus, rplus, gradplus, logu, v, j - 1, epsilon, f, joint0)
            # Choose which subtree to propagate a sample up from.
            if (np.random.uniform() < (float(nprime2) / max(float(int(nprime) + int(nprime2)), 1.))):
                thetaprime = thetaprime2[:]
                gradprime = gradprime2[:]
                logpprime = logpprime2
            # Update the number of valid points.
            nprime = int(nprime) + int(nprime2)
            # Update the stopping criterion.
            sprime = int(sprime and sprime2 and stop_criterion(thetaminus, thetaplus, rminus, rplus))
            # Update the acceptance probability statistics.
            alphaprime = alphaprime + alphaprime2
            nalphaprime = nalphaprime + nalphaprime2

    return thetaminus, rminus, gradminus, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alphaprime, nalphaprime


def nuts6(f, M, Madapt, theta0, delta=0.25):
    """
    Implements the No-U-Turn Sampler (NUTS) algorithm 6 from from the NUTS
    paper (Hoffman & Gelman, 2011).

    Runs Madapt steps of burn-in, during which it adapts the step size
    parameter epsilon, then starts generating samples to return.

    Note the initial step size is tricky and not exactly the one from the
    initial paper.  In fact the initial step size could be given by the user in
    order to avoid potential problems

    INPUTS
    ------
    epsilon: float
        step size
        see nuts8 if you want to avoid tuning this parameter

    f: callable
        it should return the log probability and gradient evaluated at theta
        logp, grad = f(theta)

    M: int
        number of samples to generate.

    Madapt: int
        the number of steps of burn-in/how long to run the dual averaging
        algorithm to fit the step size epsilon.

    theta0: ndarray[float, ndim=1]
        initial guess of the parameters.

    KEYWORDS
    --------
    delta: float
        targeted acceptance fraction

    OUTPUTS
    -------
    samples: ndarray[float, ndim=2]
    M x D matrix of samples generated by NUTS.
    note: samples[0, :] = theta0
    """
        
    if len(np.shape(theta0)) > 1:
        raise ValueError('theta0 is expected to be a 1-D array')

    D = len(theta0)
    samples = np.empty((M + Madapt, D), dtype=float)
    lnprob = np.empty(M + Madapt, dtype=float)

    logp, grad = f(theta0)
    samples[0, :] = theta0
    lnprob[0] = logp

    # Choose a reasonable first epsilon by a simple heuristic.
    epsilon = find_reasonable_epsilon(theta0, grad, logp, f)

    # Parameters to the dual averaging algorithm.
    gamma = 0.05
    t0 = 10
    kappa = 0.75
    mu = log(10. * epsilon)

    # Initialize dual averaging algorithm.
    epsilonbar = 1
    Hbar = 0

    for m in range(1, M + Madapt):
        # Resample momenta.
        r0 = np.random.normal(0, 1, D)

        #joint lnp of theta and momentum r
        joint = logp - 0.5 * np.dot(r0, r0.T)

        # Resample u ~ uniform([0, exp(joint)]).
        # Equivalent to (log(u) - joint) ~ exponential(1).
        logu = float(joint - np.random.exponential(1, size=1))

        # if all fails, the next sample will be the previous one
        samples[m, :] = samples[m - 1, :]
        lnprob[m] = lnprob[m - 1]

        # initialize the tree
        thetaminus = samples[m - 1, :]
        thetaplus = samples[m - 1, :]
        rminus = r0[:]
        rplus = r0[:]
        gradminus = grad[:]
        gradplus = grad[:]

        j = 0  # initial heigth j = 0
        n = 1  # Initially the only valid point is the initial point.
        s = 1  # Main loop: will keep going until s == 0.
        
        while (s == 1):
            # Choose a direction. -1 = backwards, 1 = forwards.
            v = int(2 * (np.random.uniform() < 0.5) - 1)

            # Double the size of the tree.
            if (v == -1):
                thetaminus, rminus, gradminus, _, _, _, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha = build_tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, f, joint)
            else:
                _, _, _, thetaplus, rplus, gradplus, thetaprime, gradprime, logpprime, nprime, sprime, alpha, nalpha = build_tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, f, joint)

            # Use Metropolis-Hastings to decide whether or not to move to a
            # point from the half-tree we just generated.
            _tmp = min(1, float(nprime) / float(n))
            if (sprime == 1) and (np.random.uniform() < _tmp):
                samples[m, :] = thetaprime[:]
                lnprob[m] = logpprime
                logp = logpprime
                grad = gradprime[:]
            # Update number of valid points we've seen.
            n += nprime
            
            # Decide if it's time to stop.
            s = sprime and stop_criterion(thetaminus, thetaplus, rminus, rplus) and (n < 50) # (n<50) EOW EDIT
                
            # Increment depth.
            j += 1

        # Do adaptation of epsilon if we're still doing burn-in.
        eta = 1. / float(m + t0)
        Hbar = (1. - eta) * Hbar + eta * (delta - alpha / float(nalpha))
        if (m <= Madapt):
            epsilon = exp(mu - sqrt(m) / gamma * Hbar)
            epsilon = np.minimum(np.maximum(epsilon, 0.001),1)
            eta = m ** -kappa
            epsilonbar = exp((1. - eta) * log(epsilonbar) + eta * log(epsilon))
        else:
            epsilon = epsilonbar
                
    samples = samples[Madapt:, :]
    lnprob = lnprob[Madapt:]
    return samples, lnprob, epsilon