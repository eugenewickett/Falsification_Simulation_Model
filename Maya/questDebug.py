# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 16:49:53 2020

@author: eugen
"""
import numpy as np

def invlogit(beta):
    return np.array(np.exp(beta)/(np.exp(beta)+1))

def Est_BernMLEProjection(A,X): #MLE OF BERNOULLI VARIABLE
    # USING ITERATIVELY REWEIGHTED LEAST SQUARES, SEE WIKIPEDIA FOR NOTATION
    A = np.array(A)
    X = np.array(X)
    np.seterr(all='print')
    currGap = 10
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
        currGap = np.linalg.norm(w_k-w_k1)
        w_k = np.copy(w_k1)
    # Now our importer SF rates are calculated; figure out variance + Wald statistics
    covarMat_Bern = np.linalg.inv(np.dot(A.T,np.dot(S_k,A)))
    w_Var = np.diag(covarMat_Bern)
    wald_stats = []
    for j in range(m):
        wald_stats.append(float((w_k[j]**2)/w_Var[j]))
    
    # Convert to intermediate and end node estimates
    intProj = np.ndarray.tolist(invlogit(w_k.T.tolist()[0]))
    errs_Bern = np.subtract(X,mu_k)
    endProj = errs_Bern.T.tolist()[0]    
    return intProj, endProj, covarMat_Bern, wald_stats


Amat = np.array([[0.6,0.2,0.2],
                 [0.0,0.0,1.0],
                 [0.2,0.7,0.1],
                 [0.2,0.5,0.3],
                 [0.9,0.1,0.0],
                 [0.1,0.5,0.4]])
Xvec = np.array([0.2,0.1,0.4,0.3,0.1,0.5])





iProj, eProj, covmat, wStats = Est_BernMLEProjection(Amat,Xvec)










