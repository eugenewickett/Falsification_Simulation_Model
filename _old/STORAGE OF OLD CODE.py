# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 15:01:54 2020

@author: eugen
"""

'''
def myloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = ((1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ))
    probs_dirJ = -invlogit_grad(betaJ) * np.squeeze(np.array(A @ invlogit(betaI))) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0])).T
    ilK = np.array([invlogit_grad(betaI),]*(betaJ.shape[0]))
    probs_dirI =  np.multiply(np.multiply((1- ilJ ) , A) , np.squeeze(ilK)).T 
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  np.array(probz_dirI @ (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz))).T)
    negloglikJ =  np.squeeze(np.array((probz_dirJ * (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz)))).T))
    return -np.squeeze(np.concatenate([negloglikI,negloglikJ]))

def Est_LinearProjection(A,X): # Linear Projection
    # Uses the (estimated) transition matrix, A, and the (estimated) percentage SF
    # at each end node, X
    intProj = np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,X))
    endProj = np.subtract(X,np.dot(A,intProj))
    return np.ndarray.tolist(intProj.T)[0], np.ndarray.tolist(endProj.T)[0]

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

def PlumleeEstimates(ydata, numsamples, A, sens, spec, rglrWt = 0.1):
    ydata = np.array(ydata)
    numsamples = np.array(numsamples)
    beta0 = -6 * np.ones(A.shape[1]+A.shape[0])
    beta0 = beta0 + np.random.uniform(-1,1,np.size(beta0))
    def invlogit_INTERIOR(beta):
        return np.exp(beta)/(np.exp(beta)+1)
    def mynegloglik_INTERIOR(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit_INTERIOR(betaJ)) * np.matmul(A,invlogit_INTERIOR(betaI)) + invlogit_INTERIOR(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
            + rglrWt*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bds = spo.Bounds(beta0-8, beta0+8)
    opval = spo.minimize(mynegloglik_INTERIOR, beta0+1,
                         args=(ydata, numsamples, A, sens, spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bds)
    return invlogit_INTERIOR(opval.x)[0:A.shape[1]].tolist(), invlogit_INTERIOR(opval.x)[A.shape[1]:].tolist()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    '''