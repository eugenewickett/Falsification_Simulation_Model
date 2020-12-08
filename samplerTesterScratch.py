# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 16:04:43 2020

@author: eugen
"""

import numpy as np
#import scipy.stats as spstat
import scipy.special as sps
import matplotlib.pyplot as plt
'''
Necessary NUTS sampler functions are included first, followed by gradient functions and
some small working examples
'''

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


'''
POSTERIORS AND OTHER FUNCTIONS
'''
def invlogit(beta):
    return sps.expit(beta)
    
def invlogit_grad(beta):
    return (np.exp(beta)/((np.exp(beta)+1) ** 2))

def UNTRACKED_Post(beta,N,Y,sens,spec,Q,regwt=0.):
    #prior distribution
    prior = - 0.1 * np.sum((beta + 3) ** 2)# -0.25*np.sum(np.abs(beta + 3))
    
    n,m = Q.shape
    th = beta[:m]
    py = beta[m:]
    betaInitial = -6*np.ones(m+n)
    pVec = invlogit(py)+(1-invlogit(py))*np.matmul(Q,invlogit(th))
    pVecTilde = sens*pVec + (1-spec)*(1-pVec)
    
    L = np.sum(np.multiply(Y,np.log(pVecTilde))+np.multiply(np.subtract(N,Y),\
               np.log(1-pVecTilde))) - regwt*np.sum(np.abs(py-betaInitial[m:]))
    
    return prior + L

def UNTRACKED_Post_Grad(beta,N,Y,sens,spec,Q,regwt=0.):
    #prior gradient
    prior_grad = -0.2 * np.sum(beta + 3)
    
    n,m = Q.shape
    th = beta[:m]
    py = beta[m:]
    betaInitial = -6*np.ones(m+n)
    pVec = invlogit(py)+(1-invlogit(py))*np.matmul(Q,invlogit(th))
    pVecTilde = sens*pVec + (1-spec)*(1-pVec)
    
    #Grab importers partials first, then outlets
    impPartials = np.sum(Y[:,None]*Q*invlogit_grad(th)*(sens+spec-1)*\
                     np.array([(1-invlogit(py))]*m).transpose()/pVecTilde[:,None]\
                     - (N-Y)[:,None]*Q*invlogit_grad(th)*(sens+spec-1)*\
                     np.array([(1-invlogit(py))]*m).transpose()/(1-pVecTilde)[:,None]\
                     ,axis=0)
    outletPartials = Y*(1-np.matmul(Q,invlogit(th)))*invlogit_grad(py)*\
                        (sens+spec-1)/pVecTilde - (N-Y)*invlogit_grad(py)*\
                        (sens+spec-1)*(1-np.matmul(Q,invlogit(th)))/(1-pVecTilde)\
                        - regwt*np.squeeze(1*(py >= betaInitial[m:]) - 1*(py <= betaInitial[m:]))
    
    like_jac = np.concatenate((impPartials,outletPartials))
    
    return prior_grad + like_jac
    
def UNTRACKED_Post_Hess(beta,N,Y,sens,spec,Q):
    n,m = Q.shape
    th = beta[:m]
    py = beta[m:]
    
    zVec = invlogit(py)+(1-invlogit(py))*np.matmul(Q,invlogit(th))
    zVecTilde = sens*zVec+(1-spec)*(1-zVec)
    sumVec = np.matmul(Q,invlogit(th))
    
    #initialize a Hessian matrix
    hess = np.zeros((n+m,n+m))
    # get off-diagonal entries first; importer-outlet entries
    for triRow in range(n):
        for triCol in range(m):
            outBeta,impBeta = py[triRow],th[triCol]
            outP,impP = invlogit(outBeta),invlogit(impBeta)
            s,r=sens,spec
            c1 = Q[triRow,triCol]*(s+r-1)*invlogit_grad(impBeta)
            yDat,nSam = Y[triRow],N[triRow]
            elem = c1*(1-outP)*(yDat*( (s+r-1)*(-sumVec[triRow]*(outP**2-outP) - outP + outP**2) )\
                    /( s*(sumVec[triRow]*(1 - outP) + outP) +\
                   (-r + 1)*(-sumVec[triRow]*(1 - outP) + 1 - outP) )**2 -\
                    (nSam - yDat)*((-r + 1-s)*(-sumVec[triRow]*(-outP + outP**2)-outP+outP**2))\
                     /(-s*(sumVec[triRow]*(1 - outP) + outP) - (1-r)*(-sumVec[triRow]*(1 - outP) +\
                   1 - outP) + 1)**2) +\
                    c1*(yDat/(s*(sumVec[triRow]*(1 - outP) + outP) + (-r + 1)*(-sumVec[triRow]*(1 - outP) +\
                   1 - outP)) - (nSam - yDat)/( -s*(sumVec[triRow]*(1 - outP) +\
                   outP) - (-r + 1)*(-sumVec[triRow]*(1 - outP) + 1 - outP) + 1))*( outP**2 - outP)
            hess[m+triRow,triCol] = elem
            hess[triCol,m+triRow] = elem
    # get off-diagonals for importer-importer entries
    for triCol in range(m-1):
        for triCol2 in range(triCol+1,m):
            elem = 0
            for i in range(n):
                nextPart = (sens+spec-1)*Q[i,triCol]*(1-invlogit(py[i]))*invlogit_grad(th[triCol])*\
                (-Y[i]*(sens+spec-1)*(1-invlogit(py[i]))*Q[i,triCol2]*(invlogit(th[triCol2]) - invlogit(th[triCol2])**2)            /\
                 (zVecTilde[i]**2)
                - (N[i]-Y[i])*(sens+spec-1)*(1-invlogit(py[i]))*Q[i,triCol2]*(invlogit(th[triCol2]) - invlogit(th[triCol2])**2) /\
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
            outP,impP = invlogit(outBeta),invlogit(impBeta)
            s,r=sens,spec                      
            c1 = Q[outlet,imp]*(s+r-1)*(1-outP)            
            c3 = (1-outP)*Q[outlet,imp]
            yDat,nSam = Y[outlet],N[outlet]
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
        outP = invlogit(outBeta)
        s,r=sens,spec
        c1 = Q[outlet] @ invlogit(th)
        c2 = (r + s - 1)
        yDat,nSam = Y[outlet],N[outlet]
        currPartial = (1-c1)*(yDat/(zVecTilde[outlet]) -\
                    (nSam - yDat)/(1-zVecTilde[outlet]))*c2*(outP -\
                    3*(outP**2) + 2*(outP**3)) + (1-c1)*(outP - outP**2 )*\
                    (yDat*(   -c2*(c1*(-outP + outP**2 ) + outP - outP**2 ) )/\
                    (zVecTilde[outlet])**2 - (nSam - yDat)*(c2*(c1*(-outP + outP**2) +\
                     outP - outP**2 ))/( -s*(c1*(1 - outP) +\
                     outP) - (1-r)*(1-c1*(1 - outP)  - outP) + 1 )**2)*c2
        outletPartials[outlet] = currPartial
    
    diags = np.diag(np.concatenate((impPartials,outletPartials)))
    
    hess = (hess + diags) - 0.002 * np.eye(hess.shape[0]) # include prior term   
    return hess

def TRACKED_Post(beta, N, Y, sens, spec, regwt=0.):
    #prior distribution
    prior = - 0.1 * np.sum((beta + 3) ** 2)# -0.25*np.sum(np.abs(beta + 3))
    
    #likelihood
    n,m = N.shape
    th = beta[:m]
    py = beta[m:]
    betaInitial = -6*np.ones(m+n)
    pMat = np.array([invlogit(th)]*n)+np.array([(1-invlogit(th))]*n)*\
            np.array([invlogit(py)]*m).transpose()
    pMatTilde = sens*pMat+(1-spec)*(1-pMat)
    L = np.sum(np.multiply(Y,np.log(pMatTilde))+np.multiply(np.subtract(N,Y),\
               np.log(1-pMatTilde))) - regwt*np.sum(np.abs(py-betaInitial[m:]))
    return prior + L

def TRACKED_Post_Grad(beta, N, Y, sens, spec, regwt=0.):
    #prior gradient
    prior_grad = -0.2 * np.sum(beta + 3)
    
    #likelihood Jacobian
    n,m = N.shape
    th = beta[:m]
    py = beta[m:]
    betaInitial = -6*np.ones(m+n)
    pMat = np.array([invlogit(th)]*n)+np.array([(1-invlogit(th))]*n)*\
            np.array([invlogit(py)]*m).transpose()
    pMatTilde = sens*pMat+(1-spec)*(1-pMat)
    
    #Grab importers partials first, then outlets
    impPartials = np.sum(Y*invlogit_grad(th)*(sens+spec-1)*\
                     np.array([(1-invlogit(py))]*m).transpose()/pMatTilde\
                     - (N-Y)*invlogit_grad(th)*(sens+spec-1)*\
                     np.array([(1-invlogit(py))]*m).transpose()/(1-pMatTilde)\
                     ,axis=0)
    outletPartials = np.sum((sens+spec-1)*(Y*invlogit_grad(py)[:,None]*\
                     np.array([(1-invlogit(th))]*n)/pMatTilde\
                     - (N-Y)*invlogit_grad(py)[:,None]*\
                     np.array([(1-invlogit(th))]*n)/(1-pMatTilde))\
                     ,axis=1) - regwt*np.squeeze(1*(py >= betaInitial[m:]) - 1*(py <= betaInitial[m:]))
       
    like_jac = np.concatenate((impPartials,outletPartials))
    
    return prior_grad+like_jac

def TRACKED_Post_Hess(betaVec,numMat,posMat,sens,spec):
    # betaVec should be [importers, outlets]
    n,m = numMat.shape
    th = betaVec[:m]
    py = betaVec[m:]
    
    zMat = np.array([invlogit(th)]*n)+np.array([(1-invlogit(th))]*n)*\
            np.array([invlogit(py)]*m).transpose()
    zMatTilde = sens*zMat+(1-spec)*(1-zMat)
    
    hess = np.zeros((n+m,n+m))
    # get off-diagonal entries first
    for triRow in range(n):
        for triCol in range(m):
            outBeta,impBeta = py[triRow],th[triCol]
            outP,impP = invlogit(outBeta),invlogit(impBeta)
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
            outP,impP = invlogit(outBeta),invlogit(impBeta)
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
            outP,impP = invlogit(outBeta),invlogit(impBeta)
            s,r=sens,spec
            z = outP + impP - outP*impP
            zTilde = s*z + (1-r)*(1-z)
            yDat,nSam = posMat[outlet,imp],numMat[outlet,imp]
            currElem = (1 - impP)*(yDat/zTilde-(nSam-yDat)/(1-zTilde))*\
                        (r+s-1)*(outP - 3*(outP**2) + 2*(outP**3)) +\
                        (1-impP)*(outP - outP**2 )*(s+r-1)*(yDat*\
                        ((1-r-s)*(outP-outP**2)*(1-impP) )/(zTilde**2) -\
                        (nSam-yDat)*((s+r-1)*(outP-outP**2)*(1-impP))/(1-zTilde)**2)
            currPartial += currElem
        outletPartials[outlet] = currPartial
    
    diags = np.diag(np.concatenate((impPartials,outletPartials)))
    
    hess = (hess + diags) - 0.002 * np.eye(hess.shape[0]) # include prior term
    return hess


'''
SAMPLER CALL FUNCTIONS
'''
def GeneratePostSamps_UNTRACKED(N,Y,Q,sens,spec,regWt,M,Madapt,delta):
    def UNTRACKEDtargetForNUTS(beta):
        return UNTRACKED_Post(beta,N,Y,sens,spec,Q,regWt),\
               UNTRACKED_Post_Grad(beta,N,Y,sens,spec,Q,regWt)

    beta0 = -4.5 * np.ones(Q.shape[1] + Q.shape[0])
    samples, lnprob, epsilon = nuts6(UNTRACKEDtargetForNUTS,M,Madapt,beta0,delta)
    
    return samples

def GeneratePostSamps_TRACKED(N,Y,sens,spec,regWt,M,Madapt,delta):
    def TRACKEDtargetForNUTS(beta):
        return TRACKED_Post(beta,N,Y,sens,spec,regWt),\
               TRACKED_Post_Grad(beta,N,Y,sens,spec,regWt)

    beta0 = -2 * np.ones(N.shape[1] + N.shape[0])
    samples, lnprob, epsilon = nuts6(TRACKEDtargetForNUTS,M,Madapt,beta0,delta)
    
    return samples








#SMALL EXAMPLES FOR TESTING OTHER SAMPLERS

#UNTRACKED
Sens,Spec,wt = 0.95,0.95,0.1
pImp_untr = np.array((0.01, 0.2, 0.1))
pOut_untr = np.array((0.01, 0.01, 0.2, 0.01, 0.1, 0.01))
Q = np.array([[0.5,0.2,0.3],
              [0.4,0.3,0.3],
              [0.1,0.6,0.3],
              [0.2,0.25,0.55],
              [0.15,0.55,0.3],
              [0.3,0.05,0.65]])
realproby_untr = (1-pOut_untr) * np.array((Q @ pImp_untr)) + pOut_untr #optimal testing
realprobz_untr = realproby_untr * Sens + (1-realproby_untr) * (1-Spec) #real testing
(n,m) = Q.shape
beta0Tr = -4.5*np.ones(n+m)
N_untr = (1000 * np.ones(n)).astype('int')
#np.random.seed(90) # Seeds 6, 8, 16, 18, 27 give negative Hessian diagonals
Y_untr = np.random.binomial(N_untr,realprobz_untr)

#TRACKED
pImp_tr = np.array((0.01, 0.2, 0.1))
pOut_tr = np.array((0.01, 0.01, 0.2, 0.01, 0.1, 0.01))
n = len(pOut_tr)
m = len(pImp_tr)
trackP = np.array([pImp_tr]*n)+np.array([1-pImp_tr]*n)*\
            np.array([pOut_tr]*m).transpose()
trackPtildes = Sens*trackP + (1-Spec)*(1-trackP)
beta0 = -4.5*np.ones(n+m)
numSamps = 300
N_tr = np.zeros(shape=(n,m))+numSamps
#np.random.seed(16) 
Y_tr = np.random.binomial(numSamps,trackPtildes)


M,Madapt,delta = 300,1000,0.4

#UNTRACKED samples
UNTR_samps = GeneratePostSamps_UNTRACKED(N_untr,Y_untr,Q,Sens,Spec,wt,M,Madapt,delta)

fig1=plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
n, bins, patches = ax1.hist(invlogit(UNTR_samps[:,:3])) #importers

fig2=plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
n, bins, patches = ax2.hist(invlogit(UNTR_samps[:,3:])) #outlets

# TRACKED samples
TR_samps = GeneratePostSamps_TRACKED(N_tr,Y_tr,Sens,Spec,wt,M,Madapt,delta)

fig3=plt.figure()
ax3 = fig3.add_subplot(1, 1, 1)
n, bins, patches = ax3.hist(invlogit(TR_samps[:,:3])) #importers

fig4=plt.figure()
ax4 = fig4.add_subplot(1, 1, 1)
n, bins, patches = ax4.hist(invlogit(TR_samps[:,3:])) #outlets




print(np.quantile(invlogit(TR_samps[:,]),(0.05,0.95) ,0).T)


