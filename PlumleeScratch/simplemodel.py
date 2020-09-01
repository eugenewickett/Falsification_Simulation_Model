# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 09:15:52 2020

@author: FloYd
"""
import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt
import numpy.random as npr

def PlumleeEstimates(ydata, numsamples, A, sens, spec):
    beta0 = -5 * np.ones(A.shape[1]+A.shape[0])
    
    def invlogit(beta):
        return np.exp(beta)/(np.exp(beta)+1)
    def logit(p):
        return np.log(p/(1-p)) 
    def mynegloglik(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
            + 4 * 1/2*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bounds = spo.Bounds(beta0-2, beta0+8)
    opval = spo.minimize(mynegloglik, beta0+5,
                         args=(ydata, numsamples, A, sens, spec),
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
sens = 0.99
spec = 0.99
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
numsamples = (100 * np.ones(A.shape[0])).astype('int')
ydata =np.random.binomial(numsamples,realprobz)

importerhat, outlethat = PlumleeEstimates(ydata, numsamples, A, sens, spec)
#print(importerhat)
#print(outlethat)


yd = np.array([ 0,  0,  1,  0,  0,  0,  1,  0,  0,  2,  0,  1,  1,  0,  0,  1,  0,
        0,  0,  0,  0,  1,  0,  0,  3,  0,  0,  2,  1,  0,  0,  1,  1,  0,
        0,  1,  1,  2,  0,  0,  0,  0,  0,  0,  0,  1,  5,  4,  3,  0,  0,
        0, 16,  0,  1,  2,  0,  1,  0,  2,  3,  0,  0,  1,  3,  1,  0,  2,
        0,  2,  1,  1,  0,  2,  0,  0,  0,  3,  1,  4,  0,  4,  0,  1,  1,
        0,  1,  1,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  3,  2,  1,
        4,  0,  0,  1])
    
nsamp = np.array([17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 16, 17, 17, 17, 17, 17,
       17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
       16, 16, 17, 17, 17, 17, 17, 15, 16, 17, 17, 15, 18, 18, 18, 18, 18,
       18, 18, 18, 18, 18, 18, 18, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17,
       17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 17, 17, 17,
       17, 17, 17, 17, 17, 17, 17, 16, 16, 17, 16, 17, 17, 17, 15, 17, 17,
       17, 17, 17, 17])

A = np.matrix([[0.        , 0.        , 0.13410753, 0.04653763, 0.        ,
        0.81935484, 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.12231559, 0.        , 0.        , 0.        ,
        0.07480029, 0.        , 0.        , 0.80288412, 0.        ],
       [0.        , 1.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.98690101, 0.        , 0.        ,
        0.        , 0.01309899, 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.96197938, 0.03802062,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 1.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 1.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.0740818 , 0.9259182 , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 1.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.52206835, 0.        , 0.20125536, 0.        ,
        0.        , 0.        , 0.0505131 , 0.2261632 , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.        , 1.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.13193963, 0.        , 0.        , 0.        ,
        0.86806037, 0.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.83773799, 0.        , 0.        , 0.        ,
        0.        , 0.01288245, 0.14937956, 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.00418479, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.99581521, 0.        ],
       [0.        , 0.        , 0.04306341, 0.42632779, 0.        ,
        0.        , 0.        , 0.        , 0.5306088 , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 1.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.25108295, 0.        ,
        0.74891705, 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 1.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 1.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.        , 0.        ],
       [0.        , 0.        , 0.        , 0.7792489 , 0.        ,
        0.2207511 , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.85242147,
        0.14757853, 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.39669049, 0.        ,
        0.        , 0.60330951, 0.        , 0.        , 0.        ],
       [0.        , 0.62679788, 0.        , 0.        , 0.37320212,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.62263703,
        0.        , 0.26694614, 0.        , 0.11041683, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        1.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.62795531, 0.        , 0.        , 0.        ,
        0.        , 0.37204469, 0.        , 0.        , 0.        ],
       [0.        , 0.62395642, 0.09112778, 0.        , 0.        ,
        0.        , 0.        , 0.28491581, 0.        , 0.        ],
       [0.        , 0.09407402, 0.14357073, 0.        , 0.        ,
        0.12737456, 0.        , 0.34049799, 0.        , 0.29448271],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.78904927, 0.        , 0.        , 0.        , 0.        ,
        0.21095073, 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.78808446, 0.21191554, 0.        ],
       [0.        , 0.        , 0.        , 0.9587776 , 0.        ,
        0.        , 0.        , 0.        , 0.0412224 , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.00709284, 0.        , 0.99290716],
       [0.26276553, 0.        , 0.        , 0.        , 0.        ,
        0.09955912, 0.        , 0.63767535, 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.0103211 ,
        0.7761115 , 0.        , 0.        , 0.        , 0.2135674 ],
       [0.        , 0.08205917, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.91794083, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.84865666,
        0.15134334, 0.        , 0.        , 0.        , 0.        ],
       [0.66798956, 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.33201044, 0.        ],
       [0.        , 0.        , 0.06583428, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.93416572, 0.        ],
       [0.19248957, 0.        , 0.        , 0.29235049, 0.        ,
        0.        , 0.        , 0.14770515, 0.        , 0.3674548 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.98781952, 0.        , 0.        , 0.01218048, 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.79107098, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.20892902, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.93330292, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.06669708],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.00517534, 0.92993975, 0.06488491, 0.        ],
       [0.        , 0.65084893, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.34915107, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.04949966, 0.03468031, 0.91582003, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.26684476, 0.        , 0.56567231, 0.        ,
        0.        , 0.00304245, 0.05246364, 0.        , 0.11197685],
       [0.82492038, 0.        , 0.17507962, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.88720879, 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.11279121, 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.47033024,
        0.        , 0.        , 0.12633299, 0.40333677, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.79872965, 0.        , 0.20127035, 0.        , 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 1.        , 0.        , 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.74472256, 0.13336852,
        0.        , 0.12190893, 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        1.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 1.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.02960383, 0.        , 0.59982586,
        0.        , 0.        , 0.37057031, 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.99319667, 0.        , 0.00680333, 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.88326046, 0.03324808, 0.08349146, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 1.        ],
       [1.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 1.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 1.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 1.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ]])

sens = 0.99
spec = 0.99

### CHANGE INITIAL EPSILON FOR A
epsVec = [0.001,1e-4,1e-5,1e-6]
importerHatVec = []
for e in epsVec:
    A = A + e
    importerhat, outlethat = PlumleeEstimates(yd, nsamp, A, sens, spec)
    importerHatVec.append(importerhat.tolist())
#print(importerhat)
#print(outlethat)

importerHatVec

fig = plt.figure()
ax = fig.add_axes([0,0,1,0.5])
#ax.bar(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[0],
#       align='center',ecolor='black',
#       capsize=5,color='lightcoral',edgecolor='firebrick')
ax.set_xlabel('Intermediate Node',fontsize=16)
ax.set_ylabel('Est. falsification %',fontsize=16)
line1, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[0])
line2, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[1])
line3, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[2])
line4, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[3])
line1.set_label(epsVec[0])
line2.set_label(epsVec[1])
line3.set_label(epsVec[2])
line4.set_label(epsVec[3])
ax.legend()


### CHANGE INITIAL BETA

betaVec = [-1,0,1,2,3,4,5]
importerHatVec = []
for b in betaVec:
    A = A + 0.001
    def PlumleeEstimates(ydata, numsamples, A, sens, spec):
        beta0 = -5 * np.ones(A.shape[1]+A.shape[0])
        
        def invlogit(beta):
            return np.exp(beta)/(np.exp(beta)+1)
        def logit(p):
            return np.log(p/(1-p)) 
        def mynegloglik(beta, ydata, numsamples, A, sens, spec):
            betaI = beta[0:A.shape[1]]
            betaJ = beta[A.shape[1]:]
            probs = (1-invlogit(betaJ)) * np.array(A @ invlogit(betaI)) + invlogit(betaJ)
            probsz = probs*sens + (1-probs) * (1-spec)
            return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
                + 4 * 1/4*np.sum(np.abs((betaJ - beta[A.shape[1]:]))) #have to regularize to prevent problems
        
        bounds = spo.Bounds(beta0-2, beta0+8)
        opval = spo.minimize(mynegloglik, beta0,
                             args=(ydata, numsamples, A, sens, spec),
                             method='L-BFGS-B',
                             options={'disp': False},
                             bounds=bounds)
        return invlogit(opval.x)[0:A.shape[1]], invlogit(opval.x)[A.shape[1]:]
    
    importerhat, outlethat = PlumleeEstimates(yd, nsamp, A, sens, spec)
    importerHatVec.append(importerhat.tolist())

importerHatVec

fig = plt.figure()
ax = fig.add_axes([0,0,1,0.5])
#ax.bar(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[0],
#       align='center',ecolor='black',
#       capsize=5,color='lightcoral',edgecolor='firebrick')
ax.set_xlabel('Intermediate Node',fontsize=16)
ax.set_ylabel('Est. falsification %',fontsize=16)
line1, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[0])
line2, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[1])
line3, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[2])
line4, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[3])
line5, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[4])
line6, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[5])
line7, = ax.plot(['2', '3', '4', '5', '6', '7', '8', '9', '10', '11'], importerHatVec[6])
line1.set_label('beta add term: ' + str(betaVec[0]))
line2.set_label('beta add term: ' + str(betaVec[1]))
line3.set_label('beta add term: ' + str(betaVec[2]))
line4.set_label('beta add term: ' + str(betaVec[3]))
line5.set_label('beta add term: ' + str(betaVec[4]))
line6.set_label('beta add term: ' + str(betaVec[5]))
line7.set_label('beta add term: ' + str(betaVec[6]))
ax.legend()







#test case
A = np.matrix([[0.5, 0.25, 0.25], 
    [0.001, 0.499, 0.5],
   [0.75, 0.001, 0.249],
    [0.499, 0.001, 0.5],
    [0.001, 0.249, 0.75],
    [0.001, 0.001, 0.998]])
pI = np.array((0.001, 0.2, 0.001))
pJ = np.array((0.001, 0.001, 0.2, 0.001, 0.001, 0.001))
sens = 0.99
spec = 0.99
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
numsamples = (100 * np.ones(A.shape[0])).astype('int')
ydata =np.random.binomial(numsamples,realprobz)


# METROPOLIS-HASTINGS
p = 0.10 # initialize w presumed SF rate
lapScale = 1 # initialize w 1; LOOK TO MODIFY
jumpDistStDev = 0.1 # standard deviation for the jumping distribution; LOOK TO MODIFY
x0 = npr.laplace(logit(p),lapScale,A.shape[0]+A.shape[1]) 
xinv = invlogit(x0)
tol = 0.02

ptArr = []
ptArr.append(x0)
numGens = 10000 # how many points to generate?
numRej = 0 # counter for rejections
for i in range(numGens):
    xCurr = ptArr[-1] #Last submitted point
    xProp = xCurr + npr.normal(loc=0.0,scale=jumpDistStDev,size=len(xCurr)) #Proposed point
    aRatio = mynegloglik(xProp, ydata, numsamples, A, sens, spec) \
         - mynegloglik(xCurr, ydata, numsamples, A, sens, spec)
    if npr.uniform(size=1) < np.exp(aRatio)-tol:
        ptArr.append(xProp)
        xCurr = np.copy(xProp)
    else:
        ptArr.append(np.copy(xCurr))
        numRej += 1
  
#plot of resulting distribution

testVec = []
for rw in ptArr:
    testVec.append(rw[4])
plt.plot(invlogit(testVec))
plt.hist(invlogit(testVec))
ptArr[99999]


























