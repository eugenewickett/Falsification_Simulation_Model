# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 18:04:01 2020

@author: eugen
"""
import numpy as np
import matplotlib.pyplot as plt
import Falsification_Sim_Modules as simModules
import random
import csv
import cProfile
import time

'''
# Define the parameters of our smaller example
n = 15 # Number of outlets
m = 3 # Number of importers
# Generate a random transition matrix A
A = np.zeros(shape=(n,m))
for i in range(n):
    newRow = [np.floor(np.random.uniform()*100.0)]
    for j in range(m-2):
        newVal = np.floor(np.random.uniform()*(100.0-np.sum(newRow)))
        newRow.append(newVal)
    newRow.append(100.0-np.sum(newRow))
    newRow = [x / 100 for x in newRow]
    random.shuffle(newRow)
    A[i] = newRow
'''

A = np.array([[0.65, 0.09, 0.26],
               [0.24, 0.21, 0.55],
               [0.22, 0.13, 0.65],
               [0.26, 0.55, 0.19],
               [0.15, 0.52, 0.33],
               [0.48, 0.37, 0.15],
               [0.76, 0.08, 0.16],
               [0.02, 0.49, 0.49],
               [0.2 , 0.32, 0.48],
               [0.99, 0.01, 0.  ],
               [0.03, 0.68, 0.29],
               [0.55, 0.37, 0.08],
               [0.42, 0.5 , 0.08],
               [0.46, 0.17, 0.37],
               [0.  , 0.39, 0.61]])
pI = np.array((0.10, 0.30, 0.50))
pJ = np.array((0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0))
sens = 0.99
spec = 0.99
realproby = (1-pJ) * np.array((A @ pI)) + pJ #optimal testing
realprobz = realproby * sens + (1-realproby) * (1-spec) #real testing
nsamp = (500 * np.ones(A.shape[0])).astype('int')
ydata = np.squeeze(np.random.binomial(nsamp,realprobz))

'''
Other estimate methods:
X = np.array([ydata[i]/nsamp[i] for i in range(len(ydata))])
imps_Lin, outs_Lin = simModules.Est_LinearProjection(A,X)
imps_Bern, outs_Bern, d1, d2 = simModules.Est_BernMLEProjection(A,X)
imps_MLE, outs_MLE = simModules.PlumleeEstimates(ydata, nsamp, A, sens, spec, rglrWt = 0.1)

'''
# Questions:
#   effect of M,Madapt,delta,beta0 on time?
#   " on quality of samples (if there exists a difference)?
#   what different epsilons are being generated WRT Madapt?

beta0 = -2 * np.ones(A.shape[1]+A.shape[0])

M_Vec = [50]
Madapt_Vec = [50]
delta_Vec = [0.1]

def exampletargetfornuts(beta):
    """
    Example of a target distribution that could be sampled from using NUTS.
    (Although of course you could sample from it more efficiently)
    Doesn't include the normalizing constant.
    """
    return simModules.mylogpost(beta,ydata, nsamp, A, sens, spec), simModules.mylogpost_grad(beta,ydata, nsamp, A, sens, spec)

'''
fields=["Madapt", "M", "delta", "TOTAL TIME", "EPSILON"]
with open(r'NUTSAnalysis.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(fields)
'''

    #writer.writerow(["Madapt", "M", "delta", "TOTAL TIME", "EPSILON"]) # Initialize the header
for rep in range(10):
    for lklhdEst_Madapt in Madapt_Vec:
        for lklhdEst_M in M_Vec:
            for lklhdEst_delta in delta_Vec:
                startTime = time.time()
                samples, lnprob, epsilon = simModules.nuts6(exampletargetfornuts, lklhdEst_M, lklhdEst_Madapt, beta0, lklhdEst_delta)
                totalTime = time.time() - startTime
                reportRow = [lklhdEst_Madapt,lklhdEst_M,lklhdEst_delta,totalTime,epsilon]
                avgsVec = [np.mean(samples[:,x]) for x in range(len(samples[0]))]
                avgsVec = simModules.invlogit(avgsVec)
                reportRow = reportRow + avgsVec.tolist()
                with open(r'NUTSAnalysis.csv', 'a', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(reportRow)

''' 
with open('NUTSAnalysis.csv', newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

# Store run time,epsilon by delta
deltaArr_runTime = [[],[],[],[]]
deltaArr_eps = [[],[],[],[]]
deltaArr_impEsts = [[],[],[],[]]
for ind,d in delta_Vec:
    for dataRow in data:
        if dataRow[2]==d:
            deltaArr_runTime[ind].append(dataRow[3])
            deltaArr_eps[ind].append(dataRow[4])
            deltaArr_impEsts[ind].append(dataRow[5:8])



fig, axs = plt.subplots(3, 1,figsize=(17,13))
fig.suptitle('Run Time',fontsize=18)
for subP in range(3):  
    axs[subP].set_ylabel('Run Time (s)',fontsize=12)
    axs[subP].set_ylim([0,0.6])
# Delta
axs[0].set_title('Linear projection',fontweight='bold')        
axs[0].bar(End_Plot_x,estEndFalsePercList,color='peachpuff',edgecolor='red')
axs[0].set_xlabel('Delta',fontsize=12)
# Bernoulli MLE projection
axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
axs[1].bar(End_Plot_x,estEndFalsePercList_Bern,color='khaki',edgecolor='goldenrod')      
# MLE w Nonlinear optimizer        
axs[2].set_title('MLE w/ Nonlinear Optimizer',fontweight='bold')
axs[2].bar(End_Plot_x,estEndFalsePercList_Plum,color='lightcyan',edgecolor='teal')   
plt.tight_layout()
plt.subplots_adjust(top=0.94)
plt.show()
'''
    

''' 
currOutputLine = {'inputParameterDictionary':inputParameterDictionary,
                          'rootConsumption':List_RootConsumption,
                          'intDemandResults':List_demandResultsInt,
                          'endDemandResults':List_demandResultsEnd,
                          'testResults':TestReportTbl,
                          'intFalseEstimates':estIntFalsePercList,
                          'endFalseEstimates':estEndFalsePercList,
                          'intFalseEstimates_Bern':estIntFalsePercList_Bern,
                          'endFalseEstimates_Bern':estEndFalsePercList_Bern,
                          'intFalseEstimates_Plum':estIntFalsePercList_Plum,
                          'endFalseEstimates_Plum':estEndFalsePercList_Plum,
                          'falsePerc_LklhdSamples':estFalsePerc_LklhdSamples,
                          'intSFTrueValues':intSFVec,
                          'simStartTime':startTime,
                          'simRunTime':totalRunTime
                          }
                          
outputDict[rep] = currOutputLine

lklhdEst_M, lklhdEst_Madapt, lklhdEst_delta = 10, 10, 0.2

plt.hist(simModules.invlogit(samples[:,3]))
'''





