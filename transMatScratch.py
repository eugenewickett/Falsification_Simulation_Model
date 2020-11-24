# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 23:44:36 2020

@author: eugen
"""

import os # for directories
import pickle # for saving/loading objects in Python
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

import Falsification_Sim_Modules as simMods

directory = r'C:\Users\eugen\OneDrive\Documents\EAGER Project\Simulator'+\
         '\Sim Model Files TOO BIG FOR REPO\SIXNOV_RUNS\Consolidated Files'

OPFileNames = os.listdir(directory)

DFdata = []
for item in OPFileNames:
    fileName = directory + '\\'+item
    outputDict = pickle.load(open(fileName, 'rb'))
    for i in range(len(outputDict.keys())):
        currRep = outputDict[i]
        sampBudget = currRep['inputParameterDictionary']['samplingBudget']
        sens = currRep['inputParameterDictionary']['diagnosticSensitivity']
        spec = currRep['inputParameterDictionary']['diagnosticSpecificity']
        globDem = currRep['inputParameterDictionary']['globalDemandLevel']
        currQ = simMods.GenerateTransitionMatrix(currRep['dynTestResults'])
        detQtQ = np.linalg.det(currQ.transpose() @ currQ)
        trQtQ = np.trace(currQ.transpose() @ currQ)
        newLine = [sampBudget, sens,spec,globDem,detQtQ,trQtQ]
        DFdata.append(newLine)

headers = ['SamplingBudget','Sensitivity','Specificity','GlobalDemandLevel','Det(QtQ)','tr(QtQ)']
pdDF = pd.DataFrame(DFdata,columns=headers)
for h in headers[:4]:
    for j in headers[4:]:
        plt.figure(figsize=(13,7))
        plt.suptitle(h+' vs. '+j,fontsize=16)
        #plt.ylim(0.,1)
        ax = sns.boxplot(y=j,x=h,data=pdDF,palette='bright')
        ax.set_xlabel(h,fontsize=12)
        ax.set_ylabel(j,fontsize=12)
        #ax.set_xticklabels(xTickLabels,rotation='vertical',fontsize=14)
        #plt.setp(ax.get_legend().get_texts(), fontsize='12') # for legend text
        #plt.setp(ax.get_legend().get_title(), fontsize='14') # for legend title
        plt.show()
'''      
DFdata = [] # We will grow a list of tuples containing [dictionary,calc method, deviation]
for dictInd,currDict in enumerate(dictNamesVec):
    block1 = list(zip(itertools.cycle([currDict]),\
                      itertools.cycle(['Linear Projection']),\
                      truePos_Lin[dictInd]))
    #block2 = list(zip(itertools.cycle([currDict]),\
    #                  itertools.cycle(['Bernoulli Projection']),\
    #                  truePos_Bern[dictInd]))
    block3 = list(zip(itertools.cycle([currDict]),\
                      itertools.cycle(['Untracked MLE']),\
                      truePos_MLE[dictInd]))
    block4 = list(zip(itertools.cycle([currDict]),\
                      itertools.cycle(['Untracked Posterior Sample Means']),\
                      truePos_NUTS[dictInd]))
    block5 = list(zip(itertools.cycle([currDict]),\
                      itertools.cycle(['Tracked MLE']),\
                      truePos_MLEtr[dictInd]))
    for tup in block1:
        DFdata.append(tup)
    #for tup in block2:
    #    DFdata.append(tup)
    for tup in block3:
        DFdata.append(tup)
    for tup in block4:
        DFdata.append(tup)
    for tup in block5:
        DFdata.append(tup) 
    
AbsDevsDF = pd.DataFrame(DFdata,columns=headCol)
# Build boxplot
plt.figure(figsize=(13,7))
plt.suptitle('True Positive Rates: Threshold: '+r"$\bf{" + str(threshold) + "}$",fontsize=18)
plt.ylim(0.,1)
ax = sns.boxplot(y='devVal',x='dict',data=AbsDevsDF,palette='bright',\
                  hue='calcMethod')
ax.set_xlabel('Output Dictionary',fontsize=16)
ax.set_ylabel('True Positive Rate',fontsize=16)
ax.set_xticklabels(xTickLabels,rotation='vertical',fontsize=14)
plt.setp(ax.get_legend().get_texts(), fontsize='12') # for legend text
plt.setp(ax.get_legend().get_title(), fontsize='14') # for legend title
plt.show()
'''