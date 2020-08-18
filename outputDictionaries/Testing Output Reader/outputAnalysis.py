# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""

import numpy as np
import matplotlib.pyplot as plt
import os # for directories
import pickle # for saving/loading objects in Python
import Falsification_Sim_Modules as simModules

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\outputDictionaries\\Testing Output Reader'
os.chdir(dirStr) # Set directory    

strVec = ['1','2','3','4','5','6']

outputDicts = []
for strEnd in strVec:
    fileName = dirStr + '\\OPDICT_' + strEnd  
    currDict = pickle.load(open(fileName, 'rb'))
    outputDicts.append(currDict)

simModules.SimSFEstimateOutput(outputDicts)




#Compare different epsilon values
epsResultList = []
for epsStr in epsStrVec:
    
    if epsStr == 'TS':
        fileName = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Thompson Sampling\\TS'
        OPdict = pickle.load(open(fileName, 'rb'))
    elif epsStr[:2] == 'EG':
        fileName = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Epsilon-Greedy\\EPSGREEDY_'+epsStr[2:]
        OPdict = pickle.load(open(fileName, 'rb'))
    else:
        fileName = dirStr + '\\EXPDECAY_' + epsStr
        OPdict = pickle.load(open(fileName, 'rb'))
        
    rootConsumptionVec = [] # Store the root consumption data as a list
    for rep in OPdict.keys():
        currDictEntry = OPdict[rep]['rootConsumption']
        rootConsumptionVec.append(currDictEntry)
    
    intDemandVec = [] # Store the intermediate node demand data as a list
    for rep in OPdict.keys():
        currDictEntry = OPdict[rep]['intDemandResults']
        intDemandVec.append(currDictEntry)
        
    endDemandVec = [] # Store the end node demand data as a list
    for rep in OPdict.keys():
        currDictEntry = OPdict[rep]['endDemandResults']
        endDemandVec.append(currDictEntry)
        
    testResultVec = [] # Store the testing data as a list of lists
    for rep in OPdict.keys():
        currDictEntry = OPdict[rep]['testResults']
        testResultVec.append(currDictEntry)
        
    # Generate a vector of average root consumption percentages    
    avgFalseConsumedVec = []
    for item in rootConsumptionVec:
        avgFalseConsumedVec.append(item[1]/(item[0]+item[1]))
    avgFalseConsumedVec_unsort = avgFalseConsumedVec
    
    # Generate summaries of our test results 
    avgFalseTestedVec = []
    avgStockoutTestedVec = []
    avgGoodTestedVec = []
    for item in testResultVec:
        currNumFalse = 0
        currNumStockout = 0
        currNumGood = 0
        currTotal = len(item)
        for testResult in item:
            if testResult[2] == 1:
                currNumFalse += 1
            elif testResult[2] == -1:
                currNumStockout += 1
            elif testResult[2] == 0:
                currNumGood += 1
                
        avgFalseTestedVec.append(currNumFalse/currTotal)
        avgStockoutTestedVec.append(currNumStockout/currTotal)
        avgGoodTestedVec.append(currNumGood/currTotal)
    
    zip_object = zip(avgFalseTestedVec, avgFalseConsumedVec)
    testConsumeDiff = []
    for list1_i, list2_i in zip_object:
        testConsumeDiff.append(list1_i-list2_i)
    epsResultList.append(testConsumeDiff)


# Generate plots of differences for each epsilon used
lowErrInd = int(np.floor(0.05*len(OPdict.keys())))
upErrInd = int(np.ceil(0.95*len(OPdict.keys())))-1 
    
means = [np.mean(i) for i in epsResultList]
lowErrVec = []
upErrVec = []
for row in epsResultList:
    rowCopy = row.copy()
    rowCopy.sort()
    lowErrVec.append(rowCopy[lowErrInd])
    upErrVec.append(rowCopy[upErrInd])

error = [lowErrVec, upErrVec] 
fig = plt.figure()
ax = fig.add_axes([0,0,2,1])
ax.bar(epsStrVec, means,yerr=error,align='center',ecolor='black',
       capsize=5,color='bisque',edgecolor='darkorange')
ax.set_xlabel('Dynamic Sampling Policy',fontsize=16)
ax.set_ylabel('Difference from underlying SF rate',fontsize=16)
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
plt.show()








