# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 16:35:26 2020

@author: eugen
"""

'''
These Quest dictionaries were run in batches of 50 replications each, with each
dictionary taking up ~30MBs of space.
'''

import os # for directories
import pickle # for saving/loading objects in Python
import matplotlib.pyplot as plt
import numpy as np

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

#import Falsification_Sim_Modules as simModules

directory = r'C:\Users\eugen\OneDrive\Documents\EAGER Project\Simulator'+\
         '\Sim Model Files TOO BIG FOR REPO\Checking Runs'
# BAD INTERMEDIATE NODE
OPFileNames = os.listdir(directory)
OPrepStarts = []
OPrepRTs = []
for item in OPFileNames:    
    fileName1 = directory + '\\'+item
    outputDict = pickle.load(open(fileName1, 'rb'))
    for repNum in outputDict.keys():
        currRepDict = outputDict[repNum]
        OPrepStarts.append(currRepDict['simStartTime'])
        OPrepRTs.append(currRepDict['simRunTime'])

OPrepRunTimes = [OPrepRTs[i][0] for i in range(len(OPrepRTs))]

OPrepStarts[0:10]
OPrepRunTimes[0:10]

plt.plot(OPrepStarts) 
plt.hist(OPrepRunTimes)
np.mean(OPrepRunTimes)
plt.plot(OPrepStarts,OPrepRunTimes,'o') 


600*50/3600 #max seconds times num reps



longReps = []
for item in OPFileNames:    
    fileName1 = directory + '\\'+item
    outputDict = pickle.load(open(fileName1, 'rb'))
    for repNum in outputDict.keys():
        if outputDict[repNum]['simRunTime'][0] > 500:
            longReps.append(outputDict[repNum]['inputParameterDictionary'])
            
itemList = []
for inputDict in longReps:
    itemList.append(inputDict['RglrWt']) 
            
plt.hist(itemList)
            
            
            
