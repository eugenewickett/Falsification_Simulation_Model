# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 23:55:58 2020

@author: eugen
"""

'''
This script reads an CSV list of PMS testing results and returns a list of
estimation intervals for SFP rates for the Linear, Untracked, and Tracked
methods, as well as posterior distribution samples for the Untracked and
Tracked methods.
'''

import numpy as np
import csv
import DRA_EstimationMethods as estMethods

### PUT DATA FILE NAME HERE; IT MUST BE LOCATED IN THE SAME DIRECTORY AS THIS FILE
fileName = 'testResults.csv'
### ENTER DIAGNOSTIC SENSITIVITY AND SPECIFICITY
diagSens = 0.90
diagSpec = 0.99


dataTbl = [] #Initialize list for raw data
with open(fileName,newline='') as file:
    reader = csv.reader(file)
    for row in reader:
        dataTbl.append(row)

# Convert results into integers
for row in dataTbl:
    row[2] = int(row[2])

# Grab list of unique outlet and importer names
outletNames = []
importerNames = []
for row in dataTbl:
    if row[0] not in outletNames:
        outletNames.append(row[0])
    if row[1] not in importerNames:
        importerNames.append(row[1])
outletNames.sort()
importerNames.sort()

outletNum = len(outletNames)
importerNum = len(importerNames)

# Build N + Y matrices 
N = np.zeros(shape=(outletNum,importerNum))
Y = np.zeros(shape=(outletNum,importerNum))
for row in dataTbl:
    outInd = outletNames.index(row[0])
    impInd = importerNames.index(row[1])
    N[outInd,impInd] += 1
    Y[outInd,impInd] += row[2]

# Form Tracked estimates
outputTrackedDict = estMethods.Est_TrackedMLE(N,Y,diagSens,diagSpec)
# Form posterior samples
outputPostSamps = estMethods.GeneratePostSamps_TRACKED(N,Y,diagSens,diagSpec,\
                                                       regWt=0.,
                                                       M=1000,
                                                       Madapt=5000,
                                                       delta=0.4,
                                                       usePriors=1.)

# Generate output list of estimates and intervals
outputTrackedDict.keys()

# Generate output for posterior samples


# Generate plots of estimates and intervals



# Generate plots of posterior distributions


















