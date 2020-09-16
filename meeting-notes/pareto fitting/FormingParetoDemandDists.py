# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 17:58:28 2020

@author: eugen
"""

'''
Aim to modify preference lists of outlets so that a Pareto-like
demand is generated WRT importers
'''

# Iterate through each outlet in order, updating preferences based
# on the prior outlets rankings
import numpy as np
import random
import matplotlib.pyplot as plt
import csv

paretoConstant = 1.1 # larger values cause larger density towards a few importers
prefDensity = 0.5 # density of preference lists relative to number of importers
endNum = 106
impNum = 10

random.seed(987)
np.random.seed(123)
newPrefLists = [random.sample(range(impNum),np.random.binomial(impNum,prefDensity))]
rankSumSoFar = np.zeros(impNum)+4.5 # to initialize
rankCountSoFar = np.zeros(impNum)+1
for i in range(len(newPrefLists[0])):
    currImp = newPrefLists[0][i]
    rankSumSoFar[currImp] += i
    rankCountSoFar[currImp] += 1
    
rankAvgSoFar = [rankSumSoFar[j]/rankCountSoFar[j] for j in range(impNum)]
# First pick new set of importers, then assign order based on prior average rankings
for ind in range(endNum-1):
    newImps = random.sample(range(impNum),np.random.binomial(impNum,prefDensity))
    newImps.sort()
    currNewPrefList = []
    while len(newImps)>1:
        probsList = [1/(rankAvgSoFar[newImps[i]]**paretoConstant) for i in range(len(newImps))]
        probsListSum = np.sum(probsList)
        probsList = [np.sum(probsList[0:i+1])/probsListSum for i in range(len(newImps))]
        rUnif = random.uniform(0,1)
        for j in range(len(newImps)):
            if probsList[j] > rUnif:
                currNewPrefList.append(newImps[j])
                newImps.remove(newImps[j])
                break
    currNewPrefList.append(newImps[0]) # add our last importer
    for i in range(len(currNewPrefList)): # update rankings
        currImp = currNewPrefList[i]
        rankSumSoFar[currImp] += i
        rankCountSoFar[currImp] += 1
    rankAvgSoFar = [rankSumSoFar[j]/rankCountSoFar[j] for j in range(impNum)]
    newPrefLists.append(currNewPrefList)

# make plot of first choices
firstPrefCount = np.zeros(impNum)
for l in newPrefLists:
    firstPref = l[0]
    firstPrefCount[firstPref] += 1

firstPrefCount = firstPrefCount.tolist()
firstPrefCount.sort()
# plot
plt.bar(np.arange(0,10).tolist(),firstPrefCount)

# make plot of second choices
secondPrefCount = np.zeros(impNum)
for l in newPrefLists:
    if len(l)>1:
        secondPref = l[1]
        secondPrefCount[secondPref] += 1

secondPrefCount = secondPrefCount.tolist()
secondPrefCount.sort()
# plot
plt.bar(np.arange(0,10).tolist(),secondPrefCount)

# make plot of third choices
thirdPrefCount = np.zeros(impNum)
for l in newPrefLists:
    if len(l)>2:
        thirdPref = l[2]
        thirdPrefCount[thirdPref] += 1

thirdPrefCount = thirdPrefCount.tolist()
thirdPrefCount.sort()
# plot
plt.bar(np.arange(0,10).tolist(),thirdPrefCount)



### Putting the new pref lists into a usable form for the simulation
# Add 2 to each importer number, to match the ID
# Add '1' at the end, to denote the SFP node
for ind,l in enumerate(newPrefLists):
    new_l = [l[i]+2 for i in range(len(l))]
    new_l = new_l + [1]
    newPrefLists[ind] = new_l

# Now shuffle the lists, as there are other node characteristics that are ordered numerically
random.shuffle(newPrefLists)

consumpImpVec = [IntermediateReportTbl[i][1] for i in range(impNum)]
consumpImpVec.sort() 
plt.bar(np.arange(0,impNum),consumpImpVec)


# Print to a csv file
csvTbl = []
for l in newPrefLists:
    newRow = []
    for imp in range(1,12):
        if imp in l:
            newRow.append(l.index(imp)+1)
        else:
            newRow.append(0)
    csvTbl.append(newRow)

with open('paretoOutput.csv','w',newline='') as f:
    writer = csv.writer(f)
    writer.writerows(csvTbl)

    

  
    
