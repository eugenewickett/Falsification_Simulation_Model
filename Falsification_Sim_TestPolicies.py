# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:48:57 2020

@author: eugen

This file contains possible static and dynamic testing policies for sampling
from end nodes. Static policies are called once at the beginning of the
simulation replication, while dynamic policies are called either every day
or on an interval basis. Each function takes the following inputs:
    1) resultsList: A list with rows corresponding to each end node, with each
                    row having the following format:[Node ID, Num Samples, 
                    Num Positive, Positive Rate, [IntNodeSourceCounts]]
    2) totalSimDays=1000: Total number of days in the simulation
    3) numDaysRemain=1000: Total number of days left in the simulation (same as
                           totalSimDays if a static policy)
    4) totalBudget=1000: Total sampling budget for the simulation run
    5) numBudgetRemain=1000: Total budget left, in number of samples (same as
                             totalBudget if a static policy)
    6) policyParamList=[0]: List of different policy parameters that might be
                            called by different policy functions
And outputs a single list, sampleSchedule, with the following elements in each entry:
    1) Day: Simulation day of the scheduled test
    2) Node: Which node to test on the respective day
"""

import numpy as np
import random
import Falsification_Sim_Modules as simModules

def testPolicyHandler(polType,resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    '''
    Takes in a testing policy choice, calls the respective function, and
    returns the generated testing schedule
    '''
    polStr = ['Static_Deterministic','Static_Random','Dyn_EpsGreedy',\
              'Dyn_EpsExpDecay','Dyn_EpsFirst','Dyn_ThompSamp','Dyn_EveryOther',\
              'Dyn_EpsSine','Dyn_TSwithNUTS','Dyn_ExploreWithNUTS']
    if polType not in polStr:
        raise ValueError("Invalid policy type. Expected one of: %s" % polStr)
    
    if polType == 'Static_Random':
        sampleSchedule = Pol_Stat_Deterministic(resultsList,totalSimDays,numDaysRemain,\
                                                totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Static_Random':
        sampleSchedule = Pol_Stat_Random(resultsList,totalSimDays,numDaysRemain,\
                                         totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_EpsGreedy':
        sampleSchedule = Pol_Dyn_EpsGreedy(resultsList,totalSimDays,numDaysRemain,\
                                           totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_EpsExpDecay':
        sampleSchedule = Pol_Dyn_EpsExpDecay(resultsList,totalSimDays,numDaysRemain,\
                                             totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_EpsFirst':
        sampleSchedule = Pol_Dyn_EpsFirst(resultsList,totalSimDays,numDaysRemain,\
                                          totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_ThompSamp':
        sampleSchedule = Pol_Dyn_ThompSamp(resultsList,totalSimDays,numDaysRemain,\
                                           totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_EveryOther':
        sampleSchedule = Pol_Dyn_EveryOther(resultsList,totalSimDays,numDaysRemain,\
                                            totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_EpsSine':
        sampleSchedule = Pol_Dyn_EpsSine(resultsList,totalSimDays,numDaysRemain,\
                                         totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_TSwithNUTS':
        sampleSchedule = Pol_Dyn_TSwithNUTS(resultsList,totalSimDays,numDaysRemain,\
                                            totalBudget,numBudgetRemain,policyParamList)
    elif polType == 'Dyn_ExploreWithNUTS':
        sampleSchedule = Pol_Dyn_ExploreWithNUTS(resultsList,totalSimDays,numDaysRemain,\
                                            totalBudget,numBudgetRemain,policyParamList)
    
    
    return sampleSchedule
    
    
    
    

def Pol_Stat_Deterministic(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Deterministic policy that rotates through each end node in numerical order
    until the sampling budget is exhausted, such that Day 1 features End Node 1,
    Day 2 features End Node 2, etc.
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    endNodes = []
    for nodeInd in range(len(resultsList)):
        endNodes.append(resultsList[nodeInd][0])
    # Generate a sampling schedule iterating through each end node
    nodeCount = 0
    currNode = endNodes[nodeCount]
    lastEndNode = endNodes[-1]        
    for samp in range(totalBudget):
        day = np.mod(samp,totalSimDays)
        sampleSchedule.append([day,currNode])
        if currNode == lastEndNode:
            nodeCount = 0
            currNode = endNodes[nodeCount]
        else:
            nodeCount += 1
            currNode = endNodes[nodeCount]
            
    sampleSchedule.sort(key=lambda x: x[0]) # Sort our schedule by day before output
    return sampleSchedule

def Pol_Stat_Random(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Random policy that selects random nodes on each day until the sampling
    budget is exhausted
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    endNodes = []
    for nodeInd in range(len(resultsList)):
        endNodes.append(resultsList[nodeInd][0])
    numEndNodes = len(endNodes)
    # Generate a sampling schedule randomly sampling the list of end nodes         
    for samp in range(totalBudget):
        day = np.mod(samp,totalSimDays)
        currEndInd = int(np.floor(np.random.uniform(low=0,high=numEndNodes,size=1)))
        currNode = endNodes[currEndInd]
        sampleSchedule.append([day,currNode])
        
    sampleSchedule.sort(key=lambda x: x[0]) # Sort our schedule by day before output
    return sampleSchedule

def Pol_Dyn_EpsGreedy(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Epsilon-greedy policy, where the first element of policyParamList is the
    desired exploration ratio, epsilon
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    eps = policyParamList[0] # Our explore parameter
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    # First grab the pool of highest SF rate nodes
    maxSFRate = 0
    maxIndsList = []
    for rw in resultsList:
        if rw[3] > maxSFRate:
            maxSFRate = rw[3]
    for currInd in range(len(resultsList)):
        if resultsList[currInd][3] == maxSFRate:
            maxIndsList.append(currInd)
    for testNum in range(numToTest):
        # Explore or exploit?
        if np.random.uniform() < 1-eps: # Exploit
            exploitBool = True
        else:
            exploitBool = False
        # Based on the previous dice roll, generate a sampling point
        if exploitBool:
            testInd = np.random.choice(maxIndsList)
            NodeToTest = resultsList[testInd][0]
        else:
            testInd = np.random.choice(len(resultsList))
            NodeToTest = resultsList[testInd][0]

        sampleSchedule.append([nextTestDay,NodeToTest])
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule
    
def Pol_Dyn_EpsExpDecay(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Similar to the epsilon-greedy strategy, except that the value of epsilon
    decays exponentially over time, resulting in more exploring at the start and
    more exploiting at the end; initial epsilon is drawn from the parameter list
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    eps = np.exp(-1*(nextTestDay/totalSimDays)/policyParamList[0])
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                    min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    # First grab the pool of highest SF rate nodes
    maxSFRate = 0
    maxIndsList = []
    for rw in resultsList:
        if rw[3] > maxSFRate:
            maxSFRate = rw[3]
    for currInd in range(len(resultsList)):
        if resultsList[currInd][3] == maxSFRate:
            maxIndsList.append(currInd)
    for testNum in range(numToTest):
        # Explore or exploit?
        if np.random.uniform() < 1-eps: # Exploit
            exploitBool = True
        else:
            exploitBool = False
        # Based on the previous dice roll, generate a sampling point
        if exploitBool:
            testInd = np.random.choice(maxIndsList)
            NodeToTest = resultsList[testInd][0]
        else:
            testInd = np.random.choice(len(resultsList))
            NodeToTest = resultsList[testInd][0]

        sampleSchedule.append([nextTestDay,NodeToTest])
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_EpsFirst(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Epsilon is now the fraction of our budget we devote to exploration before
    moving to pure exploitation
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    eps = policyParamList[0] # Our exploit parameter
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                    min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    # First grab the pool of highest SF rate nodes
    maxSFRate = 0
    maxIndsList = []
    for rw in resultsList:
        if rw[3] > maxSFRate:
            maxSFRate = rw[3]
    for currInd in range(len(resultsList)):
        if resultsList[currInd][3] == maxSFRate:
            maxIndsList.append(currInd)
    for testNum in range(numToTest):
        # Explore or exploit?
        if nextTestDay > (1-eps)*totalBudget: # Exploit
            exploitBool = True
        else:
            exploitBool = False
        # Based on the previous dice roll, generate a sampling point
        if exploitBool:
            testInd = np.random.choice(maxIndsList)
            NodeToTest = resultsList[testInd][0]
        else:
            testInd = np.random.choice(len(resultsList))
            NodeToTest = resultsList[testInd][0]

        sampleSchedule.append([nextTestDay,NodeToTest])
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_ThompSamp(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                      totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Thompson sampling, using the testing results achieved thus far
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                            min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    for testNum in range(numToTest):
        # Iterate through each end node, generating an RV according to the beta distribution of samples + positives
        betaSamples = []
        for rw in resultsList:
            alphaCurr = 1 + rw[2]
            betaCurr = 1 + (rw[1]-rw[2])
            sampleCurr = np.random.beta(alphaCurr,betaCurr)
            betaSamples.append(sampleCurr)
        # Select the highest variable
        maxSampleInd = betaSamples.index(max(betaSamples))
        NodeToTest = resultsList[maxSampleInd][0]
        sampleSchedule.append([nextTestDay,NodeToTest])
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_EveryOther(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                       totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Every-other sampling, where we exploit on even days, explore on odd days
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                    min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    # First grab the pool of highest SF rate nodes
    maxSFRate = 0
    maxIndsList = []
    for rw in resultsList:
        if rw[3] > maxSFRate:
            maxSFRate = rw[3]
    for currInd in range(len(resultsList)):
        if resultsList[currInd][3] == maxSFRate:
            maxIndsList.append(currInd)
    for testNum in range(numToTest):
        # Explore or exploit?
        if nextTestDay%2  == 1: # Exploit if we are on an odd sampling schedule day
            exploitBool = True
        else:
            exploitBool = False
        # Based on the previous dice roll, generate a sampling point
        if exploitBool:
            testInd = np.random.choice(maxIndsList)
            NodeToTest = resultsList[testInd][0]
        else:
            testInd = np.random.choice(len(resultsList))
            NodeToTest = resultsList[testInd][0]

        sampleSchedule.append([nextTestDay,NodeToTest])
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_EpsSine(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                       totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Epsilon follows a sine function of the number of days that have elapsed
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = totalSimDays - numDaysRemain # The day we are generating a schedule for
    eps = (np.sin(12.4*nextTestDay)) # Our exploit parameter
    numToTest = int(np.floor(numBudgetRemain / numDaysRemain)) +\
                    min(numBudgetRemain % numDaysRemain,1) # How many samples to conduct in the next day
    # Generate a sampling schedule using the current list of results
    # First grab the pool of highest SF rate nodes
    maxSFRate = 0
    maxIndsList = []
    for rw in resultsList:
        if rw[3] > maxSFRate:
            maxSFRate = rw[3]
    for currInd in range(len(resultsList)):
        if resultsList[currInd][3] == maxSFRate:
            maxIndsList.append(currInd)
    for testNum in range(numToTest):
        # Explore or exploit?
        if 0 < eps: # Exploit
            exploitBool = True
        else:
            exploitBool = False
        # Based on the previous dice roll, generate a sampling point
        if exploitBool:
            testInd = np.random.choice(maxIndsList)
            NodeToTest = resultsList[testInd][0]
        else:
            testInd = np.random.choice(len(resultsList))
            NodeToTest = resultsList[testInd][0]

        sampleSchedule.append([nextTestDay,NodeToTest])
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_TSwithNUTS(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                       totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Grab intermediate and end node distribtuions via NUTS, then project onto
    end nodes for different samples from the resulting distribution; pick
    the largest projected SF estimate
    policyParamList = [number days to plan for, sensitivity, specificity, M,
                            Madapt, delta]
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    # How many days to plan for?
    numDaysToSched = min(policyParamList[0],numDaysRemain)
    usedBudgetSoFar = 0
    firstTestDay = totalSimDays - numDaysRemain
    
    if numDaysRemain == totalSimDays: # Our initial schedule should just be a distrubed exploration
        currNode = resultsList[0][0]
        for currDay in range(numDaysToSched):
            numToTest = int(np.floor((numBudgetRemain-usedBudgetSoFar) / (numDaysRemain-currDay))) +\
                        min((numBudgetRemain-usedBudgetSoFar) % (numDaysRemain-currDay),1) # How many samples to conduct in the next day
            for testInd in range(numToTest): # Iterate through our end nodes
                if currNode > resultsList[len(resultsList)-1][0]:
                    currNode = resultsList[0][0]
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    currNode += 1
                else:                        
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    currNode += 1
                usedBudgetSoFar += 1
        
    else: # Generate NUTS sample using current results and use it to generate a new schedule
        ydata = []
        nSamp = []
        for rw in resultsList:
            ydata.append(rw[2])
            nSamp.append(rw[1])
        A = simModules.GenerateTransitionMatrix(resultsList)      
        sens, spec, M, Madapt, delta = policyParamList[1:]
        NUTSsamples = simModules.GenerateNUTSsamples(ydata,nSamp,A,sens,spec,M,Madapt,delta)
        # Now pick from these samples to generate projections
        for currDay in range(numDaysToSched):
            numToTest = int(np.floor((numBudgetRemain-usedBudgetSoFar) / (numDaysRemain-currDay))) +\
                        min((numBudgetRemain-usedBudgetSoFar) % (numDaysRemain-currDay),1) # How many samples to conduct in the next day
            for testInd in range(numToTest):    
                currSample = simModules.invlogit(NUTSsamples[random.randrange(len(NUTSsamples))])
                probs = currSample[A.shape[1]:] + np.matmul(A,currSample[:A.shape[1]])
                # Normalize? Or just pick largest value
                highInd = [i for i,j in enumerate(probs) if j == max(probs)]
                currNode = resultsList[highInd[0]][0]
                sampleSchedule.append([firstTestDay+currDay,currNode])
                usedBudgetSoFar += 1
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

def Pol_Dyn_ExploreWithNUTS(resultsList,totalSimDays=1000,numDaysRemain=1000,\
                            totalBudget=1000,numBudgetRemain=1000,policyParamList=[0]):
    """
    Grab intermediate and end node distribtuions via NUTS. Identify intermediate node
    sample variances. Pick an intermediate node, weighed towards picking those 
    with higher sample variances. Pick an outlet from this intermediate node's
    column in the transition matrix A, again by a weighting (where 0% nodes 
    have a non-zero probability of being selected). [log((p/1-p) + eps)?]
    policyParamList = [number days to plan for, sensitivity, specificity, M,
                       Madapt, delta]
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    # How many days to plan for?
    numDaysToSched = min(policyParamList[0],numDaysRemain)
    usedBudgetSoFar = 0
    firstTestDay = totalSimDays - numDaysRemain
    
    if numDaysRemain == totalSimDays: # Our initial schedule should just be a distrubed exploration
        currNode = resultsList[0][0]
        for currDay in range(numDaysToSched):
            numToTest = int(np.floor((numBudgetRemain-usedBudgetSoFar) / (numDaysRemain-currDay))) +\
                        min((numBudgetRemain-usedBudgetSoFar) % (numDaysRemain-currDay),1) # How many samples to conduct in the next day
            for testInd in range(numToTest): # Iterate through our end nodes
                if currNode > resultsList[len(resultsList)-1][0]:
                    currNode = resultsList[0][0]
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    currNode += 1
                else:                        
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    currNode += 1
                usedBudgetSoFar += 1
        
    else: # Generate NUTS sample using current results and use it to generate a new schedule
        ydata = []
        nSamp = []
        for rw in resultsList:
            ydata.append(rw[2])
            nSamp.append(rw[1])
        A = simModules.GenerateTransitionMatrix(resultsList)      
        sens, spec, M, Madapt, delta = policyParamList[1:]
        NUTSsamples = simModules.GenerateNUTSsamples(ydata,nSamp,A,sens,spec,M,Madapt,delta)
        # Store sample variances for intermediate nodes
        NUTSintVars = []
        for intNode in range(A.shape[1]):
            currVar = np.var(simModules.invlogit(NUTSsamples[:,intNode]))
            NUTSintVars.append(currVar)
        # Normalize sum of all variances to 1
        NUTSintVars = NUTSintVars/np.sum(NUTSintVars)
        
          
        # Now pick from these samples to generate projections
        for currDay in range(numDaysToSched):
            numToTest = int(np.floor((numBudgetRemain-usedBudgetSoFar) / (numDaysRemain-currDay))) +\
                        min((numBudgetRemain-usedBudgetSoFar) % (numDaysRemain-currDay),1) # How many samples to conduct in the next day
            for testInd in range(numToTest):    
                # Pick an intermediate node to "target", with more emphasis on higher sample variances
                rUnif = random.uniform(0,1)
                for intInd in range(A.shape[1]):
                    if rUnif < np.sum(NUTSintVars[0:(intInd+1)]):
                        targIntInd = intInd
                        break
                # Go through the same process with the column of A
                # pertaining to this target intermediate node
                AtargCol = [row[targIntInd] for row in A]
                # Add a small epsilon, for 0 values, and normalize
                AtargCol = np.add(AtargCol,1e-3)
                AtargCol = AtargCol/np.sum(AtargCol)
                rUnif = random.uniform(0,1)
                for intEnd in range(A.shape[0]):
                    if rUnif < np.sum(AtargCol[0:(intEnd+1)]):
                        currInd = intEnd
                        break
                currNode = resultsList[currInd][0]
                sampleSchedule.append([firstTestDay+currDay,currNode])
                usedBudgetSoFar += 1
        
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

