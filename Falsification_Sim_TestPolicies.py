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
    more exploiting at the end
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
    eps = np.exp(-1*(nextTestDay/int_totalDays)/arr_PolicyParameter[0])
    numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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









def dynamicTestingGenerator(resultsList,
                             int_totalDays = 1000,
                             int_numDaysRemain = 1000,
                             int_totalBudget = 1000,
                             int_sampleBudgetRemain = 1000,
                             int_PolicyType = 0, 
                             arr_PolicyParameter = [0]
                             ):
   
  
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    if int_PolicyType == 0: # Basic epsilon-greedy algorithm
        pass
    elif int_PolicyType==1:
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[nextTestDay] *arr_PolicyParameter[nextTestDay] *arr_PolicyParameter[nextTestDay]  # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
    
    elif int_PolicyType==2: #Exponential-decay epsilon
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = np.exp(-1*(nextTestDay/int_totalDays)/arr_PolicyParameter[0])
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
        ### END Exponential-decay epsilon policy ###
        
    elif int_PolicyType==3: #EPSILON FIRST POLICY
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[0] # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
            if nextTestDay > (1-eps)*int_totalBudget: # Exploit
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
        ### EPSILON FIRST POLICY
    
    elif int_PolicyType == 4: # THOMPSON-SAMPLING
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
        
        ### END THOMPSON SAMPLING ###
    elif int_PolicyType == 5: #Every-other sampling. Explore on even days, exploit on odd days
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = arr_PolicyParameter[0] # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
        
    elif int_PolicyType==6: #epsilon decreasing by r-curve (take square root of epsilon values)
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = np.sqrt(arr_PolicyParameter[nextTestDay]) # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
    elif int_PolicyType==7: #epsilon changing by sin
        nextTestDay = int_totalDays - int_numDaysRemain # The day we are generating a schedule for
        eps = (np.sin(12.4*nextTestDay)) # Our exploit parameter
        numToTest = int(np.floor(int_sampleBudgetRemain / int_numDaysRemain)) + min(int_sampleBudgetRemain % int_numDaysRemain,1) # How many samples to conduct in the next day
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
    
    elif int_PolicyType == 8: # TS w NUTS, version 1
        # Grab intermediate and end node distros, then project onto end nodes
        # for different samples from the distro. Pick the largest projected 
        # SF estimate
        # arr_PolicyParameter = [# days to plan for, sensitivity, specificity, M,
        #                        Madapt, delta]
        # How many days to plan for?
        numDaysToSched = min(arr_PolicyParameter[0],int_numDaysRemain)
        usedBudgetSoFar = 0
        firstTestDay = int_totalDays - int_numDaysRemain
        
        if int_numDaysRemain == int_totalDays: # Our initial schedule should just be a distrubed exploration
            currNode = resultsList[0][0]
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
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
            sens, spec, M, Madapt, delta = arr_PolicyParameter[1:]
            NUTSsamples = simModules.GenerateNUTSsamples(ydata,nSamp,A,sens,spec,M,Madapt,delta)
            # Now pick from these samples to generate projections
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
                for testInd in range(numToTest):    
                    currSample = simModules.invlogit(NUTSsamples[random.randrange(len(NUTSsamples))])
                    probs = currSample[A.shape[1]:] + np.matmul(A,currSample[:A.shape[1]])
                    # Normalize? Or just pick largest value
                    highInd = [i for i,j in enumerate(probs) if j == max(probs)]
                    currNode = resultsList[highInd[0]][0]
                    sampleSchedule.append([firstTestDay+currDay,currNode])
                    usedBudgetSoFar += 1
  
        
    elif int_PolicyType == 9: # Enhancing exploration w NUTS
        # Grab intermediate and end node distros. Identify intermediate node
        # sample variances. Pick an intermediate node according to a weighting.
        # Pick an outlet from this int node's column in the transition matrix
        # A, again by a weighting (where 0% nodes have a non-zero probability
        # of being selected). log((p/1-p) + eps)?
        # How many days to plan for?
        numDaysToSched = min(arr_PolicyParameter[0],int_numDaysRemain)
        usedBudgetSoFar = 0
        firstTestDay = int_totalDays - int_numDaysRemain
        
        if int_numDaysRemain == int_totalDays: # Our initial schedule should just be a distrubed exploration
            currNode = resultsList[0][0]
            for currDay in range(numDaysToSched):
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
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
            sens, spec, M, Madapt, delta = arr_PolicyParameter[1:]
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
                numToTest = int(np.floor((int_sampleBudgetRemain-usedBudgetSoFar) / (int_numDaysRemain-currDay))) +\
                            min((int_sampleBudgetRemain-usedBudgetSoFar) % (int_numDaysRemain-currDay),1) # How many samples to conduct in the next day
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
    
    else:
        print('Error generating the sampling schedule.')
    
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule