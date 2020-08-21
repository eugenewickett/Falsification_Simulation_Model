# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:04:36 2019

@author: Eugene Wickett

Stores modules for use with 'SC Simulator.py'
"""
import numpy as np
import scipy.optimize as spo
import scipy.special as sps
import csv
import os
import sys
import pickle
import matplotlib.pyplot as plt
import winsound
import Falsification_Sim_Classes as simClasses # modules for the simulation

def CEOTTKbeep():
    winsound.Beep(int(32.7032 * (2**3)*(1.059463094**10)),400)
    winsound.Beep(int(32.7032 * (2**3)*(1.059463094**12)),400)
    winsound.Beep(int(32.7032 * (2**3)*(1.059463094**8)),400)
    winsound.Beep(int(32.7032 * (2**2)*(1.059463094**8)),400)
    winsound.Beep(int(32.7032 * (2**3)*(1.059463094**3)),400)  

def generateNodeListsFromFile(nodeInputFileString='',
                              arcPreferencesFileString='',
                              arcLTsFileString='',
                              arcRsFileString='',
                              NumSimDays=0,     
                              ):
    """
    Processes node and arc files and returns node lists 
    """
    List_RootNode = []
    List_IntermediateNode = []
    List_EndNode = []
    # Read the arc and node files as matrices/lists
    with open(nodeInputFileString) as nodeFile, open(arcPreferencesFileString) as arcPrefsFile, open(arcLTsFileString) as arcLTsFile, open(arcRsFileString) as arcRsFile:
        nodeReader = csv.reader(nodeFile,delimiter=',')
        nodeList = list(nodeReader)
        arcPreferencesReader = csv.reader(arcPrefsFile,delimiter=',')
        arcPrefsList = list(arcPreferencesReader)
        arcPreferencesMatrix = np.array(arcPrefsList).astype("float") #Preferences in matrix form
        arcLTsReader = csv.reader(arcLTsFile,delimiter=',')
        arcLTsList = list(arcLTsReader)
        arcLTsMatrix = np.array(arcLTsList).astype("float") #LTs in matrix form
        arcRsReader = csv.reader(arcRsFile,delimiter=',')
        arcRsList = list(arcRsReader)
        arcRsMatrix = np.array(arcRsList).astype("float") #LTs in matrix form
        
        
    # Generate objects, lists, etc. for the nodes
    nodeListHeader = nodeList[0] # Remove the headers from the list
    nodeList = nodeList[1:] # Retain the node rows only
    nodeNum = len(nodeList)
    for currNode in range((len(nodeList))):
        currRootTF = nodeList[currNode][0]            # Root node?
        currEndTF = nodeList[currNode][1]             # End node?
        if currRootTF == 'FALSE':
            currReorderPoint = int(nodeList[currNode][2])
            currReorderAmount = int(nodeList[currNode][3])
        if currEndTF == 'TRUE': # We have an end node
            currDemandDistribution = int(nodeList[currNode][4])
            currDemandParamList = []
            if not nodeList[currNode][5] == '':
                currDemandParam1 = float(nodeList[currNode][5])
                currDemandParamList.append(currDemandParam1)
            if not nodeList[currNode][6] == '':
                currDemandParam2 = float(nodeList[currNode][6])
                currDemandParamList.append(currDemandParam2)
            if not nodeList[currNode][7] == '':
                currDemandParam3 = float(nodeList[currNode][7])
                currDemandParamList.append(currDemandParam3)
        
        # Add root node to List_RootNode; otherwise generate the preferred supplier list for this non-root node
        if currRootTF == 'TRUE': # Currently a root node
            List_RootNode.append(currNode) # List_Root only has integers
        
        else: # Generate preferred supplier list
            currPreferredSuppliers = [] # Initialize our preferred supplier list
            currPreferredSuppliersLTs = []
            currPreferredSuppliersRs = []
            numSuppliers = int(max(arcPreferencesMatrix.T[currNode])) # How many suppliers are ranked
            for indSupplier in range(numSuppliers): # Add these suppliers by rank
                findNum = indSupplier+1
                for currSupplier in range(nodeNum): # Find the next ranked supplier
                    if arcPreferencesMatrix.T[currNode][currSupplier] == findNum:
                        currPreferredSuppliers.append(currSupplier) #Add the supplier's node number
                        currPreferredSuppliersLTs.append(int(arcLTsMatrix[currSupplier][currNode])) # Add the lead time from the supplier to the current Node
                        currPreferredSuppliersRs.append(int(arcRsMatrix[currSupplier][currNode])) # Add the lead time from the supplier to the current Node
            # Supplier lists are now set            
            
            # Generate Node objects, including demand schedules if the current node is a root node
            if currEndTF == 'FALSE': # Currently an intermediate node
                # Pull the falsification procurement probability
                currFalsifierProbability = float(nodeList[currNode][9])
                # Add a node object to the List_IntermediateNode
                newIntermediateNode = simClasses.Node(nodeNum=currNode, reorderPoint=currReorderPoint, reorderAmount=currReorderAmount, preferredSupplierVec=currPreferredSuppliers, preferredSupplierLTsVec=currPreferredSuppliersLTs, preferredSupplierRsVec=currPreferredSuppliersRs, endNodeTF = False, falseProb = currFalsifierProbability)
                List_IntermediateNode.append(newIntermediateNode)
            else: # Currently an end node
                # Add a node object to the List_End
                newDemandSchedule = demandScheduleGenerator_EndNode(int_numDays=NumSimDays, int_DistType=currDemandDistribution, arr_DistParameter=currDemandParamList)
                newEndNode = simClasses.Node(nodeNum=currNode, reorderPoint=currReorderPoint, reorderAmount=currReorderAmount, preferredSupplierVec=currPreferredSuppliers, preferredSupplierLTsVec=currPreferredSuppliersLTs, preferredSupplierRsVec=currPreferredSuppliersRs, endNodeTF=True,demandSched=newDemandSchedule)
                List_EndNode.append(newEndNode)    
    ### END NODE FOR LOOP ###
    return List_RootNode, List_IntermediateNode, List_EndNode, nodeListHeader, nodeList, nodeNum, arcPreferencesMatrix, arcLTsMatrix, arcRsMatrix


def demandScheduleGenerator_EndNode(int_numDays=1000,
                             int_DistType = 0, 
                             arr_DistParameter = [4]
                             ):
    """
    Randomly generates a demand schedule for a desired number of time units, 
    with the following parameters entered:
            1) int_numDays: Number of batches desired
            2) int_DistType: Distribution for the demand of each time period
                0 = Constant; parameter array should be [Constant value]
                1 = Uniform; parameter array should be [Low Bound, High Bound]
                2 = Triangular; parameter array should be [Left,Mode,Right]
                3 = Poisson; parameter array should be [Mean]
                ...
            3) arr_DistParameter: Array of the parameters required for generatng
                the demands
                0 = Constant; parameter array should be [Constant value]
                1 = Uniform; paramter array should be [Low Bound, High Bound]
                2 = Triangular; parameter array should be [Left,Mode,Right]
                3 = Poisson; parameter array should be [Mean]
                ...
    
    Outputs a Python list with the following elements within each entry:
            1) Demand: Integer denoting the demand for each simulation day
    """
    #Initialize our output, a list with the above mentioned outputs
    demandSchedule = []
    
    if int_DistType == 0: #Constant
        #Use the first value of the parameter array
        demandSchedule = np.repeat(arr_DistParameter[0],int_numDays)
    elif int_DistType == 1: #Uniform 
        #Use the first two values of the parameter array
        demandSchedule = np.round(np.random.uniform(low=arr_DistParameter[0],high=arr_DistParameter[1],size=int_numDays))      
    elif int_DistType == 2: #Triangular
        #Use the first two values of the parameter array
        demandSchedule = np.round(np.random.triangular(left=arr_DistParameter[0],mode=arr_DistParameter[1],right=arr_DistParameter[2],size=int_numDays))
    elif int_DistType == 3: #Poisson
        #Use the first value of the parameter array
        demandSchedule = np.round(np.random.poisson(lam=arr_DistParameter[0],size=int_numDays))
    else:
        print('Error generating the demand schedule.')
    
    return demandSchedule

 ### END "demandScheduleGenerator_EndNode" ###       

def testingScheduleGenerator(nodes = [], int_numDays=1000,
                             int_sampleBudget = 1000,
                             int_PolicyType = 0, 
                             arr_PolicyParameter = [0]
                             ):
    """
    Generates a testing schedule list for a desired number of time units, 
    with the following parameters entered:
            0) nodes: The node list generated from the generateNodeListsFromFile function
            1) int_numDays: Number of schedule days desired
            2) int_sampleBudget: Budget amount, in number of samples
            3) int_PolicyType: Desired policy type, one of the following
                0 = Deterministic-End Node; samples 1 end node each day
                ...
            4) arr_PolicyParameter: Array of the parameters required for generatng
                the test schedules
                ...
    
    Outputs a Python list with the following elements within each entry:
            1) Day: Simulation day of the test
            2) Node: Which node to test that day
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    if int_PolicyType == 0: # Deterministic End-Node policy
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        # Generate a sampling schedule iterating through each end node
        nodeCount = 0
        currNode = endNodes[nodeCount]
        lastEndNode = endNodes[-1]            
        for i in range(int_numDays): 
            sampleSchedule.append([i,currNode])
            if currNode == lastEndNode:
                nodeCount = 0
                currNode = endNodes[nodeCount]
            else:
                nodeCount += 1
                currNode = endNodes[nodeCount]
        ### END DETERMINISTIC END-NODE POLICY
    
    elif int_PolicyType == 1: # Deterministic End-Node policy with multiple tests per day, set via int_sampleBudget
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        # Generate a sampling schedule iterating through each end node
        nodeCount = 0
        currNode = endNodes[nodeCount]
        lastEndNode = endNodes[-1]        
        for i in range(int_sampleBudget):
            day = np.mod(i,int_numDays)
            sampleSchedule.append([day,currNode])
            if currNode == lastEndNode:
                nodeCount = 0
                currNode = endNodes[nodeCount]
            else:
                nodeCount += 1
                currNode = endNodes[nodeCount]
        ### END DETERMINISTIC END-NODE POLICY WITH MULTIPLE TESTS
    
    elif int_PolicyType == 2: # Randomized End-Node policy with multiple tests per day, set via int_sampleBudget
        # Identify the list of end nodes
        endNodes = []
        for nodeInd in range(len(nodes)):
            if nodes[nodeInd][1] == 'TRUE': #We have an end node
                endNodes.append(nodeInd)
        numEndNodes = len(endNodes)
        
        # Generate a sampling schedule randomly sampling the list of end nodes
        for i in range(int_sampleBudget):
            day = np.mod(i,int_numDays)
            currEndInd = int(np.floor(np.random.uniform(low=0,high=numEndNodes,size=1)))
            currNode = endNodes[currEndInd]
            sampleSchedule.append([day,currNode])
        ### END RANODMIZED END-NODE POLICY WITH MULTIPLE TESTS
    
    else:
        print('Error generating the sampling schedule.')
    
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule

 ### END "testingScheduleGenerator" ###

def dynamicTestingGenerator(resultsList,
                             int_totalDays = 1000,
                             int_numDaysRemain = 1000,
                             int_totalBudget = 1000,
                             int_sampleBudgetRemain = 1000,
                             int_PolicyType = 0, 
                             arr_PolicyParameter = [0]
                             ):
    """
    Generates a dynamic testing schedule list for the following day, using the following:
            0) resultsList: A list of testing results so far, with each entry formatted [Node ID, Num Samples, Num Positive, Positive Rate]
            1) int_totalDays: Number of simulation days started with
            2) int_numDaysRemain: Number of simulation days remaining
            3) int_totalBudget: Budget amount, in number of samples
            4) int_sampleBudgetRemain: Remaining budget, in number of samples
            5) int_PolicyType: Desired policy type, one of the following:
                0 = Epsilon-Greedy; with probability epsilon, exploit one of the
                nodes with the highest positive rate
                1 = Epsilon-Decreasing; (epsilon cubed) begin mostly exploring nodes, then overtime have a greater chance of exploitation node with highest positive rate (J curve of epsilon values) 
                3 = Epsilon-First; fully explore for eps*int_totalDays. Then exploit for the remaining (1-eps)*int_totalDays
                5 = Every other; explore when creating even sampling schedules, exploit when creating odd sampling schedules. 
                6 = Epsilon-Decreasing; (square root epsilon) begin mostly exploring nodes, then overtime have a greater chance of exploitation node with highest positive rate (rcurve of epsilon values) 
            6) arr_PolicyParameter: Array of the parameters required for generatng
                the test schedules
                ...
    
    Outputs a Python list with the following elements within each entry:
            1) Day: Simulation day of the test
            2) Node: Which node to test that day
    """
    #Initialize our output, a list with the above mentioned outputs
    sampleSchedule = []
    
    if int_PolicyType == 0: # Basic epsilon-greedy algorithm
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
        ### END EPSILON-DECREASING STRATEGY    
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


    else:
        print('Error generating the sampling schedule.')
    
    # Need to sort this list before passing it through
    sampleSchedule.sort(key=lambda x: x[0])
    return sampleSchedule


 ### END "dynamicTestingGenerator" ###

 ########################### SF RATE ESTIMATORS ###########################
def Est_LinearProjection(A,X): # Linear Projection
    # Uses the (estimated) transition matrix, A, and the (estimated) percentage SF
    # at each end node, X
    intProj = np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,X))
    endProj = np.subtract(X,np.dot(A,intProj))
    return np.ndarray.tolist(intProj.T), np.ndarray.tolist(endProj.T)

def Est_BernMLEProjection(A,X): #MLE OF BERNOULLI VARIABLE
    # USING ITERATIVELY REWEIGHTED LEAST SQUARES, SEE WIKIPEDIA FOR NOTATION
    A = np.array(A)
    X = np.array(X)
    np.seterr(all='print')
    currGap = 10
    tol = 1e-2
    n = A.shape[0] # Number of end nodes
    m = A.shape[1] # Number of intermediate nodes
    X = np.reshape(X,(n,1))
    w_k = np.zeros([m,1])
    while currGap > tol:
        mu_k = []
        for i in range(n):
            mu_k.append(float(1/(1+np.exp(-1*((np.dot(w_k.T,A[i])))))))
        Sdiag = []
        for i in range(n):
            Sdiag.append(mu_k[i]*(1-mu_k[i]))            
        mu_k = np.reshape(mu_k,(n,1))
        S_k = np.diag(Sdiag)
        w_k1 = np.dot(np.linalg.inv(np.dot(A.T,np.dot(S_k,A))),np.dot(A.T, np.subtract(np.add(np.dot(np.dot(S_k,A),w_k),X),mu_k)))
        currGap = np.linalg.norm(w_k-w_k1)
        w_k = np.copy(w_k1)
    # Now our importer SF rates are calculated; figure out variance + Wald statistics
    covarMat_Bern = np.linalg.inv(np.dot(A.T,np.dot(S_k,A)))
    w_Var = np.diag(covarMat_Bern)
    wald_stats = []
    for j in range(m):
        wald_stats.append(float((w_k[j]**2)/w_Var[j]))
    
    # Convert to intermediate and end node estimates
    intProj = np.ndarray.tolist(invlogit(w_k.T.tolist()[0]))
    errs_Bern = np.subtract(X,mu_k)
    endProj = errs_Bern.T.tolist()[0]    
    return intProj, endProj, covarMat_Bern, wald_stats

def PlumleeEstimates(ydata, numsamples, A, sens, spec, rglrWt = 0.1):
    ydata = np.array(ydata)
    numsamples = np.array(numsamples)
    beta0 = -6 * np.ones(A.shape[1]+A.shape[0])
    def invlogit_INTERIOR(beta):
        return np.exp(beta)/(np.exp(beta)+1)
    def mynegloglik_INTERIOR(beta, ydata, numsamples, A, sens, spec):
        betaI = beta[0:A.shape[1]]
        betaJ = beta[A.shape[1]:]
        probs = (1-invlogit_INTERIOR(betaJ)) * np.matmul(A,invlogit_INTERIOR(betaI)) + invlogit_INTERIOR(betaJ)
        probsz = probs*sens + (1-probs) * (1-spec)
        return -np.sum(ydata * np.log(probsz) + (numsamples-ydata) * np.log(1-probsz)) \
            + rglrWt*np.sum(np.abs((betaJ - beta0[A.shape[1]:]))) #have to regularize to prevent problems
    
    bds = spo.Bounds(beta0-8, beta0+8)
    opval = spo.minimize(mynegloglik_INTERIOR, beta0+1,
                         args=(ydata, numsamples, A, sens, spec),
                         method='L-BFGS-B',
                         options={'disp': False},
                         bounds=bds)
    return invlogit_INTERIOR(opval.x)[0:A.shape[1]].tolist(), invlogit_INTERIOR(opval.x)[A.shape[1]:]

########################### END SF RATE ESTIMATORS ###########################

def SimReplicationOutput(OPdict):
    """
    Generates output tables and plots for a given output dictionary, OPdict.
    Each element of the output dictionary should be a dictionary for a given 
    simulation replication containing the following keys:
            0) 'rootConsumption': The root node consumption list
            1) 'intDemandResults': Intermediate nodes demand results
            2) 'endDemandResults': End nodes demand results
            3) 'testResults': The list of test results, which comprises entries
                where each entry is [simDay,testedNode,testResult], where the 
                testResult is stored as the genesis node for the sample
                procured; a testResult of -1 means there were no samples
                available when the test was conducted.
            4) 'intFalseEstimates': The list of calculated falsification
                estimates for the intermediate nodes, using p=((A'A)^(-1))A'X,
                where A is the estimated transition matrix between end nodes and
                intermediate nodes, X is the observed falsification rate at the
                end nodes, and p is the estimate for intermediate nodes
            5) 'endFalseEstimates': X, as given in 4)
    """

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
    
    # Generate summaries of our falsification estimates
    intFalseEstVec = []
    endFalseEstVec = []
    intFalseEstVec_Plum = []
    endFalseEstVec_Plum = []
    for rep in OPdict.keys():
        currIntVec = OPdict[rep]['intFalseEstimates']
        intFalseEstVec.append(currIntVec)
        currEndVec = OPdict[rep]['endFalseEstimates']
        endFalseEstVec.append(currEndVec)
        currIntVec_Plum = OPdict[rep]['intFalseEstimates_Plum']
        intFalseEstVec_Plum.append(currIntVec_Plum)
        currEndVec_Plum = OPdict[rep]['endFalseEstimates_Plum']
        endFalseEstVec_Plum.append(currEndVec_Plum)
    
    
    # For our plots' x axes
    numRoots = len(OPdict[0]['rootConsumption'])
    numInts = len(OPdict[0]['intDemandResults'])
    numEnds = len(OPdict[0]['endDemandResults'])
    Root_Plot_x = []
    for i in range(numRoots):
        Root_Plot_x.append(str(i))
    Int_Plot_x = []
    for i in range(numInts):
        Int_Plot_x.append(str(i+numRoots))
    End_Plot_x = []
    for i in range(numEnds):
        End_Plot_x.append(str(i+numRoots+numInts))
    
    
    
    '''
    currOutputLine = {'rootConsumption':List_RootConsumption,
                          'intDemandResults':List_demandResultsInt,
                          'endDemandResults':List_demandResultsEnd,
                          'testResults':TestReportTbl,
                          'intFalseEstimates':estIntFalsePercList,
                          'endFalseEstimates':estEndFalsePercList}
    '''
    
    
    # PLOTS
    lowErrInd = int(np.floor(0.05*len(OPdict.keys())))
    upErrInd = int(np.ceil(0.95*len(OPdict.keys())))-1
     
    avgFalseConsumedVec.sort()
    # Root node consumption
    rootNode1_mean = np.mean(avgFalseConsumedVec)
    rootNode0_mean = 1-rootNode1_mean
    # Calculate the standard deviation
    rootNode1_lowErr = rootNode1_mean-avgFalseConsumedVec[int(np.floor(0.05*len(avgFalseConsumedVec)))] 
    rootNode1_upErr = avgFalseConsumedVec[int(np.ceil(0.95*len(avgFalseConsumedVec)))-1]-rootNode1_mean 
    rootNode0_lowErr = rootNode1_upErr 
    rootNode0_upErr = rootNode1_lowErr
    # Define positions, bar heights and error bar heights
    means = [rootNode0_mean, rootNode1_mean]
    error = [[rootNode0_lowErr,rootNode1_lowErr], [rootNode0_upErr,rootNode1_upErr]] 
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,0.3,0.5])
    ax.bar(Root_Plot_x, means,
           yerr=error,
           align='center',
           ecolor='black',
           capsize=10,
           color='thistle',edgecolor='indigo')
    ax.set_xlabel('Root Node',fontsize=16)
    ax.set_ylabel('Percentage consumption',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.show()
    
    # Intermediate node stockouts
    intNode_SOs = []
    for i in range(numInts):
        repsAvgVec = []
        for rep in intDemandVec:
            newRow = rep[i]
            newSOPerc = newRow[1]/(newRow[0]+newRow[1])
            repsAvgVec.append(newSOPerc)
        repsAvgVec.sort() 
        intNode_SOs.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    means = [np.mean(x) for x in intNode_SOs] 
    error = [[np.mean(impVec)-impVec[lowErrInd] for impVec in intNode_SOs], 
              [impVec[upErrInd]-np.mean(impVec) for impVec in intNode_SOs]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,0.5])
    ax.bar(Int_Plot_x, means,yerr=error,align='center',ecolor='black',
           capsize=5,color='bisque',edgecolor='darkorange')
    ax.set_xlabel('Intermediate Node',fontsize=16)
    ax.set_ylabel('Percentage stocked out',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.show()
    
    # End node stockouts
    endNode_SOs = []
    for i in range(numEnds):
        repsAvgVec = []
        for rep in endDemandVec:
            newRow = rep[i]
            newSOPerc = newRow[1]/(newRow[0]+newRow[1])
            repsAvgVec.append(newSOPerc)
        repsAvgVec.sort() 
        endNode_SOs.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    endNode_means = [np.mean(x) for x in endNode_SOs] 
    endNode_err = [[np.mean(endVec)-endVec[lowErrInd] for endVec in endNode_SOs], 
              [endVec[upErrInd]-np.mean(endVec) for endVec in endNode_SOs]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,3,0.5])
    ax.bar(End_Plot_x, endNode_means,yerr=endNode_err,align='center',
           ecolor='black',capsize=2,
           color='mintcream',edgecolor='mediumseagreen')
    ax.set_xlabel('End Node',fontsize=16)
    ax.set_ylabel('Percentage stocked out',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.xticks(rotation=90)
    plt.show()
    
    # Testing results
    #### NEED TO DO LATER/FIGURE OUT
    #################################
    
    # Intermediate nodes falsification estimates
    intNodeFalseEsts = []
    for i in range(numInts):
        repsAvgVec = []
        for rep in intFalseEstVec:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        repsAvgVec.sort() 
        intNodeFalseEsts.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    intEst_means = [np.mean(x) for x in intNodeFalseEsts] 
    intEst_err = [[np.mean(intEstVec)-intEstVec[lowErrInd] for intEstVec in intNodeFalseEsts], 
              [intEstVec[upErrInd]-np.mean(intEstVec) for intEstVec in intNodeFalseEsts]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,0.5])
    ax.bar(Int_Plot_x, intEst_means,yerr=intEst_err,
           align='center',ecolor='black',
           capsize=5,color='lightcoral',edgecolor='firebrick')
    ax.set_xlabel('Intermediate Node',fontsize=16)
    ax.set_ylabel('Est. falsification %',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.show()
    
    # Intermediate nodes falsification estimates - Plumlee model
    intNodeFalseEsts_Plum = []
    for i in range(numInts):
        repsAvgVec = []
        for rep in intFalseEstVec_Plum:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        repsAvgVec.sort() 
        intNodeFalseEsts_Plum.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    intEstPlum_means = [np.mean(x) for x in intNodeFalseEsts_Plum] 
    intEstPlum_err =   [[np.mean(intEstVec)-intEstVec[lowErrInd] for intEstVec in intNodeFalseEsts_Plum], 
                       [intEstVec[upErrInd]-np.mean(intEstVec) for intEstVec in intNodeFalseEsts_Plum]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,0.5])
    ax.bar(Int_Plot_x, intEstPlum_means,yerr=intEstPlum_err,
           align='center',ecolor='black',
           capsize=5,color='navajowhite',edgecolor='darkorange')
    ax.set_xlabel('Intermediate Node',fontsize=16)
    ax.set_ylabel('Est. falsification %',fontsize=16)
    #vals = ax.get_yticks()
    #ax.set_yticklabels(['{:,.0}'.format(x) for x in vals])
    plt.show()
    
    # End nodes falsification estimates
    endNodeFalseEsts = []
    for i in range(numEnds):
        repsAvgVec = []
        for rep in endFalseEstVec:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        repsAvgVec.sort() 
        endNodeFalseEsts.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    endEst_means = [np.mean(x) for x in endNodeFalseEsts] 
    endEst_err = [[np.mean(endEstVec)-endEstVec[lowErrInd] for endEstVec in endNodeFalseEsts], 
              [endEstVec[upErrInd]-np.mean(endEstVec) for endEstVec in endNodeFalseEsts]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,3,0.5])
    ax.bar(End_Plot_x, endEst_means,yerr=endEst_err,
           align='center',ecolor='black',
           capsize=1,color='aliceblue',edgecolor='dodgerblue')
    ax.set_xlabel('End Node',fontsize=16)
    ax.set_ylabel('Est. falsification %',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.xticks(rotation=90)
    plt.show()
    
    # End nodes falsification estimates - Plumlee model
    endNodeFalseEsts_Plum = []
    for i in range(numEnds):
        repsAvgVec = []
        for rep in endFalseEstVec_Plum:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        repsAvgVec.sort() 
        endNodeFalseEsts_Plum.append(repsAvgVec)
    # Define positions, bar heights and error bar heights
    endEstPlum_means = [np.mean(x) for x in endNodeFalseEsts_Plum] 
    endEstPlum_err = [[np.mean(endEstVec)-endEstVec[lowErrInd] for endEstVec in endNodeFalseEsts_Plum], 
                      [endEstVec[upErrInd]-np.mean(endEstVec) for endEstVec in endNodeFalseEsts_Plum]]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,3,0.5])
    ax.bar(End_Plot_x, endEstPlum_means,yerr=endEstPlum_err,
           align='center',ecolor='black',
           capsize=1,color='mintcream',edgecolor='forestgreen')
    ax.set_xlabel('End Node',fontsize=16)
    ax.set_ylabel('Est. falsification %',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.xticks(rotation=90)
    plt.show()
    
    
    '''
    alphaLevel = 0.8
    g1 = (avgFalseConsumedVec,avgFalseTestedVec)
    g2 = (avgFalseConsumedVec,avgStockoutTestedVec)
    
    # Plot of testing SF rate vs underlying SF rate
    color = ("red")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x, y = g1
    ax.scatter(x, y, alpha=alphaLevel, c=color, edgecolors='none', s=30)
    lims = [np.min(avgFalseConsumedVec), 
            np.max(avgFalseConsumedVec)]
    ax.plot(lims, lims, 'k-', alpha=0.25, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim([lims[0]-0.01,lims[1]+0.01])
    ax.set_ylim([lims[0]-0.01,lims[1]+0.01])
    ax.set_xlabel('True SF rate', fontsize=12)
    ax.set_ylabel('Test result SFs', fontsize=12)
    plt.title(r'Test results of SF FOUND vs. Underlying SF consumption rates', fontsize=14)
    plt.show()
    
    # Plot of testing stockout rate vs underlying SF rate
    color = ("blue")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x, y = g2
    ax.scatter(x, y, alpha=alphaLevel, c=color, edgecolors='none', s=30)
    ax.set_xlabel('True SF rate', fontsize=12)
    ax.set_ylabel('Test result stockouts', fontsize=12)
    plt.title(r'Test results of STOCKOUTS vs. Underlying SF consumption rates', fontsize=14)
    plt.show()
    '''
 ### END "SimReplicationOutput" ###

def SimSFEstimateOutput(OPdicts):
    '''
    Generates comparison tables and plots for a LIST of output dictionaries.
    Intended for comparison with the underlying "true" SF rates at the
    importer level.
    '''
    numDicts = len(OPdicts)
    scenarioList = [] # Initialize a list of possible 'true' underyling SF rates
    # Initialize deviation lists; contains lists of deviations for each replication
    avgDevList_Lin = []
    avgDevList_Bern = []
    avgDevList_MLE = []
    stdDevList_Lin = []
    stdDevList_Bern = []
    stdDevList_MLE = []    
    
    # For each output dictionary, generate deviation estimates of varying types
    for currDict in OPdicts:
        # Loop through each replication contained in the current output dictionary
        currDict_avgDevList_Lin = []
        currDict_avgDevList_Bern = []
        currDict_avgDevList_MLE = []
        currDict_stdDevList_Lin = []
        currDict_stdDevList_Bern = []
        currDict_stdDevList_MLE = []
        
        for repNum in currDict.keys():
            currTrueSFVec = currDict[repNum]['intSFTrueValues']
            for scen in currTrueSFVec:
                if not scen in scenarioList:
                    scenarioList.append(scen)
                    scenarioList.sort()
            currLinProj = currDict[repNum]['intFalseEstimates']
            currBernProj = currDict[repNum]['intFalseEstimates_Bern']
            currMLEProj = currDict[repNum]['intFalseEstimates_Plum']
            currLinProjdevs = [currLinProj[i]-currTrueSFVec[i] for i in range(len(currTrueSFVec))]
            currBernProjdevs = [currBernProj[i]-currTrueSFVec[i] for i in range(len(currTrueSFVec))]
            currMLEProjdevs = [currMLEProj[i]-currTrueSFVec[i] for i in range(len(currTrueSFVec))]
            
            currDict_avgDevList_Lin.append(np.mean(currLinProjdevs))
            currDict_avgDevList_Bern.append(np.mean(currBernProjdevs))
            currDict_avgDevList_MLE.append(np.mean(currMLEProjdevs))
            currDict_stdDevList_Lin.append(np.std(currLinProjdevs))
            currDict_stdDevList_Bern.append(np.std(currBernProjdevs))
            currDict_stdDevList_MLE.append(np.std(currMLEProjdevs))
        
        avgDevList_Lin.append(currDict_avgDevList_Lin)
        avgDevList_Bern.append(currDict_avgDevList_Bern)
        avgDevList_MLE.append(currDict_avgDevList_MLE)
        stdDevList_Lin.append(currDict_stdDevList_Lin)
        stdDevList_Bern.append(currDict_stdDevList_Bern)
        stdDevList_MLE.append(currDict_stdDevList_MLE)
    
    # Scenario-dependent looks at performance
    scenDict = {}
    ind = 0 
    for currScen in scenarioList:
        SCENavgDevList_Lin = []
        SCENavgDevList_Bern = []
        SCENavgDevList_MLE = []
        SCENstdDevList_Lin = []
        SCENstdDevList_Bern = []
        SCENstdDevList_MLE = []
        for currDict in OPdicts:
        # Loop through each replication contained in the current output dictionary
            currScenDict_avgDevList_Lin = []
            currScenDict_avgDevList_Bern = []
            currScenDict_avgDevList_MLE = []
            currScenDict_stdDevList_Lin = []
            currScenDict_stdDevList_Bern = []
            currScenDict_stdDevList_MLE = []
            
            for repNum in currDict.keys():
                currTrueSFVec = currDict[repNum]['intSFTrueValues']
                currScenInds = [i for i, val in enumerate(currTrueSFVec) if val == currScen]
                if not not currScenInds: #Only find deviations if the scenario was used (list is nonempty)
                    currScenTrueSFVec = [currTrueSFVec[i] for i in currScenInds]
                    currScenLinProj = [currDict[repNum]['intFalseEstimates'][i] for i in currScenInds]
                    currScenBernProj = [currDict[repNum]['intFalseEstimates_Bern'][i] for i in currScenInds]
                    currScenMLEProj = [currDict[repNum]['intFalseEstimates_Plum'][i] for i in currScenInds]
                    currScenLinProjdevs = [currScenLinProj[i]-currScenTrueSFVec[i] for i in range(len(currScenTrueSFVec))]
                    currScenBernProjdevs = [currScenBernProj[i]-currScenTrueSFVec[i] for i in range(len(currScenTrueSFVec))]
                    currScenMLEProjdevs = [currScenMLEProj[i]-currScenTrueSFVec[i] for i in range(len(currScenTrueSFVec))]
                    
                    currScenDict_avgDevList_Lin.append(np.mean(currScenLinProjdevs))
                    currScenDict_avgDevList_Bern.append(np.mean(currScenBernProjdevs))
                    currScenDict_avgDevList_MLE.append(np.mean(currScenMLEProjdevs))
                    currScenDict_stdDevList_Lin.append(np.std(currScenLinProjdevs))
                    currScenDict_stdDevList_Bern.append(np.std(currScenBernProjdevs))
                    currScenDict_stdDevList_MLE.append(np.std(currScenMLEProjdevs))
            
            SCENavgDevList_Lin.append(currScenDict_avgDevList_Lin)
            SCENavgDevList_Bern.append(currScenDict_avgDevList_Bern)
            SCENavgDevList_MLE.append(currScenDict_avgDevList_MLE)
            SCENstdDevList_Lin.append(currScenDict_stdDevList_Lin)
            SCENstdDevList_Bern.append(currScenDict_stdDevList_Bern)
            SCENstdDevList_MLE.append(currScenDict_stdDevList_MLE)
        
        currOutputLine = {'scenario': currScen,
                          'SCENavgDevList_Lin':SCENavgDevList_Lin,
                          'SCENavgDevList_Bern':SCENavgDevList_Bern,
                          'SCENavgDevList_MLE':SCENavgDevList_MLE,
                          'SCENstdDevList_Lin':SCENstdDevList_Lin,
                          'SCENstdDevList_Bern':SCENstdDevList_Bern,
                          'SCENstdDevList_MLE':SCENstdDevList_MLE
                          }
        scenDict[ind] = currOutputLine
        ind += 1
        
    ### Plots of our results
    #Some initial plotting needs
    #Indices
    lowErrIndVec = []
    upErrIndVec = []
    for dictNum in range(numDicts):
        lowErrIndVec.append(int(np.floor(0.05*len(OPdicts[dictNum].keys()))))
        upErrIndVec.append(int(np.ceil(0.95*len(OPdicts[dictNum].keys())))-1)
    #x-axis title
    Plot_XLabels = [str(i) for i in range(numDicts)]
    
    # Linear projection - avgs
    dictAvgDevs_Lin = []
    dictAvgLows_Lin = []
    dictAvgUps_Lin = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currAvgDevList = avgDevList_Lin[i]
        currAvgDevList.sort()
        meanDev = np.mean(currAvgDevList)
        lowDev = currAvgDevList[lowErrInd]
        upDev = currAvgDevList[upErrInd]
        dictAvgDevs_Lin.append(meanDev)
        dictAvgLows_Lin.append(lowDev)
        dictAvgUps_Lin.append(upDev)
    dictAvgErr_Lin = [[dictAvgDevs_Lin[i]-dictAvgLows_Lin[i] for i in range(numDicts)], 
              [dictAvgUps_Lin[j]-dictAvgDevs_Lin[j] for j in range(numDicts)]]
    
    # Bernoulli - avgs
    dictAvgDevs_Bern = []
    dictAvgLows_Bern = []
    dictAvgUps_Bern = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currAvgDevList = avgDevList_Bern[i]
        currAvgDevList.sort()
        meanDev = np.mean(currAvgDevList)
        lowDev = currAvgDevList[lowErrInd]
        upDev = currAvgDevList[upErrInd]
        dictAvgDevs_Bern.append(meanDev)
        dictAvgLows_Bern.append(lowDev)
        dictAvgUps_Bern.append(upDev)
    dictAvgErr_Bern = [[dictAvgDevs_Bern[i]-dictAvgLows_Bern[i] for i in range(numDicts)], 
              [dictAvgUps_Bern[j]-dictAvgDevs_Bern[j] for j in range(numDicts)]]
    
    # MLE w Optimizer - avgs
    dictAvgDevs_MLE = []
    dictAvgLows_MLE = []
    dictAvgUps_MLE = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currAvgDevList = avgDevList_MLE[i]
        currAvgDevList.sort()
        meanDev = np.mean(currAvgDevList)
        lowDev = currAvgDevList[lowErrInd]
        upDev = currAvgDevList[upErrInd]
        dictAvgDevs_MLE.append(meanDev)
        dictAvgLows_MLE.append(lowDev)
        dictAvgUps_MLE.append(upDev)
    dictAvgErr_MLE = [[dictAvgDevs_MLE[i]-dictAvgLows_MLE[i] for i in range(numDicts)], 
              [dictAvgUps_MLE[j]-dictAvgDevs_MLE[j] for j in range(numDicts)]]
    
    # Linear projection - standard deviations
    dictStdDevs_Lin = []
    dictStdLows_Lin = []
    dictStdUps_Lin = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currDevList = stdDevList_Lin[i]
        currDevList.sort()
        meanDev = np.mean(currDevList)
        lowDev = currDevList[lowErrInd]
        upDev = currDevList[upErrInd]
        dictStdDevs_Lin.append(meanDev)
        dictStdLows_Lin.append(lowDev)
        dictStdUps_Lin.append(upDev)
    dictStdErr_Lin = [[dictStdDevs_Lin[i]-dictStdLows_Lin[i] for i in range(numDicts)], 
              [dictStdUps_Lin[j]-dictStdDevs_Lin[j] for j in range(numDicts)]]
    
    # Bernoulli - standard deviations
    dictStdDevs_Bern = []
    dictStdLows_Bern = []
    dictStdUps_Bern = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currDevList = stdDevList_Bern[i]
        currDevList.sort()
        meanDev = np.mean(currDevList)
        lowDev = currDevList[lowErrInd]
        upDev = currDevList[upErrInd]
        dictStdDevs_Bern.append(meanDev)
        dictStdLows_Bern.append(lowDev)
        dictStdUps_Bern.append(upDev)
    dictStdErr_Bern = [[dictStdDevs_Bern[i]-dictStdLows_Bern[i] for i in range(numDicts)], 
              [dictStdUps_Bern[j]-dictStdDevs_Bern[j] for j in range(numDicts)]]
    
    # MLE w Optimizer - standard deviations
    dictStdDevs_MLE = []
    dictStdLows_MLE = []
    dictStdUps_MLE = []
    for i in range(numDicts):
        lowErrInd = lowErrIndVec[i]
        upErrInd = upErrIndVec[i]
        currDevList = stdDevList_MLE[i]
        currDevList.sort()
        meanDev = np.mean(currDevList)
        lowDev = currDevList[lowErrInd]
        upDev = currDevList[upErrInd]
        dictStdDevs_MLE.append(meanDev)
        dictStdLows_MLE.append(lowDev)
        dictStdUps_MLE.append(upDev)
    dictStdErr_MLE = [[dictStdDevs_MLE[i]-dictStdLows_MLE[i] for i in range(numDicts)], 
              [dictStdUps_MLE[j]-dictStdDevs_MLE[j] for j in range(numDicts)]]
    
    
    # Build plots
    # Devaition averages (biases)
    fig, axs = plt.subplots(3, 1,figsize=(9,13))
    fig.suptitle('Estimate Deviation AVERAGES',fontsize=18)
    for subP in range(3):
        axs[subP].set_xlabel('Output Dictionary',fontsize=12)
        axs[subP].set_ylabel('Est. Avg. Deviation',fontsize=12)        
        axs[subP].set_ylim([-0.4,0.4])
        #vals = axs[subP].get_yticks()
        #axs[subP].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    # Linear projection
    axs[0].set_title('Linear projection',fontweight='bold')        
    axs[0].errorbar(Plot_XLabels,dictAvgDevs_Lin,yerr=dictAvgErr_Lin,
                   ecolor='mediumpurple',capsize=5,color='indigo')
    # Bernoulli MLE projection
    axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
    axs[1].errorbar(Plot_XLabels,dictAvgDevs_Bern,yerr=dictAvgErr_Bern,
                   ecolor='seagreen',capsize=5,color='green')      
    # MLE w Nonlinear optimizer        
    axs[2].set_title('MLE w/ Nonlinear Optimizer',fontweight='bold')
    axs[2].errorbar(Plot_XLabels,dictAvgDevs_MLE,yerr=dictAvgErr_MLE,
                   ecolor='orange',capsize=5,color='darkorange')   
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()
    
    #Deviation variation
    fig, axs = plt.subplots(3, 1,figsize=(9,13))
    fig.suptitle('Estimate Deviation STANDARD DEVIATIONS',fontsize=18)
    for subP in range(3):
        axs[subP].set_xlabel('Output Dictionary',fontsize=12)
        axs[subP].set_ylabel('Est. StDev of Deviation',fontsize=12)        
        axs[subP].set_ylim([0,0.4])
        #vals = axs[subP].get_yticks()
        #axs[subP].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    # Linear projection
    axs[0].set_title('Linear projection',fontweight='bold')        
    axs[0].errorbar(Plot_XLabels,dictStdDevs_Lin,yerr=dictStdErr_Lin,
                   ecolor='mediumpurple',capsize=5,color='indigo')
    # Bernoulli MLE projection
    axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
    axs[1].errorbar(Plot_XLabels,dictStdDevs_Bern,yerr=dictStdErr_Bern,
                   ecolor='seagreen',capsize=5,color='green')      
    # MLE w Nonlinear optimizer        
    axs[2].set_title('MLE w/ Nonlinear Optimizer',fontweight='bold')
    axs[2].errorbar(Plot_XLabels,dictStdDevs_MLE,yerr=dictStdErr_MLE,
                   ecolor='orange',capsize=5,color='darkorange')   
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.show()
    
    # For each scenario; only find average deviations
    for scenInd in range(len(scenarioList)):
        currScenDict = scenDict[scenInd]
        avgDevListLin_scen = currScenDict['SCENavgDevList_Lin']
        avgDevListBern_scen = currScenDict['SCENavgDevList_Bern']
        avgDevListMLE_scen = currScenDict['SCENavgDevList_MLE']
        
        lowErrIndVec = []
        upErrIndVec = []
        for dictNum in range(numDicts):
            lowErrIndVec.append(int(np.floor(0.05*len(avgDevListLin_scen[dictNum]))))
            upErrIndVec.append(int(np.ceil(0.95*len(avgDevListLin_scen[dictNum])))-1)
        #x-axis title
        Plot_XLabels = [str(i) for i in range(numDicts)]
        
        # Linear projection - avgs
        dictAvgDevs_Lin = []
        dictAvgLows_Lin = []
        dictAvgUps_Lin = []
        for i in range(numDicts):
            lowErrInd = lowErrIndVec[i]
            upErrInd = upErrIndVec[i]
            currAvgDevList = avgDevListLin_scen[i]
            currAvgDevList.sort()
            meanDev = np.mean(currAvgDevList)
            lowDev = currAvgDevList[lowErrInd]
            upDev = currAvgDevList[upErrInd]
            dictAvgDevs_Lin.append(meanDev)
            dictAvgLows_Lin.append(lowDev)
            dictAvgUps_Lin.append(upDev)
        dictAvgErr_Lin = [[dictAvgDevs_Lin[i]-dictAvgLows_Lin[i] for i in range(numDicts)], 
                  [dictAvgUps_Lin[j]-dictAvgDevs_Lin[j] for j in range(numDicts)]]
        
        # Bernoulli - avgs
        dictAvgDevs_Bern = []
        dictAvgLows_Bern = []
        dictAvgUps_Bern = []
        for i in range(numDicts):
            lowErrInd = lowErrIndVec[i]
            upErrInd = upErrIndVec[i]
            currAvgDevList = avgDevListBern_scen[i]
            currAvgDevList.sort()
            meanDev = np.mean(currAvgDevList)
            lowDev = currAvgDevList[lowErrInd]
            upDev = currAvgDevList[upErrInd]
            dictAvgDevs_Bern.append(meanDev)
            dictAvgLows_Bern.append(lowDev)
            dictAvgUps_Bern.append(upDev)
        dictAvgErr_Bern = [[dictAvgDevs_Bern[i]-dictAvgLows_Bern[i] for i in range(numDicts)], 
                  [dictAvgUps_Bern[j]-dictAvgDevs_Bern[j] for j in range(numDicts)]]
        
        # MLE w Optimizer - avgs
        dictAvgDevs_MLE = []
        dictAvgLows_MLE = []
        dictAvgUps_MLE = []
        for i in range(numDicts):
            lowErrInd = lowErrIndVec[i]
            upErrInd = upErrIndVec[i]
            currAvgDevList = avgDevListMLE_scen[i]
            currAvgDevList.sort()
            meanDev = np.mean(currAvgDevList)
            lowDev = currAvgDevList[lowErrInd]
            upDev = currAvgDevList[upErrInd]
            dictAvgDevs_MLE.append(meanDev)
            dictAvgLows_MLE.append(lowDev)
            dictAvgUps_MLE.append(upDev)
        dictAvgErr_MLE = [[dictAvgDevs_MLE[i]-dictAvgLows_MLE[i] for i in range(numDicts)], 
                  [dictAvgUps_MLE[j]-dictAvgDevs_MLE[j] for j in range(numDicts)]]
    
        # Generate plots for each scenario
        # Devaition averages (biases)
        fig, axs = plt.subplots(3, 1,figsize=(9,13))
        fig.suptitle('Estimate Deviation AVERAGES for SF RATE OF '+str(scenarioList[scenInd]),fontsize=18)
        for subP in range(3):
            axs[subP].set_xlabel('Output Dictionary',fontsize=12)
            axs[subP].set_ylabel('Est. Avg. Deviation',fontsize=12)        
            axs[subP].set_ylim([-0.4,0.4])
            #vals = axs[subP].get_yticks()
            #axs[subP].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        # Linear projection
        axs[0].set_title('Linear projection',fontweight='bold')        
        axs[0].errorbar(Plot_XLabels,dictAvgDevs_Lin,yerr=dictAvgErr_Lin,
                       ecolor='mediumpurple',capsize=5,color='indigo')
        # Bernoulli MLE projection
        axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
        axs[1].errorbar(Plot_XLabels,dictAvgDevs_Bern,yerr=dictAvgErr_Bern,
                       ecolor='seagreen',capsize=5,color='green')      
        # MLE w Nonlinear optimizer        
        axs[2].set_title('MLE w/ Nonlinear Optimizer',fontweight='bold')
        axs[2].errorbar(Plot_XLabels,dictAvgDevs_MLE,yerr=dictAvgErr_MLE,
                       ecolor='orange',capsize=5,color='darkorange')   
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.show()
        
        
    # END OF SCENARIOS LOOP
        
        
    
    
    
 ### END "SimSFEstimateOutput" ###   
    
def setWarmUp(useWarmUpFileBool = False, warmUpRunBool = False, numReps = 1,
              currDirect = ''):
    """
    Sets up warm-up files as a function of the chosen parameters.
    Warm-up dictionaries are saved to a folder 'warm up dictionaries' in the
    current working directory.
    """
    warmUpDirectory = ''
    warmUpFileName_str = ''
    warmUpDict = {}
    if useWarmUpFileBool == True and warmUpRunBool == True:
        print('Cannot use warm up files and conduct warm up runs at the same time!')
        useWarmUpFileBool = False
        warmUpRunBool = False
        numReps = 0

    elif useWarmUpFileBool == True and warmUpRunBool == False:
        warmUpDirectory = os.getcwd() + '\\warm up dictionaries' # Location of warm-up files
        warmUpFileName_str =  os.path.basename(sys.argv[0]) # Current file name
        warmUpFileName_str = warmUpFileName_str[:-3] + '_WARM_UP' # Warm-up file name
        warmUpFileName_str = os.path.join(warmUpDirectory, warmUpFileName_str)
        if not os.path.exists(warmUpFileName_str): # Flag if this directory not found
            print('Warm up file not found.')
            numReps = 0
        else:
            with open(warmUpFileName_str, 'rb') as f:
                warmUpDict = pickle.load(f) # Load the dictionary
            
    elif useWarmUpFileBool == False and warmUpRunBool == True: # Generate warm-up runs file
        # Generate a directory if one does not already exist
        warmUpDirectory = os.getcwd() + '\\warm up dictionaries' # Location of warm-up files
        if not os.path.exists(warmUpDirectory): # Generate this folder if one does not already exist
            os.makedirs(warmUpDirectory)
        warmUpFileName_str =  os.path.basename(sys.argv[0]) # Current file name
        warmUpFileName_str = warmUpFileName_str[:-3] + '_WARM_UP' # Warm-up file name
        warmUpFileName_str = os.path.join(warmUpDirectory, warmUpFileName_str)
        if os.path.exists(warmUpFileName_str): # Generate this file if one doesn't exist
            with open(warmUpFileName_str, 'rb') as f:
                warmUpDict = pickle.load(f) # Load the dictionary
        else:
            warmUpDict = {} # Initialize the dictionary
        pickle.dump(warmUpDict, open(warmUpFileName_str,'wb'))
      
    elif useWarmUpFileBool == False and warmUpRunBool == False: # Nothing done WRT warm-ups
        pass  
    
    
    return numReps, warmUpRunBool, useWarmUpFileBool, warmUpDirectory, warmUpFileName_str, warmUpDict
    
#### Likelihood estimate functions
def invlogit(beta):
    return sps.expit(beta)

def invlogit_grad(beta):
    return np.diag(np.exp(beta)/((np.exp(beta)+1) ** 2))

def myloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = (1-invlogit(betaJ)) * np.squeeze(np.matmul(A,invlogit(betaI)))  + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return np.sum(ydata * np.log(probsz) + (np.asarray(nsamp)-np.asarray(ydata)) * np.log(1-probsz))

'''
def myloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    probs = ((1-invlogit(betaJ)) * np.squeeze(np.array(A @ invlogit(betaI)))  + invlogit(betaJ))
    probs_dirJ = -invlogit_grad(betaJ) * np.squeeze(np.array(A @ invlogit(betaI))) + invlogit_grad(betaJ)
    ilJ = np.array([invlogit(betaJ).T,]*(betaI.shape[0])).T
    ilK = np.array([invlogit_grad(betaI),]*(betaJ.shape[0]))
    probs_dirI =  np.multiply(np.multiply((1- ilJ ) , A) , np.squeeze(ilK)).T 
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  np.array(probz_dirI @ (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz))).T)
    negloglikJ =  np.squeeze(np.array((probz_dirJ * (-(np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz)))).T))
    return -np.squeeze(np.concatenate([negloglikI,negloglikJ]))
'''

def mynegloglik_grad(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]
    iliJ = invlogit(betaJ)
    iliJg = np.diag(iliJ * (1-iliJ))
    iliI = invlogit(betaI)
    AiliI = A * iliI
    AiliIg = AiliI * (1-iliI)
    AiliI = np.sum(AiliI,1)
    probs = (1-iliJ) * AiliI + iliJ
    
    probs_dirJ = -iliJg * AiliI + iliJg
    probs_dirI =  AiliIg.T * (1- iliJ)
    
    probz_dirI = probs_dirI*sens - (probs_dirI) * (1-spec)
    probz_dirJ = probs_dirJ*sens - (probs_dirJ) * (1-spec)
    probsz = probs*sens + (1-probs) * (1-spec)
    negloglikI =  -np.sum((np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz)) *  (probz_dirI),1) 
    negloglikJ =  -np.sum((np.asarray(ydata) / probsz  - (np.asarray(nsamp)-np.asarray(ydata)) / (1-probsz)) *  (probz_dirJ),1) 

    return np.concatenate((negloglikI,negloglikJ))

def mynegloglik(beta, ydata, nsamp, A, sens, spec):
    betaI = beta[0:A.shape[1]]
    betaJ = beta[A.shape[1]:]   
    probs = (1-invlogit(betaJ)) * np.matmul(A,invlogit(betaI)) + invlogit(betaJ)
    probsz = probs*sens + (1-probs) * (1-spec)
    return -np.sum(np.asarray(ydata) * np.log(probsz) + (np.asarray(nsamp)-np.asarray(ydata)) * np.log(1-probsz))

def mylogprior(beta, ydata, nsamp, A, sens, spec):
    #betaJ = beta[A.shape[1]:]
    return -0.25*np.sum(np.abs(beta + 3))

def mylogprior_grad(beta, ydata, nsamp, A, sens, spec):
    #betaI = beta[0:A.shape[1]]
    #betaJ = beta[A.shape[1]:]
    return -0.25*np.squeeze(1*(beta >= -3) - 1*(beta <= -3))

def mylogpost(beta, ydata, nsamp, A, sens, spec):
    return mylogprior(beta, ydata, nsamp, A, sens, spec)+myloglik(beta, ydata, nsamp, A, sens, spec)

def mylogpost_grad(beta, ydata, nsamp, A, sens, spec):
    return mylogprior_grad(beta, ydata, nsamp, A, sens, spec)-mynegloglik_grad(beta, ydata, nsamp, A, sens, spec) #modified by EOW

def exampletargetfornuts(beta):
    """
    Example of a target distribution that could be sampled from using NUTS.
    (Although of course you could sample from it more efficiently)
    Doesn't include the normalizing constant.
    """
    return mylogpost(beta,ydata, nsamp, A, sens, spec), mylogpost_grad(beta,ydata, nsamp, A, sens, spec) 


#### Necessary NUTS module ####
"""
This package implements the No-U-Turn Sampler (NUTS) algorithm 6 from the NUTS
paper (Hoffman & Gelman, 2011).

Content
-------

The package mainly contains:
  nuts6                     return samples using the NUTS
  test_nuts6                example usage of this package

and subroutines of nuts6:
  build_tree                the main recursion in NUTS
  find_reasonable_epsilon   Heuristic for choosing an initial value of epsilon
  leapfrog                  Perfom a leapfrog jump in the Hamiltonian space
  stop_criterion            Compute the stop condition in the main loop


A few words about NUTS
----------------------

Hamiltonian Monte Carlo or Hybrid Monte Carlo (HMC) is a Markov chain Monte
Carlo (MCMC) algorithm that avoids the random walk behavior and sensitivity to
correlated parameters, biggest weakness of many MCMC methods. Instead, it takes
a series of steps informed by first-order gradient information.

This feature allows it to converge much more quickly to high-dimensional target
distributions compared to simpler methods such as Metropolis, Gibbs sampling
(and derivatives).

However, HMC's performance is highly sensitive to two user-specified
parameters: a step size, and a desired number of steps.  In particular, if the
number of steps is too small then the algorithm will just exhibit random walk
behavior, whereas if it is too large it will waste computations.

Hoffman & Gelman introduced NUTS or the No-U-Turn Sampler, an extension to HMC
that eliminates the need to set a number of steps.  NUTS uses a recursive
algorithm to find likely candidate points that automatically stops when it
starts to double back and retrace its steps.  Empirically, NUTS perform at
least as effciently as and sometimes more effciently than a well tuned standard
HMC method, without requiring user intervention or costly tuning runs.

Moreover, Hoffman & Gelman derived a method for adapting the step size
parameter on the fly based on primal-dual averaging.  NUTS can thus be used
with no hand-tuning at all.

In practice, the implementation still requires a number of steps, a burning
period and a stepsize. However, the stepsize will be optimized during the
burning period, and the final values of all the user-defined values will be
revised by the algorithm.

reference: arXiv:1111.4246
"The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte
Carlo", Matthew D. Hoffman & Andrew Gelman
"""

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



def find_reasonable_epsilon(theta0, grad0, logp0, f, epsilonLB = 0.05, epsilonUB = 0.5):
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

    print("find_reasonable_epsilon=", epsilon)

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
            s = sprime and stop_criterion(thetaminus, thetaplus, rminus, rplus) and (n < 50)
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
        
        ### EOW ADDITION
        if np.mod(m,10)==0:
            print('Round '+str(m)+' of nuts6 finished' )
        ### END EOW
        
    samples = samples[Madapt:, :]
    lnprob = lnprob[Madapt:]
    return samples, lnprob, epsilon






