# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:04:36 2019

@author: Eugene Wickett

Stores modules for use with 'SC Simulator.py'
"""

import numpy as np
import random
import csv
import matplotlib.pyplot as plt

import Falsification_Sim_Classes as simClasses # modules for the simulation

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
                # Add a node object to the List_IntermediateNode
                newIntermediateNode = simClasses.Node(nodeNum=currNode, reorderPoint=currReorderPoint, reorderAmount=currReorderAmount, preferredSupplierVec=currPreferredSuppliers, preferredSupplierLTsVec=currPreferredSuppliersLTs, preferredSupplierRsVec=currPreferredSuppliersRs, endNodeTF = False)
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
    
    
    
    
    
    
    
    
    
    
    
    # How different in stockout rates is lowest int node
    #avgIntStockout1 = []
    #for item in intDemandVec1:
    #    currLowRate = 1
    #    rateVec = []
    #    for intNode in item:
    #        rate = intNode[1]/(intNode[1]+intNode[0])
    #        rateVec.append(rate)
    #        if rate < currLowRate:
    #            currLowRate = rate
    #    currDiff = currLowRate - np.mean(rateVec)
    #    avgIntStockout1.append(currDiff)
        
    
    
    # Variability in tested falsification rates
    
    
    
    # PLOTS 
    # Histogram of SF consumption rates across replications
    plt.title(r'Histogram of SF consumption rates')
    plt.hist(avgFalseConsumedVec)
    plt.show()
    #plt.title(r'Histogram Intermediate stockout rates')
    #plt.hist(avgIntStockout1)
    
    
    alphaLevel = 0.8
    g1 = (avgFalseConsumedVec,avgFalseTestedVec)
    g2 = (avgFalseConsumedVec,avgStockoutTestedVec)
    
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
    
    color = ("blue")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    x, y = g2
    ax.scatter(x, y, alpha=alphaLevel, c=color, edgecolors='none', s=30)
    ax.set_xlabel('True SF rate', fontsize=12)
    ax.set_ylabel('Test result stockouts', fontsize=12)
    plt.title(r'Test results of STOCKOUTS vs. Underlying SF consumption rates', fontsize=14)
    plt.show()

 ### END "SimReplicationOutput" ###
