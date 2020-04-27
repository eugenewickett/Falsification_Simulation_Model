# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 17:04:36 2019

@author: Eugene Wickett

Stores modules for use with 'SC Simulator.py'
"""

import numpy as np
import random
import csv
import os
import sys
import pickle
from tabulate import tabulate # for making outputs
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
    for rep in OPdict.keys():
        currIntVec = OPdict[rep]['intFalseEstimates']
        intFalseEstVec.append(currIntVec)
        currEndVec = OPdict[rep]['endFalseEstimates']
        endFalseEstVec.append(currEndVec)
    
    
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
    # Root node consumption
    rootNode1_mean = np.mean(avgFalseConsumedVec)
    rootNode0_mean = 1-rootNode1_mean
    # Calculate the standard deviation
    rootNode1_std = np.std(avgFalseConsumedVec)
    rootNode0_std = rootNode1_std
    # Define positions, bar heights and error bar heights
    means = [rootNode0_mean, rootNode1_mean]
    error = [1.6*rootNode0_std, 1.6*rootNode1_std]
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
    intNode_means = []
    intNode_stds = []
    for i in range(numInts):
        repsAvgVec = []
        for rep in intDemandVec:
            newRow = rep[i]
            newSOPerc = newRow[1]/(newRow[0]+newRow[1])
            repsAvgVec.append(newSOPerc)
        intNode_means.append(np.mean(repsAvgVec))
        intNode_stds.append(np.std(repsAvgVec))
    # Define positions, bar heights and error bar heights
    means = intNode_means
    error = [x*1.6 for x in intNode_stds]
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
    endNode_means = []
    endNode_stds = []
    for i in range(numEnds):
        repsAvgVec = []
        for rep in endDemandVec:
            newRow = rep[i]
            newSOPerc = newRow[1]/(newRow[0]+newRow[1])
            repsAvgVec.append(newSOPerc)
        endNode_means.append(np.mean(repsAvgVec))
        endNode_stds.append(np.std(repsAvgVec))
    # Define positions, bar heights and error bar heights
    endNode_stds = [x*1.6 for x in endNode_stds]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,3,0.5])
    ax.bar(End_Plot_x, endNode_means,yerr=endNode_stds,align='center',
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
    intNodeFalse_means = []
    intNodeFalse_stds = []
    for i in range(numInts):
        repsAvgVec = []
        for rep in intFalseEstVec:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        intNodeFalse_means.append(np.mean(repsAvgVec))
        intNodeFalse_stds.append(np.std(repsAvgVec))
    # Define positions, bar heights and error bar heights
    intNodeFalse_stds = [x*1.6 for x in intNodeFalse_stds]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,0.5])
    ax.bar(Int_Plot_x, intNodeFalse_means,yerr=intNodeFalse_stds,
           align='center',ecolor='black',
           capsize=5,color='lightcoral',edgecolor='firebrick')
    ax.set_xlabel('Intermediate Node',fontsize=16)
    ax.set_ylabel('Estimated falsification percentage',fontsize=16)
    vals = ax.get_yticks()
    ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
    plt.show()
    
    # End nodes falsification estimates
    endNodeFalse_means = []
    endNodeFalse_stds = []
    for i in range(numEnds):
        repsAvgVec = []
        for rep in endFalseEstVec:
            newItem = rep[i]
            repsAvgVec.append(newItem)
        endNodeFalse_means.append(np.mean(repsAvgVec))
        endNodeFalse_stds.append(np.std(repsAvgVec))
    # Define positions, bar heights and error bar heights
    endNodeFalse_stds = [x*1.6 for x in endNodeFalse_stds]
    # Build the plot
    fig = plt.figure()
    ax = fig.add_axes([0,0,3,0.5])
    ax.bar(End_Plot_x, endNodeFalse_means,yerr=endNodeFalse_stds,
           align='center',ecolor='black',
           capsize=1,color='aliceblue',edgecolor='dodgerblue')
    ax.set_xlabel('End Node',fontsize=16)
    ax.set_ylabel('Estimated falsification percentage',fontsize=16)
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
    
    
  