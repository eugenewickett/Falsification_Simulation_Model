# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:59:35 2019

@author: Eugene Wickett

This code is for building and running simulations of a supply chain susceptible
to falsification and substandardization.
"""

import numpy as np
import matplotlib.pyplot as plt
import random #for seeds
#import sys
import time # for time tracking
import os # for directories
from tabulate import tabulate # for making outputs
import pickle # for saving/loading objects in Python
import scipy.special as sps

#import Falsification_Sim_Classes as simClasses # our class objects and methods
import SFP_Sim_Helpers as simHelpers # module for the simulation
import SFP_Sim_EstimationMethods as simEstMethods # estimation methods
import SFP_Sim_TestPolicies as simTestPolicies # testing policies

# Run supporting files
currDirectory = os.path.dirname(os.path.realpath(__file__)) #Run this command to get the current working directory string
os.chdir(currDirectory) # Set directory

##### GRAPH INPUT FILES
# Read input lists: node list, arc preferences, and arc lead times
nodeInputFileString = 'LIB_Nodes_1.csv'
arcPreferencesFileString = 'LIB_Arcs_Preferences_1.csv'
arcLTsFileString = 'LIB_Arcs_LTs_1.csv'
arcRsFileString = 'LIB_Arcs_Rs_1.csv'

##### SIMULATION PARAMETERS
# Enter the length of the simulation and the sampling budget
NumSimDays = 600
samplingBudget = NumSimDays*1
diagnosticSensitivity = 0.95 # Tool sensitivity
diagnosticSpecificity = 0.98 # Tool specificity
globalDemand = 0. #Level of global demand increase across all outlets, in mean demand/simulation day
numReplications = 1

alertIter = 7 # How frequently we're alerted of a set of replications being completed
printOutput = True # Whether individual replication output should be displayed
storeOutput = False # Do we store the output in an output dictionary file?
OPfilename = "OP_Static_"+str(NumSimDays)+'_'+str(samplingBudget/NumSimDays)+'_'+\
            str(diagnosticSensitivity) + '_' + str(diagnosticSpecificity) + '_' + str(globalDemand)
intSFscenario_bool = True # Are we randomly generating some importer SF rates for scenario testing?
endSFscenario_bool = True # Are we randomly generating some outlet SF rates for scenario testing?
saveTestResults = False # Store the testing results to the output dictionary? (Greatly increases file size if True)
'''
testPolicy should be one of: ['Static_Deterministic','Static_Random','Dyn_EpsGreedy',
              'Dyn_EpsExpDecay','Dyn_EpsFirst','Dyn_ThompSamp','Dyn_EveryOther',
              'Dyn_EpsSine','Dyn_TSwithNUTS','Dyn_ExploreWithNUTS',
              'Dyn_ExploreWithNUTS_2','Dyn_ThresholdWithNUTS']
'''
testPolicy = 'Static_Deterministic'
testPolicyParam = [[200,300,400],0.30] # Set testing policy parameter list here

# Diffusion level is set by modifying importer lead times from their suppliers,
# affecting importer stockouts
intLTvar = 0. # Values other than 0. are interpreted as the variance of a log-normal distribution with mean as listed in the imported LT arc list
endLTvar = 0. # Variance in LT for end node procuring from an intermediate node

usePrior = 1. #If not 1., will use regularization instead
optRegularizationWeight = 0.2 # Regularization weight to use with the MLE nonlinear optimizer
lklhdBool = False #Generate the estimates using the likelihood estimator + NUTS (takes time)
lklhdEst_M, lklhdEst_Madapt, lklhdEst_delta = 500, 5000, 0.4 #NUTS parameters

burnInDays_End = 25 # No end-node demand or testing until after this day
burnInDays_Int = 10 # No intermediate demand until after this day

if testPolicy in ['Dyn_TSwithNUTS','Dyn_ExploreWithNUTS']:
    # Sanity checks on the entered parameters for these particular testing policies
    if testPolicyParam[0] > NumSimDays:
        testPolicyParam[0] = NumSimDays # Testing interval can't be more than the total number of days
    if testPolicyParam[0] < 1:
        testPolicyParam[0] = 1 # Minimum planning interval is 1 day
    lklhdBool = True
    testPolicyParam.extend((diagnosticSensitivity,diagnosticSpecificity,lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta))
if testPolicy in ['Dyn_ExploreWithNUTS_2']:
    # Sanity checks for Dyn_ExploreWithNUTS_2
    for elem in testPolicyParam:
        if elem > NumSimDays:
            print('Check the recalculation schedule.')
            break
    lklhdBool = True
    testPolicyParam = [testPolicyParam] + [diagnosticSensitivity,diagnosticSpecificity,lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta]
if testPolicy in ['Dyn_ThresholdWithNUTS']:
    # Sanity checks for Dyn_ThresholdWithNUTS
    for elem in testPolicyParam[0]:
        if elem > NumSimDays:
            print('Check the recalculation schedule.')
            break
    lklhdBool = True
    testPolicyParam = [testPolicyParam[0]] + [testPolicyParam[1]] + [diagnosticSensitivity,diagnosticSpecificity,lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta]




inputParameterDictionary = {'NumSimDays':NumSimDays,'samplingBudget':samplingBudget,
                            'testPolicy':testPolicy,'testPolicyParam':testPolicyParam,
                            'intLTvar':intLTvar,'endLTvar':endLTvar,'globalDemandLevel':globalDemand,
                            'diagnosticSensitivity':diagnosticSensitivity,
                            'diagnosticSpecificity':diagnosticSpecificity,'RglrWt':optRegularizationWeight,
                            'lklhdBool':lklhdBool,'lklhdEst_M':lklhdEst_M,
                            'lklhdEst_Madapt':lklhdEst_Madapt, 'lklhdEst_delta':lklhdEst_delta,
                            'usePrior':usePrior}

useWarmUpFile = False # True if we're going to pull a warm-up file for replications
warmUpRun = False # True if this is a warm-up run to generate bootstrap samples
warmUpIterationGap = 1000 # How often, in sim days, to store the current object lists
# If true, the file name needs to be given here, and its location needs to be in a 'warm up dictionaries' file
# Establish any warm-up settings
numReplications, warmUpRun, useWarmUpFile, warmUpDirectory, warmUpFileName, warmUpDict =\
                        simHelpers.setWarmUp(useWarmUpFileBool = useWarmUpFile,\
                                             warmUpRunBool=warmUpRun, numReps=numReplications,\
                                             currDirect = currDirectory)

##### MAIN SIMULATION CODE #####
# Initialize output collection dictionary
outputDict = {}

# Iterate through each simulation replication
for rep in range(numReplications):

    # Generate the lists of root, intermediate, and end nodes; also keep the list of node headers
    List_RootNode, List_IntermediateNode, List_EndNode, nodeListHeader, nodeList, nodeNum, arcPreferencesMatrix, arcLTsMatrix, arcRsMatrix =\
                                        simHelpers.generateNodeListsFromFile(nodeInputFileString,arcPreferencesFileString,arcLTsFileString,arcRsFileString, NumSimDays, globalDemand)
    rootNum = len(List_RootNode)
    intermediateNum = len(List_IntermediateNode)
    endNum = len(List_EndNode)

    # Simulation statistic lists
    List_RootConsumption = []
    for rootInd in range(len(List_RootNode)):
        List_RootConsumption.append(0)

    # Initialize the drug packet object list
    List_DP = []

    ### TESTING POLICIES GENERATED HERE ###
    # Generate sampling schedule for current graph, as well as a report table shell
    # Depends on whether the sampling plan is static or dynamic
    # Initialize a results list for testing from each end node
    BudgetRemaining = samplingBudget
    List_TestResults = []
    for indEnd in range(endNum):
        List_TestResults.append([rootNum+intermediateNum+indEnd,0,0,1.0,\
                                       np.zeros(intermediateNum,np.int8).tolist()])

    List_TestingSchedule = simTestPolicies.testPolicyHandler(polType=testPolicy,\
                                                             resultsList=List_TestResults,\
                                                             totalSimDays=NumSimDays,\
                                                             numDaysRemain=NumSimDays,\
                                                             totalBudget=samplingBudget,\
                                                             numBudgetRemain=BudgetRemaining,\
                                                             policyParamList=testPolicyParam,\
                                                             startDay=burnInDays_End)
    # Initialize our testing results reporting table
    TestReportTbl = []
    SampleReportTbl = []
    ### THINGS TO CHANGE AND PLAY WITH ONCE THE LISTS ARE GENERATED ###
    # Freeze end node demand to give intermediate nodes a chance to order and stock their shelves
    # We do this by setting the first (maxLT) days at each node equal to 0
    maxLT = 0
    for indInt in range(intermediateNum): # Find the highest lead time at intermediate nodes
        currIntermediate = List_IntermediateNode[indInt]
        if max(currIntermediate.PreferenceLTsList) > maxLT:
            maxLT = max(currIntermediate.PreferenceLTsList)
    for indEnd in range(endNum): # Change demand at the end nodes to zero for these first few days
        currEnd = List_EndNode[indEnd]
        currEnd.demandSched[0:burnInDays_End+1] = 0
    # Also freeze testing for this number of days plus the max number of LT days for end nodes
    List_TestingSchedule = [testRow for testRow in List_TestingSchedule if testRow[0] > burnInDays_End]

    ######### MODIFICATION LOOPS ######### ######### #########

    # Generate an importer SF scenario if the boolean is active
    intSFVec = []
    if intSFscenario_bool == True:
        for indEnd in range(2):
            intSFVec.append(0.01)
        for indEnd in range(2):
            intSFVec.append(0.1)
        for indEnd in range(2):
            intSFVec.append(0.25)
        for indEnd in range(2):
            intSFVec.append(0.5)
        for indEnd in range(1):
            intSFVec.append(0.75)
        for indEnd in range(1):
            intSFVec.append(0.9)
    
    else:
        for indInt in range(intermediateNum):
            intSFVec.append(0)
    random.shuffle(intSFVec)
    for indInt in range(intermediateNum):
        currIntermediate = List_IntermediateNode[indInt]
        currIntermediate.FalsifierProbability = intSFVec[indInt]

    # Do similarly for end node SF scenarios if the boolean is active
    endSFVec = []
    if endSFscenario_bool == True:
        for indEnd in range(76):
            endSFVec.append(0.01)
        for indEnd in range(10):
            endSFVec.append(0.1)
        for indEnd in range(10):
            endSFVec.append(0.25)
        for indEnd in range(5):
            endSFVec.append(0.5)
        for indEnd in range(5):
            endSFVec.append(0.75)
    else:
        for indEnd in range(endNum):
            endSFVec.append(0)
    #print(endSFVec)
    random.shuffle(endSFVec)
    for indEnd in range(endNum):
        currEnd = List_EndNode[indEnd]
        currEnd.FalsifierProbability = endSFVec[indEnd]


    # End nodes - put something inside this loop
    for indEnd in range(endNum):
        currEnd = List_EndNode[indEnd]

    # Intermediate nodes - put something inside this loop
    for indInt in range(intermediateNum):
        currIntermediate = List_IntermediateNode[indInt]


    ######### END MODIFICATION LOOPS ######### ######### #########

    if useWarmUpFile == True: # Randomly select a path from the dictionary
        numScenarios = len(warmUpDict.keys())
        currScenario = int(np.floor(np.random.uniform(0,numScenarios)))
        lastDay = list(warmUpDict[currScenario].keys())[-1] # Retrieve the last day recorded in the scenario's dictionary
        currObjDict = warmUpDict[currScenario][lastDay]
        List_RootNode = currObjDict['rootList']
        List_IntermediateNode, List_EndNode = currObjDict['intList'], currObjDict['endList']
        nodeListHeader, nodeList, nodeNum = currObjDict['nodeHeader'], currObjDict['nodeList'], currObjDict['nodeNumber']
        List_RootConsumption = currObjDict['rootConsumption']
        List_DP = currObjDict['DPList']

    # Initialize simulation time
    startTime = time.time()

    # Main simulation code; iterate over each day
    for today in range(NumSimDays):
        # Intermediate nodes process arriving orders, moving incoming orders remaining days up by one day
        for indInt in range(intermediateNum):
            currIntermediate = List_IntermediateNode[indInt]
            currIntermediate.ProcessIncomingOrder()
        # End nodes process arriving orders, moving incoming orders remaining days up by one day
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            currEnd.ProcessIncomingOrder()

        # End nodes process incoming demand and update consumption records, once burn-in day reached
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            currStockoutDays = currEnd.demandResults[1] # For tracking if today is a stockout day
            List_DP, List_RootNode, List_RootConsumption = currEnd.ProcessScheduleDemand(simDay=today,DPList=List_DP,\
                                                                                         RootList=List_RootNode,\
                                                                                         rootConsumptionVec=List_RootConsumption)
            if currEnd.demandResults[1] > currStockoutDays: # We have a stockout day
                currEnd.DaysStockedOut += 1

        # Run sampling according to sorted schedule
        if List_TestingSchedule != []:
            while List_TestingSchedule[0][0] == today:
                currSample = List_TestingSchedule.pop(0)
                currTestNodeID = currSample[1]
                if nodeList[currTestNodeID][1] == 'FALSE': #Intermediate node
                    endNodeBoolean = False
                else: #End node
                    endNodeBoolean = True
                if not endNodeBoolean:
                    # Find the node in the intermediate list
                    for indInt in List_IntermediateNode:
                        if indInt.id == currTestNodeID:
                            currNode = indInt
                else:
                    # Find the node in the end list
                    for indEnd in List_EndNode:
                        if indEnd.id == currTestNodeID:
                            currNode = indEnd
                # 'currNode' is now the node to be tested
                # Loop through each inventory pile until a sample is found
                madeTest = False #Did we find a sample?
                invPileOriginList = [] # For storing where inventory piles were procured from
                for invPile in range(len(currNode.InventoryDPID)):
                    if (currNode.InventoryPeriodsRemaining[invPile] == 0) and (currNode.InventoryLevel[invPile] > 0): #Test
                        # Do a scan of the inventory present at the tested node, and
                        # record the level of each inventory with respect to the last
                        # node it visited, as long as that node was not node 1 (the falsifier)
                        currPileLevel = currNode.InventoryLevel[invPile] # Greater than 0
                        currPileDPID = currNode.InventoryDPID[invPile]
                        for indDP in List_DP:
                            if currPileDPID == indDP.id:
                                currPilePriorNode = indDP.nodeIDList[-2]
                        invPileOriginList.append([currPilePriorNode,currPileLevel])

                        # If we haven't procured a sample yet, do that here
                        if madeTest == False:
                            madeTest = True
                            testDPID = currNode.InventoryDPID[invPile] #DPID for this inventory pile
                            for indDP in List_DP:
                                if testDPID == indDP.id:
                                    DPRoot = indDP.nodeIDList[0] # Source of this DP
                                    DPImpID = indDP.nodeIDList[1] # Importer of this DP
                                    if DPImpID == currNode.id: # if outlet got this straight from the falsifier, choose an importer based on A
                                        currNodeIndex = currNode.id-intermediateNum-rootNum
                                        if List_TestResults[currNodeIndex][4] !=  np.zeros(intermediateNum,np.int8).tolist():
                                            normedVec = List_TestResults[currNodeIndex][4]/np.sum(List_TestResults[currNodeIndex][4])
                                            DPImpID = rootNum + np.random.choice(range(intermediateNum),p=normedVec)
                                        else: # just pick a random importer from the outlet's preference list
                                            currPrefList = currNode.PreferenceList[:-1]
                                            DPImpID = np.random.choice(currPrefList)
                                            #DPImpID = currPrefList[0]
                                                                                     
                            currNode.InventoryLevel[invPile] -= 1 # "Buy" the sample
                            currNode.demandResults[0] += 1 # "Satisifed" demand
                            # Update available budget
                            BudgetRemaining -= 1 # Reduce our budget balance by 1

                        # END IF
                    # END IF
                if madeTest == False: # Report -1 if no sample procured due to stockout
                    DPRoot = -1
                    List_TestingSchedule.append([today+1,int(np.floor(np.random.uniform(rootNum+intermediateNum,rootNum+intermediateNum+endNum)))]) #append a random node to be tested the next day
                    List_TestingSchedule.sort()
                # END IF

                # Conduct sensitivity and specificity filters
                randUnif = random.uniform(0,1) # Generate a random uniform value
                if DPRoot == 0: # Originally from a non-falsifier
                    if randUnif > diagnosticSpecificity:
                        DPRoot = 1 # Generate a false positive
                elif DPRoot == 1: # Originally from a falsifier
                    if randUnif > diagnosticSensitivity:
                        DPRoot = 0 # Generate a false negative


                # Save the result to the reporting table
                newTestRow = [today,currTestNodeID,DPRoot,invPileOriginList]
                TestReportTbl.append(newTestRow)
                
                # Save the result to the sample-wise result table:
                #   [importer INDEX (not ID), outlet INDEX, result]
                for ind,obj in enumerate(List_EndNode):
                    if obj.id == currTestNodeID:
                        sampOutletID = ind
                for ind,obj in enumerate(List_IntermediateNode):
                    if obj.id == DPImpID:
                        sampImporterID = ind
                SampleReportTbl.append([sampImporterID,sampOutletID,DPRoot])
                
                if madeTest == True: # Recalculate the dynamic results list
                    for currInd in range(len(List_TestResults)):
                        if List_TestResults[currInd][0] == currTestNodeID: # Found matching node ID
                            List_TestResults[currInd][1] += 1 # Add 1 to the number of samples
                            if DPRoot == 1:
                                List_TestResults[currInd][2] += 1
                            List_TestResults[currInd][3] = List_TestResults[currInd][2] / List_TestResults[currInd][1]
                            # Update the inventory levels observed from intermediate nodes
                            for pile in invPileOriginList: # [intermediate node, amount]
                                pileIntNode = pile[0]
                                if not pileIntNode == 1: # Don't count falsified nodes purchases
                                    pileLevel = pile[1]
                                    List_TestResults[currInd][4][pileIntNode-rootNum] += pileLevel

                if List_TestingSchedule == []: # Check if the testing schedule is now empty
                    break

        # End nodes make orders if total inventory is leq the reorder point
        for indEnd in range(endNum):
            if today > burnInDays_Int:
                currEnd = List_EndNode[indEnd]
                List_IntermediateNode, List_DP = currEnd.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP,LeadTimeVariance=endLTvar)

        # Intermediate nodes make orders if total inventory is leq the reorder point
        for indInt in range(intermediateNum):
            currIntermediate = List_IntermediateNode[indInt]
            List_IntermediateNode, List_DP = currIntermediate.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP,LeadTimeVariance=intLTvar)
        # Update dynamic testing if relevant
        if today != NumSimDays-1 and List_TestingSchedule == []:
            List_TestingSchedule = simTestPolicies.testPolicyHandler(polType=testPolicy,\
                                                                 resultsList=List_TestResults,\
                                                                 totalSimDays=NumSimDays,\
                                                                 numDaysRemain=NumSimDays-(today+1),\
                                                                 totalBudget=samplingBudget,\
                                                                 numBudgetRemain=BudgetRemaining,\
                                                                 policyParamList=testPolicyParam)

        if np.mod(today+1,200) == 0: # For updating while running
            print('Rep ' + str(rep+1) + ', Day ' + str(today+1) + ' finished.')

        if today == NumSimDays-1 and np.mod(rep,alertIter)==0: # For updating while running
            simHelpers.CEOTTKbeep()

        if warmUpRun == True and np.mod(today,warmUpIterationGap) == 0: # For storing dictionaries during long run
            warmUpFile = open(warmUpFileName,'rb') # Read the file
            warmUpDict = pickle.load(warmUpFile)
            warmUpFile.close()
            lastKey = len(warmUpDict.keys())-1
            if today==0: # Add a new dictionary for this rep
                lastKey += 1
                warmUpDict[lastKey] = {}
            # Add the current objects, keyed by the simulation day
            currObjDict = {'rootList':List_RootNode,'intList':List_IntermediateNode,
                           'endList':List_EndNode,'nodeHeader':nodeListHeader,
                           'nodeList':nodeList,'nodeNumber':nodeNum,
                           'arcPrefs':arcPreferencesMatrix,'arcLTs':arcLTsMatrix,
                           'arcRs':arcRsMatrix,'rootConsumption':List_RootConsumption,
                           'DPList':List_DP}
            warmUpDict[lastKey][today] = currObjDict
            # Write the dictionary file back
            warmUpFile = open(warmUpFileName,'wb')
            pickle.dump(warmUpDict,warmUpFile)
            warmUpFile.close()



    # END OF SIMULATIONS DAYS LOOP

    # OUTPUT FUNCTIONS
    # Report consumption levels from each node
    # Loop through each root node
    RootReportTbl = []
    Root_Plot_x = []
    Root_Plot_y = []
    rootTotalConsumed = np.sum(List_RootConsumption)
    for indRoot in range(rootNum):
        currRoot = List_RootNode[indRoot]
        currRootConsumption = List_RootConsumption[indRoot]
        currPercConsumption = currRootConsumption / rootTotalConsumed
        currPercString = "{0:.1%}".format(currPercConsumption)
        currRootRow = [currRoot] + [currRootConsumption] + [currPercString]
        RootReportTbl = RootReportTbl + [currRootRow]
        Root_Plot_x.append(str(currRoot))
        Root_Plot_y.append(currPercConsumption)

    IntermediateReportTbl = []
    Intermediate_Plot_x = []
    Intermediate_Plot_y = []
    for indInt in range(intermediateNum):
        currInt = List_IntermediateNode[indInt]
        currInventory = np.sum(currInt.InventoryLevel)
        currTotalDemand = np.sum(currInt.demandResults)
        currPercString = "{0:.1%}".format(currInt.demandResults[1]/currTotalDemand)
        currIntRow = [currInt.id] + currInt.demandResults + [currPercString] + [currInventory]
        IntermediateReportTbl = IntermediateReportTbl + [currIntRow]
        Intermediate_Plot_x.append(str(currInt.id))
        Intermediate_Plot_y.append(currInt.demandResults[1]/currTotalDemand)

    EndReportTbl = []
    End_Plot_x = []
    End_Plot_y = []
    for indEnd in range(endNum):
        currEnd = List_EndNode[indEnd]
        currInventory = np.sum(currEnd.InventoryLevel)
        currTotalDemand = np.sum(currEnd.demandResults)
        currStockoutDays = currEnd.DaysStockedOut
        currPercString = "{0:.1%}".format(currEnd.demandResults[1]/currTotalDemand)
        currEndRow = [currEnd.id] + currEnd.demandResults + [currStockoutDays] + [currPercString] + [currInventory]
        EndReportTbl = EndReportTbl + [currEndRow]
        End_Plot_x.append(str(currEnd.id))
        End_Plot_y.append(currEnd.demandResults[1]/currTotalDemand)

    # Summarize sampling report
    testNodeList = []
    Test_Plot_x = []
    Test_Plot_y1 = []
    Test_Plot_y2 = []
    Test_Plot_y3 = []
    # Initialize all end nodes in testNodeList
    for indEnd in range(endNum):
        currEnd = List_EndNode[indEnd]
        testNodeList.append(currEnd.id)
    
    TestSummaryTbl = []
    InvCheckSummaryTbl = [] # For counting inventory checks during sampling
    tempInvCheckVec = []
    InvCheckHeadersVec = []
    estTransitionMatrix = np.zeros([endNum,intermediateNum]) # For storing transition data
    estFalseVector = np.zeros([endNum,1]) # For storing falsifieds found at each

    for indInt in range(intermediateNum): # Generate headers for printing
        InvCheckHeadersVec.append('N'+str(List_IntermediateNode[indInt].id))
        tempInvCheckVec.append(0)
    for testedNode in testNodeList:
        TestSummaryTbl.append([testedNode,0,0,0]) # Times tested, # falsified, # stocked out
        InvCheckSummaryTbl.append([testedNode] + tempInvCheckVec)
    for sample in TestReportTbl:
        indNodeList = testNodeList.index(sample[1]) # The row of the summary table
        TestSummaryTbl[indNodeList][1] += 1
        if sample[2] == 1: # From falsified root
            TestSummaryTbl[indNodeList][2] += 1
        if sample[2] == -1: # Stocked out; remove 1 from the "samples collected" number
            TestSummaryTbl[indNodeList][3] += 1
            TestSummaryTbl[indNodeList][1] -= 1
        currInvCheckList = sample[3]
        for pile in currInvCheckList:
            if not pile[0] == 1: # Can't be from the falsifier node
                InvCheckSummaryTbl[indNodeList][pile[0]-1] += pile[1]
    j=0
    for testedNode in testNodeList: # For storing the number of falsifieds found at each end node
        if TestSummaryTbl[j][1] > 0:
            percFalse = TestSummaryTbl[j][2] / TestSummaryTbl[j][1]
        else:
            percFalse = 0.
        estFalseVector[j] = percFalse
        j += 1

    i=0
    for row in InvCheckSummaryTbl:
        currSum = np.sum(row[1:])
        for col in range(len(row)):
            if (not col == 0) and (currSum>0):
                row[col] = row[col]/currSum
                estTransitionMatrix[i][col-1] = row[col]
                row[col] = "{0:.1%}".format(row[col])
        TestSummaryTbl[i] = TestSummaryTbl[i] + row[1:]
        i += 1
    TestOverallSummary = [0,0,0]
    for testedNodeRow in TestSummaryTbl:
        TestOverallSummary[0] += testedNodeRow[1] # Total tests
        TestOverallSummary[1] += testedNodeRow[2] # Total falsifieds
        TestOverallSummary[2] += testedNodeRow[3] # Total stock-outs

    counter = 0
    for testedNode in testNodeList: # for plots
        Test_Plot_x.append(testedNode)
        Test_Plot_y1.append(TestSummaryTbl[counter][1]) # Times tested
        Test_Plot_y2.append(TestSummaryTbl[counter][2]) # Times found falsified
        Test_Plot_y3.append(TestSummaryTbl[counter][3]) # Times found stocked out
        counter += 1

        # END OF RX REPORTING LOOP

    # Form regression estimates of suspected bad intermediate nodes
    X = estFalseVector
    A = estTransitionMatrix

    # Get required arguments from the testing summary table
    ydata = []
    numSamples = []
    for endNodeTestRow in TestSummaryTbl:
        ydata.append(endNodeTestRow[2])
        numSamples.append(endNodeTestRow[1])
    ydata = np.array(ydata)
    numSamples = np.array(numSamples)
    # These matrices are for our tracked estimation methods
    Nmat, Ymat = simHelpers.GenerateMatrixForTracked(SampleReportTbl,intermediateNum,endNum)

    #LINEAR PROJECTION
    try:
        output_Lin = simEstMethods.Est_LinearProjection(\
                        A,PosData=ydata,NumSamples=numSamples,\
                        Sens=diagnosticSensitivity,\
                        Spec=diagnosticSpecificity)
        estIntFalsePercList = output_Lin['intProj']
        estEndFalsePercList = output_Lin['endProj']
    except:
        print("Couldn't generate the LINEAR PROJECTION estimates")
        output_Lin = {}
        estIntFalsePercList = []
        estEndFalsePercList = []

    #BERNOULLI MLE
    try:
        output_Bern = simEstMethods.Est_BernoulliProjection(\
                        A,PosData=ydata,NumSamples=numSamples,\
                        Sens=diagnosticSensitivity,\
                        Spec=diagnosticSpecificity)
        estIntFalsePercList_Bern = output_Bern['intProj']
        estEndFalsePercList_Bern = output_Bern['endProj']
    except:
        #print("Couldn't generate the BERNOULLI MLE estimates")
        output_Bern = {}
        estIntFalsePercList_Bern = []
        estEndFalsePercList_Bern = []

    #UNTRACKED POSTERIOR SAMPLING
    if lklhdBool == True:
        try:
            estFalsePerc_PostSampsUNTRACKED = simEstMethods.GeneratePostSamps_UNTRACKED(numSamples,ydata,A,diagnosticSensitivity,\
                                                                       diagnosticSpecificity,optRegularizationWeight,\
                                                                       lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta,usePrior)
                       
        except:
            print("Couldn't generate the UNTRACKED POSTERIOR SAMPLES")
            estFalsePerc_PostSampsUNTRACKED = []
            
    else:
        estFalsePerc_PostSampsUNTRACKED = []
    ### END UNTRACKED POSTERIOR SAMPLING ###
    
    #TRACKED POSTERIOR SAMPLING
    if lklhdBool == True:
        try:
            estFalsePerc_PostSampsTRACKED = simEstMethods.GeneratePostSamps_TRACKED(Nmat,Ymat,diagnosticSensitivity,\
                                                                       diagnosticSpecificity,optRegularizationWeight,\
                                                                       lklhdEst_M,lklhdEst_Madapt,lklhdEst_delta,usePrior)
           
        except:
            print("Couldn't generate the TRACKED POSTERIOR SAMPLES")
            estFalsePerc_PostSampsTRACKED = []
            
    else:
        estFalsePerc_PostSampsTRACKED = []
    ### END TRACKED POSTERIOR SAMPLING ###
    
    #UNTRACKED MLE
    try:
        k = 10
        if len(estFalsePerc_PostSampsUNTRACKED)>0: #Use a random selection of the posterior samples for the MLE maximization if available
            randInds = random.sample(range(len(estFalsePerc_PostSampsUNTRACKED)),min(k,len(estFalsePerc_PostSampsUNTRACKED)))
            randList = estFalsePerc_PostSampsUNTRACKED[randInds]
        else:
            randList = []
            
        output_Untracked = simEstMethods.Est_UntrackedMLE(\
                                   A,PosData=ydata,NumSamples=numSamples,\
                                   Sens=diagnosticSensitivity,\
                                   Spec=diagnosticSpecificity,\
                                   RglrWt=optRegularizationWeight,\
                                   beta0_List=randList,
                                   usePrior=usePrior)
        estIntMLE_Untracked = output_Untracked['intProj']
        estEndMLE_Untracked = output_Untracked['endProj']
    except:
        print("Couldn't generate the UNTRACKED MLE estimates")
        output_Untracked = {}
        estIntMLE_Untracked = []
        estEndMLE_Untracked = []
    ### END UNTRACKED MLE ###
    
    #TRACKED MLE
    try:
        k = 10
        if len(estFalsePerc_PostSampsTRACKED)>0: #Use a random selection of the posterior samples for the MLE maximization if available
            randInds = random.sample(range(len(estFalsePerc_PostSampsTRACKED)),min(k,len(estFalsePerc_PostSampsTRACKED)))
            randList = estFalsePerc_PostSampsTRACKED[randInds]
        else:
            randList = []
            
        output_Tracked = simEstMethods.Est_TrackedMLE(\
                                   N=Nmat,Y=Ymat,\
                                   Sens=diagnosticSensitivity,\
                                   Spec=diagnosticSpecificity,\
                                   RglrWt=optRegularizationWeight,\
                                   beta0_List=randList,
                                   usePrior=usePrior)
        estIntMLE_Tracked = output_Tracked['intProj']
        estEndMLE_Tracked = output_Tracked['endProj']
    except:
        print("Couldn't generate the TRACKED MLE estimates")
        output_Tracked = {}
        estIntMLE_Tracked = []
        estEndMLE_Tracked = []
    ### END TRACKED MLE ###
    
    # Simulation end time
    totalRunTime = time.time() - startTime
    
    if printOutput == True: # Printing of tables and charts
        # PRINT RESULTS TABLES
        print('*'*100)
        print('ROOT NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(RootReportTbl, headers=['Root ID', 'End Consumption', 'Overall Percentage']))
        print('*'*100)
        print('INTERMEDIATE NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(IntermediateReportTbl, headers=['Intermediate Node ID', 'Satisfied Demand', 'Stockout Demand','Stockout %','Inventory Available']))
        print('*'*100)
        print('END NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(EndReportTbl, headers=['End Node ID', 'Satisfied Demand', 'Stockout Demand', 'Stockout Days','Stockout %', 'Inventory Available']))
        print('*'*100)
        print('TESTING SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(TestSummaryTbl, headers=['Node ID', '# of Samples', '# SFP', '# Stocked Out']+InvCheckHeadersVec))
        print('*'*100)
        print('OVERALL SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate([[totalRunTime]+TestOverallSummary],headers=['Simulation Run Time','Total Tests','SFPs Found','Stocked Out Found']))

        # GENERATE PLOTS
        # Root consumption %'s
        fig = plt.figure()
        ax = fig.add_axes([0,0,0.3,0.5])
        ax.set_xlabel('Root Node',fontsize=16)
        ax.set_ylabel('Percentage consumption',fontsize=16)
        ax.bar(Root_Plot_x,Root_Plot_y)
        vals = ax.get_yticks()
        ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        plt.show()

        # Intermediate stockout %'s
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,0.5])
        ax.set_xlabel('Intermediate Node',fontsize=16)
        ax.set_ylabel('Percentage stocked out',fontsize=16)
        #vals = ax.get_yticks()
        #ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        ax.bar(Intermediate_Plot_x,Intermediate_Plot_y)
        plt.show()

        # End stockout %'s
        fig = plt.figure()
        ax = fig.add_axes([0,0,3,0.5])
        ax.set_xlabel('End Node',fontsize=16)
        ax.set_ylabel('Percentage stocked out',fontsize=16)
        #vals = ax.get_yticks()
        #ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        ax.bar(End_Plot_x,End_Plot_y)
        plt.show()

        # Testing results
        Test_Plot_x = np.array(Test_Plot_x)
        fig = plt.figure()
        ax = fig.add_axes([0,0,3,0.5])
        ax.set_xlabel('Tested Node',fontsize=16)
        ax.set_ylabel('Result Amount',fontsize=16)
        ax.bar(Test_Plot_x+0.00, Test_Plot_y1, color='black', width=0.25)
        ax.bar(Test_Plot_x+0.25, Test_Plot_y2, color='red', width=0.25)
        ax.bar(Test_Plot_x+0.50, Test_Plot_y3, color='gray', width=0.25)
        plt.legend(('Times tested','Times falsified','Times stocked out'))
        plt.xticks(rotation=90)
        plt.show()

        # Intermediate nodes
        fig, axs = plt.subplots(3, 1,figsize=(9,13))
        fig.suptitle('Intermediate Node SF % Estimates',fontsize=18)
        for subP in range(3):
            axs[subP].set_xlabel('Intermediate Node',fontsize=12)
            axs[subP].set_ylabel('Est. SF %',fontsize=12)
            axs[subP].set_ylim([0.0,1.0])
            vals = axs[subP].get_yticks()
            axs[subP].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        # Linear projection
        axs[0].set_title('Linear projection',fontweight='bold')
        axs[0].bar(Intermediate_Plot_x,estIntFalsePercList,color='thistle',edgecolor='indigo')
        # Bernoulli MLE projection
        #axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
        #axs[1].bar(Intermediate_Plot_x,estIntFalsePercList_Bern,color='mediumspringgreen',edgecolor='green')
        # MLE w Nonlinear optimizer
        axs[1].set_title('Untracked MLE',fontweight='bold')
        axs[1].bar(Intermediate_Plot_x,estIntMLE_Untracked,color='navajowhite',edgecolor='darkorange')
        # Sample MLE
        axs[2].set_title('Tracked MLE',fontweight='bold')
        axs[2].bar(Intermediate_Plot_x,estIntMLE_Tracked,color='deepskyblue',edgecolor='darkcyan')
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.show()

        # End nodes
        fig, axs = plt.subplots(3, 1,figsize=(17,13))
        fig.suptitle('End Node SF % Estimates',fontsize=18)
        for subP in range(3):
            axs[subP].set_xlabel('End Node',fontsize=12)
            axs[subP].set_ylabel('Est. SF %',fontsize=12)
            axs[subP].tick_params(labelrotation=90)
            axs[subP].set_ylim([-0.6,0.6])
            vals = axs[subP].get_yticks()
            axs[subP].set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        # Linear projection
        axs[0].set_title('Linear projection',fontweight='bold')
        axs[0].bar(End_Plot_x,estEndFalsePercList,color='peachpuff',edgecolor='red')
        # Bernoulli MLE projection
        #axs[1].set_title('Bernoulli MLE projection',fontweight='bold')
        #axs[1].bar(End_Plot_x,estEndFalsePercList_Bern,color='khaki',edgecolor='goldenrod')
        # Untracked MLE
        axs[1].set_title('Untracked MLE',fontweight='bold')
        axs[1].bar(End_Plot_x,estEndMLE_Untracked,color='lightcyan',edgecolor='teal')
        # Tracked MLE
        axs[2].set_title('Tracked MLE',fontweight='bold')
        axs[2].bar(End_Plot_x,estEndMLE_Tracked,color='pink',edgecolor='deeppink')
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.show()

        #Posterior samples
        if lklhdBool == True:
            #UNTRACKED POSTERIOR SAMPLES
            fig = plt.figure()
            ax = fig.add_axes([0,0,2,1])
            ax.set_xlabel('Intermediate Node',fontsize=16)
            ax.set_ylabel('UNTRACKED posterior distribution',fontsize=16)
            for i in range(intermediateNum):
                plt.hist(sps.expit(estFalsePerc_PostSampsUNTRACKED[:,i]))

            fig = plt.figure()
            ax = fig.add_axes([0,0,2,1])
            ax.set_xlabel('End Node',fontsize=16)
            ax.set_ylabel('UNTRACKED posterior distribution',fontsize=16)
            for i in range(endNum):
                plt.hist(sps.expit(estFalsePerc_PostSampsUNTRACKED[:,intermediateNum+i]))

            meanSampVec = []
            for i in range(intermediateNum):
                meanSampVec.append(np.mean(sps.expit(estFalsePerc_PostSampsUNTRACKED[:,i])))
            meanSampVec = [round(meanSampVec[i],3) for i in range(len(meanSampVec))]
            
            #TRACKED POSTERIOR SAMPLES
            fig = plt.figure()
            ax = fig.add_axes([0,0,2,1])
            ax.set_xlabel('Intermediate Node',fontsize=16)
            ax.set_ylabel('TRACKED posterior distribution',fontsize=16)
            for i in range(intermediateNum):
                plt.hist(sps.expit(estFalsePerc_PostSampsTRACKED[:,i]))

            fig = plt.figure()
            ax = fig.add_axes([0,0,2,1])
            ax.set_xlabel('End Node',fontsize=16)
            ax.set_ylabel('TRACKED posterior distribution',fontsize=16)
            for i in range(endNum):
                plt.hist(sps.expit(estFalsePerc_PostSampsTRACKED[:,intermediateNum+i]))

            meanSampVec = []
            for i in range(intermediateNum):
                meanSampVec.append(np.mean(sps.expit(estFalsePerc_PostSampsTRACKED[:,i])))
            meanSampVec = [round(meanSampVec[i],3) for i in range(len(meanSampVec))]
            
            
            

    ### END OF PRINT OUTPUT LOOP

    if warmUpRun == False:
        # Update our output dictionary
        List_demandResultsInt = [] # For intermediate node demand results
        for indInt in range(intermediateNum):
            currInt = List_IntermediateNode[indInt]
            List_demandResultsInt.append(currInt.demandResults)
        List_demandResultsEnd = [] # For end node demand results
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            List_demandResultsEnd.append(currEnd.demandResults)

        # Update the "true" underlying vector for end nodes to be the underlying rate + the proportion of orders due to stockouts
        endSFVecCombo = []
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            if currEnd.AmountProcured > 0:
                currStockoutOrderRate = currEnd.AmountProcuredDueToStockouts / currEnd.AmountProcured
            else:
                currStockoutOrderRate = 0
            endSFVecCombo.append(endSFVec[indEnd] + currStockoutOrderRate)

        #Send the testing results to the output dictionary?

        if saveTestResults == True:
            TestReportTblToSend = TestReportTbl
        else:
            TestReportTblToSend = []

        currOutputLine = {'inputParameterDictionary':inputParameterDictionary,
                          'rootConsumption':List_RootConsumption,
                          'intDemandResults':List_demandResultsInt,
                          'endDemandResults':List_demandResultsEnd,
                          'testResults':TestReportTblToSend,
                          'dynTestResults':List_TestResults,
                          'intFalseEstimates':estIntFalsePercList,
                          'endFalseEstimates':estEndFalsePercList,
                          'intFalseEstimates_Bern':estIntFalsePercList_Bern,
                          'endFalseEstimates_Bern':estEndFalsePercList_Bern,
                          'intEstMLE_Untracked':estIntMLE_Untracked,
                          'endEstMLE_Untracked':estEndMLE_Untracked,
                          'intEstMLE_Tracked':estIntMLE_Tracked,
                          'endEstMLE_Tracked':estEndMLE_Tracked,
                          'postSamps_Untracked':estFalsePerc_PostSampsUNTRACKED,
                          'postSamps_Tracked':estFalsePerc_PostSampsTRACKED,
                          'estDict_Lin':output_Lin, 'estDict_Bern':output_Bern,
                          'estDict_Untracked':output_Untracked,'estDict_Tracked':output_Tracked,
                          'intSFTrueValues':intSFVec,'endSFTrueValues':endSFVecCombo,
                          'simStartTime':startTime,
                          'simRunTime':totalRunTime
                          }

        outputDict[rep] = currOutputLine # Save to the output dictionary
    
########## END OF REPLICATION LOOP ##########

# Store the outputDict
if warmUpRun == False and storeOutput == True:
    outputFilePath  = os.getcwd() + '\\outputDictionaries'
    if not os.path.exists(outputFilePath): # Generate this folder if one does not already exist
            os.makedirs(outputFilePath)
    #outputFileName = os.path.basename(sys.argv[0])[:-3] + '_OUTPUT' # Current file name
    outputFileName = os.path.join(outputFilePath, OPfilename)
    pickle.dump(outputDict, open(outputFileName,'wb'))


# WRITE TO CSV FILE
#with open('testRxResults.csv', 'w', newline='') as file:
#    writer = csv.writer(file)
#    writer.writerow(["RxID", "GoodConsumed", "SubstandardConsumed", "FalsifiedConsumed", "StockoutConsumed"]) # Initialize the header
#    for rw in range(NumPrivateRxs):
#        writer.writerow(RxReportTbl[rw])


### WISH LIST:
# 1) ADD WARM-UP PERIOD
        # Have multiple (100?) long warm-up periods that we sample from each replication
# 2) PUT A SMALLER FALSIFIER "ORDER" AMOUNT IN THE "ROOT" SECTION OF 'MAKEORDER'
#       TRIGGERED BY CURRENT SUPPLIER NODE HAVING THE 'FALSIFIER' LABEL
#       ALSO POTENTIALLY TRIGGERED BY NUMBER OF DAYS STOCKED OUT
#       ORDER AMOUNT IS REORDER POINT/2
# 2.1) RECORD NUMBER OF TRIGGERS
# 3) RUN MULTIPLE REPLICATIONS OVER LONG-RUN SIMULATIONS
# 5) OUTPUT BATCH CONSUMPTION RATE
# 7) "SANDY" CHECKS TO ENSURE THINGS ARE RUNNING SMOOTHLY
# 8) PUT WARM-UP LINES INTO A MODULE
# 9) BUILD INVENTORY CHECK SUMMARIES TO HANDLE INTERMEDIATE NODES



