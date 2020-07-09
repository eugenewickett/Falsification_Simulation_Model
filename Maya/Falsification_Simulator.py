# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:59:35 2019

@author: Eugene Wickett

This code is for building and running simulations of a supply chain susceptible
to falsification and substandardization.
"""


import matplotlib.pyplot as plt
import numpy as np
import random #for seeds
import sys
import csv # for csv manipulation
import time # for time tracking
import os # for directories
from tabulate import tabulate # for making outputs
import pickle # for saving/loading objects in Python
#import winsound # for beeps

import Falsification_Sim_Classes as simClasses # our class objects and methods
import Falsification_Sim_Modules as simModules # modules for the simulation

# Run supporting files
currDirectory = os.path.dirname(os.path.realpath(__file__)) #Run this command to get the current working directory string
os.chdir(currDirectory) # Set directory    

##### GRAPH INPUT FILES
# Read input lists: node list, arc preferences, and arc lead times
nodeInputFileString = 'LIB_Nodes_1.csv' #says whether root, int, or end node & reorder points
arcPreferencesFileString = 'LIB_Arcs_Preferences_1.csv' #where nodes prefer to get stock from
arcLTsFileString = 'LIB_Arcs_LTs_1.csv' # number days from outlet to importer
arcRsFileString = 'LIB_Arcs_Rs_1.csv'#how much each node orders from prev node

##### SIMULATION PARAMETERS
# Enter the length of the simulation and the sampling budget
NumSimDays = 400
samplingBudget = NumSimDays*5
numReplications = 1
testPolicy = 2
printOutput = True # Whether individual replication output should be displayed
diagnosticSensitivity = 0.95 # Tool sensitivity
diagnosticSpecificity = 0.98 # Tool specificity
useWarmUpFile = False # True if we're going to pull a warm-up file for replications
warmUpRun = False # True if this is a warm-up run to generate bootstrap samples
warmUpIterationGap = 1000 # How often, in sim days, to store the current object lists
# If true, the file name needs to be given here, and its location needs to be in a 'warm up dictionaries' file
storeOutput = False # Do we store the output in an output dictionary file?
alertIter = 5 # How frequently we're alerted of a set of replications being completed

# Establish any warm-up settings 
numReplications, warmUpRun, useWarmUpFile, warmUpDirectory, warmUpFileName, warmUpDict = simModules.setWarmUp(useWarmUpFileBool = useWarmUpFile, warmUpRunBool = warmUpRun, numReps = numReplications, currDirect = currDirectory)

##### MAIN SIMULATION CODE #####
# Initialize output collection dictionary
outputDict = {}

# Iterate through each simulation replication
for rep in range(numReplications):

    # Generate the lists of root, intermediate, and end nodes; also keep the list of node headers
    List_RootNode, List_IntermediateNode, List_EndNode, nodeListHeader, nodeList, nodeNum, arcPreferencesMatrix, arcLTsMatrix, arcRsMatrix = simModules.generateNodeListsFromFile(nodeInputFileString,arcPreferencesFileString,arcLTsFileString,arcRsFileString, NumSimDays)
    rootNum = len(List_RootNode)
    intermediateNum = len(List_IntermediateNode)
    endNum = len(List_EndNode)
    
    # Simulation statistic lists
    List_RootConsumption = []
    for rootInd in range(len(List_RootNode)):
        List_RootConsumption.append(0) #MZ: list of zeros, size of list_rootnode
    
    # Initialize the drug packet object list
    List_DP = []
    
    ### TESTING POLICIES GENERATED HERE ###
    # Generate sampling schedule for current graph, as well as a report table shell
    List_TestingSchedule = simModules.testingScheduleGenerator(nodes=nodeList, int_numDays=NumSimDays, int_sampleBudget = samplingBudget, int_PolicyType = testPolicy, arr_PolicyParameter = [0])
    TestReportTbl = []
    
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
        currEnd.demandSched[0:(maxLT+2)] = 0
    # Also freeze testing for this number of days plus the max number of LT days for end nodes
    List_TestingSchedule = [testRow for testRow in List_TestingSchedule if testRow[0] > (maxLT*7)]
    
    
    ######### MODIFICATION LOOPS ######### ######### #########
    
    # Intermediate nodes - put something inside this loop
    for indInt in range(intermediateNum):
        currIntermediate = List_IntermediateNode[indInt]
        
    # End nodes - put something inside this loop
    for indEnd in range(endNum):
        currEnd = List_EndNode[indEnd]
    
    #List_IntermediateNode[0].FalsifierProbability = 0.2
    
    
    ######### MODIFICATION LOOPS ######### ######### #########
    
    
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
    testid_lst  = []  #MZL
    timestestedsum = 0 #MZL
    list_prev_id = []  #MZ
    parent_list =[] #MZ
    unique_par_list  =[] #MZ
    # Main simulation code; iterate over each day
    
    for today in range(NumSimDays):
        # Intermediate nodes process arriving orders, moving incoming orders remaining days up by one day
        for indInt in range(intermediateNum):
            currIntermediate = List_IntermediateNode[indInt]
            currIntermediate.ProcessIncomingOrder()
            #print(currIntermediate.prev_id_lst)  #MZ
        # End nodes process arriving orders, moving incoming orders remaining days up by one day
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            #print(currEnd.prev_id_int, currEnd.id) #MZ

            currEnd.ProcessIncomingOrder()
            ####MZ START
            prev_id =  currEnd.prev_id_int
            list_prev_id.append([prev_id, currEnd.id])
            #mz end 
        # End nodes process incoming demand and update consumption records    
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            currStockoutDays = currEnd.demandResults[1] # For tracking if today is a stockout day
            List_DP, List_RootNode, List_RootConsumption = currEnd.ProcessScheduleDemand(simDay=today,DPList=List_DP, RootList=List_RootNode, rootConsumptionVec=List_RootConsumption)
            if currEnd.demandResults[1] > currStockoutDays: # We have a stockout day
                currEnd.DaysStockedOut += 1
            
        # Run sampling according to sorted schedule
        
        if List_TestingSchedule != []:
            while List_TestingSchedule[0][0] == today:
                currSample = List_TestingSchedule.pop(0)
                currTestNodeID = currSample[1]
                testid_lst.append(currTestNodeID) #MZ
                
                #MZ START
                if testid_lst.count(currTestNodeID) != 0:
                    timestested = testid_lst.count(currTestNodeID)
                    timestestedsum += timestested
                
                elif testid_lst.count(currTestNodeID) == 0:
                    timestested = 0
                    timestestedsum += timestested
                    print("times tested is zero?")
                timestestedavg=timestestedsum/len(testid_lst)
                #print(timestested, timestestedsum, len(testid_lst), timestestedavg)
                #print(testid_lst)
                #MZEND
                List_Root_and_IntermediateNode =  List_RootNode+List_IntermediateNode #MZ
                if nodeList[currTestNodeID][1] == 'False':   
                    endNodeBoolean = False
                    rootNodeBoolean = False
                elif nodeList[currTestNodeID][0]=='False':
                    rootNodeBoolean = True
                else: #End node
                    endNodeBoolean = True
                    rootNodeBoolean = False
                if not endNodeBoolean:
                    
                    if rootNodeBoolean:  # Find the node in the root list
                        for rootInt in List_RootNode:
                            if rootInt == currTestNodeID:
                                currNode = rootInt

                    else:  # Find the node in the intermediate list
                        for indInt in List_Root_and_IntermediateNode:
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
                if not rootNodeBoolean:
                    for invPile in range(len(currNode.InventoryDPID)):
                        if (currNode.InventoryPeriodsRemaining[invPile] == 0) and (currNode.InventoryLevel[invPile] > 0): #Test
                        # Do a scan of the inventory present at the tested node, and
                        # record the level of each inventory with respect to the last
                        # node it visited, as long as that node was not node 1 (the falsifier)
                            currPileLevel = currNode.InventoryLevel[invPile] # Greater than 0
                            currPileDPID = currNode.InventoryDPID[invPile]
                        
                            for indDP in List_DP: #MZ maybe helpful
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
                                        #MZ start
                                        DPInt = indDP.nodeIDList[1]
                                        DP_parent = indDP.intNode
                                        parent_list.append(DP_parent)
                                        #MZ end

                                    #MZ get the prev parent DPID of this ID. the drug packet thats found needs to have a parent stored
                                    #MZ add DrugPacket property. if its a child retain who was the last id. 
                                    #MZ result of test should have that parents drug pacekt id in the testing report table
                                    #print(DPInt)
                                    #print(indDP.intNode)
                                    #MZ end
                                currNode.InventoryLevel[invPile] -= 1 # "Buy" the sample
                                currNode.demandResults[0] += 1 # "Satisifed" demand
                        # END IF
                    # END IF
                    if madeTest == False: # Report -1 if no sample procured due to stockout
                        DPRoot = -1
                    # END IF     

                    # Conduct sensitivity and specificity filters
                    randUnif = random.uniform(0,1) # Generate a random uniform value
                    if DPRoot == 0: # Originally from a non-falsifier
                        if randUnif > diagnosticSpecificity:
                            DPRoot = 1 # Generate a false positive
                    elif DPRoot == 1: # Originally from a falsifier
                        DPRoot = 0 # Generate a false negative
                        
                
                    # Save the result to the reporting table
                    newTestRow = [today,currTestNodeID,DPRoot,invPileOriginList, DP_parent] #updated
                    TestReportTbl.append(newTestRow) 
                    #MZ in output make list of how many tests i've had, how many unique batch #s have I seen, what is the distribution

                #print(TestReportTbl, parent_list)
                    if List_TestingSchedule == []: # Check if the testing schedule is now empty
                        break
               
        # End nodes make orders if total inventory is leq the reorder point
        for indEnd in range(endNum):
            currEnd = List_EndNode[indEnd]
            List_IntermediateNode, List_DP = currEnd.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP)
        
        # Intermediate nodes make orders if total inventory is leq the reorder point
        for indInt in range(intermediateNum):
            currIntermediate = List_IntermediateNode[indInt]
            List_IntermediateNode, List_DP = currIntermediate.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP)              
        
        if np.mod(today+1,200) == 0: # For updating while running
            print('Rep ' + str(rep+1) + ', Day ' + str(today+1) + ' finished.')
        
        if today == NumSimDays-1 and np.mod(rep,alertIter)==0: # For updating while running
            pass
            #winsound.Beep(int(32.7032 * (2**3)*(1.059463094**10)),400)
            #winsound.Beep(int(32.7032 * (2**3)*(1.059463094**12)),400)
            #winsound.Beep(int(32.7032 * (2**3)*(1.059463094**8)),400)
            #winsound.Beep(int(32.7032 * (2**2)*(1.059463094**8)),400)
            #winsound.Beep(int(32.7032 * (2**3)*(1.059463094**3)),400)            
        
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
    # Simulation end time
    totalRunTime = [time.time() - startTime]
     
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
        currInvPeriod  = np.sum(currInt.InventoryPeriodsRemaining) #MZ added this line
        currTotalDemand = np.sum(currInt.demandResults)
        currPercString = "{0:.1%}".format(currInt.demandResults[1]/currTotalDemand)
        currIntRow = [currInt.id] + currInt.demandResults + [currPercString] + [currInventory] + [currInvPeriod]
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
        currInvPeriod = np.sum(currEnd.InventoryPeriodsRemaining)
        currEndRow = [currEnd.id] + currEnd.demandResults + [currStockoutDays] + [currPercString] + [currInventory] +  [currInvPeriod]
        EndReportTbl = EndReportTbl + [currEndRow]
        End_Plot_x.append(str(currEnd.id))
        End_Plot_y.append(currEnd.demandResults[1]/currTotalDemand)
    
    # Summarize sampling report
    testNodeList = []
    Test_Plot_x = []
    Test_Plot_y1 = []
    Test_Plot_y2 = []
    Test_Plot_y3 = []
    for sample in TestReportTbl:
        if not sample[1] in testNodeList:
            testNodeList.append(sample[1])
    testNodeList.sort()
    TestSummaryTbl = []
    
    InvCheckSummaryTbl = [] # For counting inventory checks during sampling
    tempInvCheckVec = []
    InvCheckHeadersVec = []
    estTransitionMatrix = np.zeros([endNum,intermediateNum]) # For storing transition data 
    estFalseVector = np.zeros([endNum,1]) # For storing falsifieds found at each 
    
    for indInt in range(intermediateNum):
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
        if sample[2] == -1: # Stocked out
            TestSummaryTbl[indNodeList][3] += 1
        currInvCheckList = sample[3]
        for pile in currInvCheckList:
            if not pile[0] == 1: # Can't be from the falsifier node
                InvCheckSummaryTbl[indNodeList][pile[0]-1] += pile[1]
    j=0            
    for testedNode in testNodeList: # For storing the number of falsifieds found at each end node
        percFalse = TestSummaryTbl[j][2] / TestSummaryTbl[j][1]
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
    try:
        estIntFalsePerc = np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,X))
        estEndFalsePerc = np.subtract(X,np.dot(A,estIntFalsePerc))
    
        # Estimated intermediate falsification percentages
        # First, some stats for the plots
        estIntFalsePercList = np.ndarray.tolist(estIntFalsePerc.T[0])
        estIntFalsePerc_median = np.median(estIntFalsePercList)
        estIntFalsePerc_SD1 = np.std(estIntFalsePercList)
        estIntFalsePerc_SD2 = 2*estIntFalsePerc_SD1
        # Store these values in one place
        plot_y = (estIntFalsePerc_median,estIntFalsePerc_median+estIntFalsePerc_SD1,estIntFalsePerc_median+estIntFalsePerc_SD2)
        # For end nodes
        estEndFalsePercList = np.ndarray.tolist(estEndFalsePerc.T[0])
    except:
        print("Couldn't generate the estimated node falsification percentages")
        estIntFalsePercList = []
        estEndFalsePercList = []
    
    
    if printOutput == True: # Printing of tables and charts
        # PRINT RESULTS TABLES
        print('*'*100)
        print('ROOT NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(RootReportTbl, headers=['Root ID', 'End Consumption', 'Overall Percentage']))
        print('*'*100)
        print('INTERMEDIATE NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(IntermediateReportTbl, headers=['Intermediate Node ID', 'Satisfied Demand', 'Stockout Demand','Stockout %','Inventory Available', "Inventory Period"]))
        print('*'*100)
        print('END NODE SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(EndReportTbl, headers=['End Node ID', 'Satisfied Demand', 'Stockout Demand', 'Stockout Days','Stockout %', 'Inventory Available', 'Inventory Period']))
        print('*'*100)
        print('TESTING SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate(TestSummaryTbl, headers=['Node ID', '# of Samples', '# Falsified', '# Stocked Out']+InvCheckHeadersVec))
        print('*'*100)
        print('OVERALL SUMMARY STATISTICS')
        print('*'*100)
        print(tabulate([totalRunTime+TestOverallSummary],headers=['Simulation Run Time','Total Tests','Falsifieds Found','Stocked Out Found']))
        #MZ added
        print('Amount of tests performed: ' + str(len(parent_list)))
        
        for p in parent_list:
            if unique_par_list.count(p)==0:
                unique_par_list.append(p)
        len_unique_par_list = len(unique_par_list)
        print('Amount of unique batch numbers seen: ' + str(len_unique_par_list))



        # GENERATE PLOTS
        # Root consumption %'s
        fig = plt.figure()
        ax = fig.add_axes([0,0,0.3,0.5])
        ax.set_xlabel('Root Node',fontsize=16)
        ax.set_ylabel('Percentage consumption',fontsize=16)
        vals = ax.get_yticks()
        ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        ax.bar(Root_Plot_x,Root_Plot_y)
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
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,0.5])
        ax.set_xlabel('Intermediate Node',fontsize=16)
        ax.set_ylabel('Estimated falsification percentage',fontsize=16)
        ax.bar(Intermediate_Plot_x,estIntFalsePercList,color='thistle',edgecolor='indigo')
        plt.hlines(y=plot_y,xmin=0,xmax=(intermediateNum-1),colors=("orangered"),linestyles=['solid','dashed','dotted'])
        plt.show()
             
        # End nodes
        fig = plt.figure()
        ax = fig.add_axes([0,0,3,0.5])
        ax.set_xlabel('End Node',fontsize=16)
        ax.set_ylabel('Estimated falsification percentage',fontsize=16)
        #vals = ax.get_yticks()
        #ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
        ax.bar(End_Plot_x,estEndFalsePercList,color='peachpuff',edgecolor='red')
        plt.xticks(rotation=90)
        plt.show()
        ##MZ start      
        freq  = [0] * len_unique_par_list
        unique_par_list.sort()  
        print(unique_par_list)
        for p1 in range(len_unique_par_list):
            freq[p1] = parent_list.count(unique_par_list[p1]) 
        print(freq)
        plt.plot(unique_par_list, freq, '.')
        plt.xlabel('Parent Nodes')
        plt.ylabel('Frequency')
        plt.show()
        ##MZ end   
        """
        ### MZ HISTOGRAM#
        x_values = []
        y_values = []
        str_list_prev_id =[]
        for p in list_prev_id:
            x1 = p[0] #MZ prev_id value of DP
            y1 = p[1] #MZ Current id value of DP (the end node id of the current DP)
            x_values.append(x1)
            y_values.append(y1)
        plt.hist(x_values, bins = 100, color = 'pink')
        plt.title('Intermediate Node Frequency')
        plt.xlabel('Int Nodes')
        plt.ylabel('Frequency')
        plt.show()
        plt.hist(y_values, bins = 150,  color = 'blue')
        plt.title('End Node Frequency')
        plt.xlabel('End Nodes')
        plt.ylabel('Frequency')
        plt.show()
        plt.plot(x_values,y_values, '.', markersize=5)
        plt.xlabel('Intermediate Nodes')
        plt.ylabel('End Nodes')
        plt.title('End vs Int Nodes')
        plt.show()
        
        for sublist in list_prev_id:
            subliststr= str(sublist)
            str_list_prev_id.append(subliststr)
        #print(len(str_list_prev_id))    
        plt.hist(str_list_prev_id, bins=5000, color ='red')
        plt.title('Transistions from Previous Intermediate Node to End Node')
        #plt.rcParams.update({'font.size': 2})
        plt.xlabel('[Intermediate Node, End Node]')
        plt.ylabel('Frequency')
        plt.show()
    ##MZ END
    """
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
        
        currOutputLine = {'rootConsumption':List_RootConsumption,
                          'intDemandResults':List_demandResultsInt,
                          'endDemandResults':List_demandResultsEnd,
                          'testResults':TestReportTbl,
                          'intFalseEstimates':estIntFalsePercList,
                          'endFalseEstimates':estEndFalsePercList}
                          
        outputDict[rep] = currOutputLine # Save to the output dictionary
        
        
########## END OF REPLICATION LOOP ##########

# Store the outputDict
if warmUpRun == False and storeOutput == True:
    outputFilePath  = os.getcwd() + '\\output dictionaries'
    if not os.path.exists(outputFilePath): # Generate this folder if one does not already exist
            os.makedirs(outputFilePath)
    outputFileName = os.path.basename(sys.argv[0])[:-3] + '_OUTPUT' # Current file name
    outputFileName = os.path.join(outputFilePath, outputFileName)
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
# 4) GENERATE BATCH FILES VARYING DIFFERENT PARAMETERS
# 5) OUTPUT BATCH CONSUMPTION RATE
# 6) STREAMLINE TESTING/SAMPLING PROCESS INTO A MODULE
# 7) "SANDY" CHECKS TO ENSURE THINGS ARE RUNNING SMOOTHLY
# 8) PUT WARM-UP LINES INTO A MODULE
# 9) BUILD INVENTORY CHECK SUMMARIES TO HANDLE INTERMEDIATE NODES



