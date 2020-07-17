"""

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
NumSimDays = 500
samplingBudget = 1000
numReplications = 1
testPolicy = 3
epsilon =  .1 #MZ. relevant only for testpolicy 3
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
dynamic_list = [0]*106
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
    TestReportTbl = []
    if testPolicy  == 3: 
        dynamicbool =True
        #initialize list of falsification rates as every end node has 100% falsification rates until proven otherwise
        bud_rem = samplingBudget #initialize budget remaining to samplingBudget
        days_rem = NumSimDays #initialize sampling days remaining to NumSimDays
        bud_over_days = bud_rem/days_rem
        remainder_bud_over_days = int(bud_rem%days_rem) #get the remainder. How many days in the simulation can afford 1 extra tests / day 
        floor_bud_over_days = int(bud_rem//days_rem) #get the regular tests / day 
        List_tests_per_day = (([floor_bud_over_days]*(NumSimDays-remainder_bud_over_days) + [floor_bud_over_days+1]*remainder_bud_over_days)) 
        Epsilon =  .1
        List_Report  = [[0,0,1] for x in range(106)] 
        List_False_Rates = [1]*106
        List_Nodes_Tested=[]
        sampleSchedule = []
        invPileOriginList =[]
        i = 0
                        
        while days_rem > 0 and bud_rem > 0:
            for ii in range(List_tests_per_day[i]):  #sample for amount of tests/ day
                madeTest= False
                rand  = random.random()  
                if rand < 1-Epsilon: #exploit
                    currNode = List_EndNode[List_False_Rates.index(max(List_False_Rates))]
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
                        if randUnif > diagnosticSensitivity:
                            DPRoot = 0 # Generate a false negative

                    
                    sampleSchedule.append([i, currNode.id])
                    List_Nodes_Tested.append(currNode.id)
                    if DPRoot > -1:
                        List_Report[currNode.id-12][0] += DPRoot
                    List_Report[currNode.id-12][1]= List_Nodes_Tested.count(currNode.id)
                    List_Report[currNode.id-12][2] = List_Report[currNode.id-12][0]/List_Report[currNode.id-12][1]
                    List_False_Rates[currNode.id-12] =List_Report[currNode.id-12][2]
                else: #explore
                    currNode =random.choice(List_EndNode)
                    sampleSchedule.append([i, currNode])
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
                        if randUnif > diagnosticSensitivity:
                            DPRoot = 0 # Generate a false negative

                    
                    sampleSchedule.append([i, currNode.id])
                    List_Nodes_Tested.append(currNode.id)
                    if DPRoot > -1:
                        List_Report[currNode.id-12][0] += DPRoot
                    List_Report[currNode.id-12][1]= List_Nodes_Tested.count(currNode.id)
                    List_Report[currNode.id-12][2] = List_Report[currNode.id-12][0]/List_Report[currNode.id-12][1]
                    List_False_Rates[currNode.id-12] =List_Report[currNode.id-12][2]            

                #List_Nodes_Tested.append([currNode.id, DPRoot])
            if List_tests_per_day[i] == floor_bud_over_days:
                days_rem-=1
                bud_rem -= floor_bud_over_days
                i += 1
            elif List_tests_per_day[i] == floor_bud_over_days+1:
                days_rem-=1
                bud_rem -= floor_bud_over_days+1
                i+=1
print(List_Nodes_Tested, List_False_Rates, List_Report)
    
    
    ### TESTING POLICIES GENERATED HERE ###
    # Generate sampling schedule for current graph, as well as a report table shell