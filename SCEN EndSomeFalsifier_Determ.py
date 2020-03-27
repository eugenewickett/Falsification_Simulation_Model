# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 14:59:35 2019

@author: Eugene Wickett

This code is for building and running simulations of a supply chain susceptible
to falsification and substandardization.
"""

import numpy as np
import matplotlib.pyplot as plt
import random
import csv # for csv manipulation
import time # for time tracking
import os # for directories
from tabulate import tabulate # for making outputs


# Run supporting files
#os.getcwd() Run this command to get the current working directory string
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

import Falsification_Sim_Classes as simClasses # our class objects and methods
import Falsification_Sim_Modules as simModules # modules for the simulation

# Read input lists: node list, arc preferences, and arc lead times
nodeInputFileString = 'LIB_Nodes_1.csv'
arcPreferencesFileString = 'LIB_Arcs_Preferences_1.csv'
arcLTsFileString = 'LIB_Arcs_LTs_1.csv'
arcRsFileString = 'LIB_Arcs_Rs_1.csv'
# Enter the length of the simulation
NumSimDays = 600
samplingBudget = NumSimDays*3

scenarioStr = 'EndSomeFalsifier-Determ' # For file saving

# Generate the lists of root, intermediate, and end nodes; also keep the list of node headers
List_RootNode, List_IntermediateNode, List_EndNode, nodeListHeader, nodeList, nodeNum, arcPreferencesMatrix, arcLTsMatrix, arcRsMatrix = simModules.generateNodeListsFromFile(nodeInputFileString,arcPreferencesFileString,arcLTsFileString,arcRsFileString, NumSimDays)
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
List_TestingSchedule = simModules.testingScheduleGenerator(nodes=nodeList, int_numDays=NumSimDays, int_sampleBudget = samplingBudget, int_PolicyType = 1, arr_PolicyParameter = [0])
TestReportTbl = []


### THINGS TO CHANGE AND PLAY WITH ONCE THE LISTS ARE GENERATED ###

# Intermediate nodes - put something inside this loop
for indInt in range(intermediateNum):
    currIntermediate = List_IntermediateNode[indInt]
    currIntermediate.r = 1500 # Change from 800
    currIntermediate.R = 7500 # Change from 5000
# End nodes - put something inside this loop
for indEnd in range(endNum):
    currEnd = List_EndNode[indEnd]
    
    
# Remove tendency to procure from falsifiers for random 50% of end nodes
for indEnd in range(endNum):
    currEnd = List_EndNode[indEnd]
    if np.random.uniform() > 0.5:
        tempHolder = currEnd.PreferenceList.pop()
        tempHolder = currEnd.PreferenceLTsList.pop()
        tempHolder = currEnd.PreferenceRsList.pop()
        
    




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
            for invPile in range(len(currNode.InventoryDPID)):
                if (currNode.InventoryPeriodsRemaining[invPile] == 0) and (currNode.InventoryLevel[invPile] > 0): #Test
                    madeTest = True
                    testDPID = currNode.InventoryDPID[invPile] #DPID for this inventory pile
                    for indDP in List_DP:
                        if testDPID == indDP.id:
                            DPRoot = indDP.nodeIDList[0] # Source of this DP
                    currNode.InventoryLevel[invPile] -= 1 # "Buy" the sample
                    currNode.demandResults[0] += 1 # "Satisifed" demand
                    break
            if madeTest == False: # Report -1 if no sample procured due to stockout
                DPRoot = -1
            # Save the result to the reporting table
            newTestRow = [today,currTestNodeID,DPRoot]
            TestReportTbl.append(newTestRow) 
        
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
# END OF SIMULATIONS DAYS LOOP

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
for sample in TestReportTbl:
    if not sample[1] in testNodeList:
        testNodeList.append(sample[1])
testNodeList.sort()
TestSummaryTbl = []
for testedNode in testNodeList:
    TestSummaryTbl.append([testedNode,0,0,0]) # Times tested, # falsified, # stocked out
for sample in TestReportTbl:
    indNodeList = testNodeList.index(sample[1]) # The row of the summary table
    TestSummaryTbl[indNodeList][1] += 1
    if sample[2] == 1: # From falsified root
        TestSummaryTbl[indNodeList][2] += 1
    if sample[2] == -1: # Stocked out
        TestSummaryTbl[indNodeList][3] += 1
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

# Simulation end time
totalRunTime = [time.time() - startTime]

# END OF RX REPORTING LOOP    

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
print(tabulate(TestSummaryTbl, headers=['Node ID', 'Number of Samples', 'Number Falsified', 'Number Stocked Out']))
print('*'*100)
print('OVERALL SUMMARY STATISTICS')
print('*'*100)
print(tabulate([totalRunTime+TestOverallSummary],headers=['Simulation Run Time','Total Tests','Falsifieds Found','Stocked Out Found']))

# GENERATE PLOTS
# Root consumption %'s
fig = plt.figure()
ax = fig.add_axes([0.2,0.2,0.3,0.5])
ax.set_xlabel('Root Node',fontsize=16)
ax.set_ylabel('Percentage consumption',fontsize=16)
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.bar(Root_Plot_x,Root_Plot_y)
plt.savefig('RootConsumption_' + scenarioStr + '.png')
plt.show()


# Intermediate stockout %'s
fig = plt.figure()
ax = fig.add_axes([0.1,0.4,0.8,0.5])
ax.set_xlabel('Intermediate Node',fontsize=16)
ax.set_ylabel('Percentage stocked out',fontsize=16)
#vals = ax.get_yticks()
#ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.bar(Intermediate_Plot_x,Intermediate_Plot_y)
plt.savefig('IntermediateStockout_' + scenarioStr + '.png')
plt.show()

# End stockout %'s
fig = plt.figure()
ax = fig.add_axes([0.3,0.1,3,0.5])
ax.set_xlabel('End Node',fontsize=16)
ax.set_ylabel('Percentage stocked out',fontsize=16)
#vals = ax.get_yticks()
#ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.bar(End_Plot_x,End_Plot_y)
plt.savefig('EndStockout_' + scenarioStr + '.png')
plt.show()

# Testing results
Test_Plot_x = np.array(Test_Plot_x)
fig = plt.figure()
ax = fig.add_axes([0.2,0.2,3,0.5])
ax.set_xlabel('Tested Node',fontsize=16)
ax.set_ylabel('Result Amount',fontsize=16)
ax.bar(Test_Plot_x+0.00, Test_Plot_y1, color='black', width=0.25)
ax.bar(Test_Plot_x+0.25, Test_Plot_y2, color='red', width=0.25)
ax.bar(Test_Plot_x+0.50, Test_Plot_y3, color='gray', width=0.25)
plt.legend(('Times tested','Times falsified','Times stocked out'))
plt.savefig('TestingResults_' + scenarioStr + '.png')
plt.show()



# Code for line plots
#fig, ax = plt.subplots()
#plt.plot(Root_Plot_x, Root_Plot_x, 'red')
#fig.suptitle('Invested Capacity vs. Algorithm Iteration', fontsize=14)
#ax.set_xlabel('Iteration', fontsize=12)
#ax.set_ylabel('Invested Capacity', fontsize=12)
#plt.legend(("K1","K2","K3"))
#plt.savefig('capacity_vs_iterations.png')





# WRITE TO CSV FILE
with open('SCEN ' + scenarioStr + '.csv', 'w') as outputFile:
    outputFile.write(tabulate([totalRunTime+TestOverallSummary],headers=['Simulation Run Time','Total Tests','Falsifieds Found','Stocked Out Found']))
    outputFile.write('\n')
    outputFile.write(tabulate(TestSummaryTbl, headers=['Node ID', 'Number of Samples', 'Number Falsified', 'Number Stocked Out']))
    outputFile.write('\n')
    outputFile.write(tabulate(RootReportTbl, headers=['Root ID', 'End Consumption', 'Overall Percentage']))
    outputFile.write('\n')
    outputFile.write(tabulate(IntermediateReportTbl, headers=['Intermediate Node ID', 'Satisfied Demand', 'Stockout Demand','Stockout %','Inventory Available']))
    outputFile.write('\n')
    outputFile.write(tabulate(EndReportTbl, headers=['End Node ID', 'Satisfied Demand', 'Stockout Demand', 'Stockout Days','Stockout %', 'Inventory Available']))
    
#with open('SCEN 1 Results.csv', 'w', newline='') as file:
    #writer = csv.writer(file)
    #writer.writerow(["RxID", "GoodConsumed", "SubstandardConsumed", "FalsifiedConsumed", "StockoutConsumed"]) # Initialize the header
    #for rw in range(endNum):
    #writer.writerow(tabulate(RootReportTbl, headers=['Root ID', 'End Consumption', 'Overall Percentage']))


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




