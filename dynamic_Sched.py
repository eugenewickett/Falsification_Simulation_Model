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
#NumSimDays = 1000
#samplingBudget = NumSimDays*2 
numReplications = 1 
#testPolicy = 3
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

    
    # Initialize the drug packet object list 
List_DP = [] 
List_Nodes_Tested = []

### TESTING POLICIES GENERATED HERE ### 
# Generate sampling schedule for current graph, as well as a report table shell 

#List_TestingSchedule = simModules.testingScheduleGenerator(nodes=nodeList, int_numDays=NumSimDays, int_sampleBudget = samplingBudget, int_PolicyType = testPolicy, arr_PolicyParameter = [0]) 
TestReportTbl = [] 
    
    ### THINGS TO CHANGE AND PLAY WITH ONCE THE LISTS ARE GENERATED ### 
    # Freeze end node demand to give intermediate nodes a chance to order and stock their shelves 
    # We do this by setting the first (maxLT) days at each node equal to 0 

    
    
    
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


def DynamicTestingScheduleGenerator(nodes = [], int_numDays=100,
                             int_sampleBudget = 100,
                             int_PolicyType = 3, 
                             arr_PolicyParameter = [0],   
                             List_False_Rates  = [1]*106,
                             epsilon = .11, 
                             show_plot_bool = True
                             ): 
    """
    MZ developed this function to create a dynamic testing schedule for sampling. 
    A random number between 0-1 is chosen and if the number is greater than 1-epsilon, a random node will be sampled. 
    However, if the random number is less than 1-epsilon, the node with the highest falsification rate will be sampled.
    This is based upon the epsilon greedy strategy. 
    """
        
    # Initialize the drug packet object list 
    List_DP = [] 
    List_Nodes_Tested = []
    startTime = time.time() 
    testid_lst  = []  #MZL 
    timestestedsum = 0 #MZL 
    list_prev_id = []  #MZ 
    parent_list =[] #MZ 
    unique_par_list  =[] #MZ 
    DPRoot_lst = []
    falsifieds_detected = 0 #counts the amount of falsifieds found
    ### TESTING POLICIES GENERATED HERE ### 
    # Generate sampling schedule for current graph, as well as a report table shell 

    #List_TestingSchedule = simModules.testingScheduleGenerator(nodes=nodeList, int_numDays=NumSimDays, int_sampleBudget = samplingBudget, int_PolicyType = testPolicy, arr_PolicyParameter = [0]) 
    TestReportTbl = []

        # Main simulation code; iterate over each day 
    for rep in range(numReplications): 

    # Generate the lists of root, intermediate, and end nodes; also keep the list of node headers 
        List_RootNode, List_IntermediateNode, List_EndNode, nodeListHeader, nodeList, nodeNum, arcPreferencesMatrix, arcLTsMatrix, arcRsMatrix = simModules.generateNodeListsFromFile(nodeInputFileString,arcPreferencesFileString,arcLTsFileString,arcRsFileString, int_numDays) 
        rootNum = len(List_RootNode) 
        intermediateNum = len(List_IntermediateNode) 
        endNum = len(List_EndNode) 
        
    # Simulation statistic lists 
        List_RootConsumption = [] 
        for rootInd in range(len(List_RootNode)): 
            List_RootConsumption.append(0) #MZ: list of zeros, size of list_rootnode     
        for today in range(int_numDays): 
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
            for indEnd in range(endNum): 
                currEnd = List_EndNode[indEnd] 
                List_IntermediateNode, List_DP = currEnd.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP) 
            
            # Intermediate nodes make orders if total inventory is leq the reorder point 
            for indInt in range(intermediateNum): 
                currIntermediate = List_IntermediateNode[indInt] 
                List_IntermediateNode, List_DP = currIntermediate.MakeOrder(rootList=List_RootNode,intermediateList=List_IntermediateNode,DPList=List_DP)               
            #print(List_DP)  


        #List_DP, List_RootNode, List_RootConsumption = currEnd.ProcessScheduleDemand(simDay=today,DPList=List_DP, RootList=List_RootNode, rootConsumptionVec=List_RootConsumption) 
        #print(List_DP)  
        
        #MZ MAIN CODE BELOW FOR EPSILON GREEDY STRATEGY
        if int_PolicyType  == 3: # MZ only run if policy type is 3
            #print('hi')
            dynamicbool =True #MZ
            bud_rem = int_sampleBudget #MZ initialize budget remaining to samplingBudget
            days_rem = int_numDays #MZ initialize sampling days remaining to NumSimDays
            bud_over_days = bud_rem/days_rem #MZ
            remainder_bud_over_days = int(bud_rem%days_rem) #MZ get the remainder. How many days in the simulation can afford 1 extra tests / day 
            floor_bud_over_days = int(bud_rem//days_rem) #MZ get the regular tests / day 
            List_tests_per_day = (([floor_bud_over_days]*(int_numDays-remainder_bud_over_days) + [floor_bud_over_days+1]*remainder_bud_over_days)) #MZ
            Epsilon = epsilon #MZ
            List_Report  = [[0,0,1] for x in range(106)] #MZ
            List_False_Rates = [1]*106  #MZ initialize list of falsification rates as every end node has 100% falsification rates until proven otherwise
            
            sampleSchedule = []
            invPileOriginList =[]
            
            i = 0
                            
            while days_rem > 0 and bud_rem > 0: #MZ continue sampling until no more $ or days available
                for ii in range(List_tests_per_day[i]):  #MZ sample for amount of tests/ day
                    madeTest= False
                    rand  = random.random()  #MZ initialize random number betweeen 0 and 1
                    if rand < 1-Epsilon: #MZ exploit
                        List_index_false_rates = [] #MZ
                        max_false_rate = max(List_False_Rates) #MZ exploit the node with the greatest falsification rate
                        if List_False_Rates.count(max_false_rate) > 1: #MZ if there are multiple nodes with the dame false rates, choose a random node among the ones with the highest false rate.                        
                            for j in range(len(List_False_Rates)): #MZ
                                if List_False_Rates[j] == max_false_rate: #MZ
                                    List_index_false_rates.append(j) #MZ
                            currNode = List_EndNode[random.choice(List_index_false_rates)] #MZ
                        else: #MZ
                            currNode = List_EndNode[List_False_Rates.index(max(List_False_Rates))] #MZ
                        currNodeid = currNode.id #MZ

                        #MZ Code taken from falsification simulator to test if a falsification is present at currNode 
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
                                            DPRoot = (indDP.nodeIDList)[0] # Source of this DP

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
                        #MZ end of code taken from falisification simulator
                        
                        sampleSchedule.append([i, currNode.id]) #MZ append day and id tested to sampleSchedule
                        List_Nodes_Tested.append(currNode.id) #MZ append node tested to keep track of frequency of each node
                        if DPRoot > -1: #MZ
                            List_Report[currNode.id-12][0] += DPRoot #MZ update how many falsifiers found at this end node. -12 to account for end node id's starting at 12, but list indices staring at 0.  
                            falsifieds_detected += DPRoot #MZ keep track of total falsifiers found
                        List_Report[currNode.id-12][1]= List_Nodes_Tested.count(currNode.id) #MZ update how many times this endnode has been tested
                        List_Report[currNode.id-12][2] = List_Report[currNode.id-12][0]/List_Report[currNode.id-12][1] #MZ update the falsificaiton rate of this end node
                        List_False_Rates[currNode.id-12] =List_Report[currNode.id-12][2] #MZ update falsification rate of this end node in list_false_rates so for the next iteration, the max falsification rate is accurate. 
                    
                    else: #MZ explore
                        currNode =random.choice(List_EndNode) #MZ choose a random endnode to explore and test
                        
                        #MZ Code taken from falsification simulator to test if a falsification is present at currNode: 
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
                                            DPRoot_lst.append(DPRoot)
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
                        #MZ end of code taken from falisification simulator
                        
                        sampleSchedule.append([i, currNode.id]) #MZ append the day and this end node id to sampleSchedule 
                        List_Nodes_Tested.append(currNode.id) #MZ
                        if DPRoot > -1: #MZ
                            List_Report[currNode.id-12][0] += DPRoot #MZ how many times have we seen this from a falsified root?
                            falsifieds_detected += DPRoot #MZ
                        List_Report[currNode.id-12][1]= List_Nodes_Tested.count(currNode.id) #MZ 
                        List_Report[currNode.id-12][2] = List_Report[currNode.id-12][0]/List_Report[currNode.id-12][1] #MZ 
                        List_False_Rates[currNode.id-12] =List_Report[currNode.id-12][2]    #MZ        

                    #List_Nodes_Tested.append([currNode.id, DPRoot])
                if List_tests_per_day[i] == floor_bud_over_days: #MZ
                    days_rem-=1 #MZ subtract a day
                    bud_rem -= floor_bud_over_days #MZ subtract the budget used this day
                    i += 1 #MZ update counter
                elif List_tests_per_day[i] == floor_bud_over_days+1: #MZ only necessary if there's a remainder between budget/samplingdays
                    days_rem-=1 #MZ subtract a day
                    bud_rem -= floor_bud_over_days+1 #MZ substract the budget used this day
                    i+=1 #MZ update counter
                
           
            #print(List_False_Rates, DPRoot)
            #print(List_Report, DPRoot_lst)
            
            #print(sampleSchedule)
            if show_plot_bool == True: #only display graphs if boolean is true
                List_EndNode_ID = [0]*106 #MZ
                for jj in range(len(List_EndNode)): #MZ
                    List_EndNode_ID[jj] =  List_EndNode[jj].id #MZ make a list of end node id's:
                plt.plot(List_EndNode_ID,List_False_Rates, '.') #MZ scatter plot showing end node vs falsification rate
                plt.xlabel("End Nodes") #MZ
                plt.ylabel("Falsification Rate") #MZ
                plt.show() #MZ
                plt.hist(List_Nodes_Tested, bins=200) #MZ shows frequency of each node tested
                plt.xlabel("End Nodes") #MZ 
                plt.ylabel("Amount of times tested") #MZ
                plt.show() #MZ
            print(falsifieds_detected) #MZ
        return sampleSchedule #MZ return sample schedule for falsifiacation simulator
        #return falsifieds_detected #MZ return number of falsifieds detected with this epsilon value
            #print(DPRoot_lst)
            #print(list_prev_id)
        
        ### TESTING POLICIES GENERATED HERE ###
        # Generate sampling schedule for current graph, as well as a report table shell

