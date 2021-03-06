"""
Created on Tue Nov 19 21:50:21 2019
@author: Eugene Wickett

Stores classes for use with 'Falsification Simulator.py'
"""
import numpy as np
import random

class DrugPacket:
    _numDP = 0 # Counter for our DPs
    
    def __init__(self, size=0, nodeIDList=[0], intNode=0):
        DrugPacket._numDP += 1
        self.id = DrugPacket._numDP #Integer
        self.size = size #Integer
        self.nodeIDList = nodeIDList #List
        if len(self.nodeIDList)>1:     
            self.intNode = self.nodeIDList[1]   

        
class Node:
    def __init__(self, nodeNum=0, reorderPoint=0, reorderAmount=0, preferredSupplierVec=[0], preferredSupplierLTsVec=[0], preferredSupplierRsVec=[0], endNodeTF = False, falseProb = 0, demandSched=[], prev_id_int=0): #MZ added prev_id_lst
        self.id = nodeNum
        self.r = reorderPoint
        self.R = reorderAmount
        self.PreferenceList = preferredSupplierVec # A list of node numbers, in order of preference
        self.PreferenceLTsList = preferredSupplierLTsVec # A list of lead times
        self.PreferenceRsList = preferredSupplierRsVec
        self.FalsifierProbability = falseProb
        self.InventoryLevel = [] # List of levels corresponding to each DPID 
        self.InventoryDPID = [] # List of DP IDs
        self.InventoryPeriodsRemaining = [] # List of time periods remaining before corresponding inventory arrives; "0" means inventory is available for further distrubtion
        self.demandResults = [0,0] #List of demand results [Satisfied,Stockout] in number of units
        self.stockoutDays = 0 # Number of days stocked out
        self.EndNode = endNodeTF
        self.demandSched = demandSched # List of demands for each day; generated by demandScheduleGenerator(); only if EndNode = True
        self.DaysStockedOut = 0 # Tracker of how many simulation days the node has been stocked out
        self.prev_id_int = prev_id_int #MZ list of the id's the DP had at intermediate nodes before reaching 
    
    def ProcessIncomingOrder(self): # Update incoming order time periods remaining by 1
        if np.sum(self.InventoryPeriodsRemaining) > 0: # We have an order in transit
            # Iterate through the list and drop each days remaining by one
            for invPile in range(len(self.InventoryPeriodsRemaining)):
                if self.InventoryPeriodsRemaining[invPile] > 0:
                    self.InventoryPeriodsRemaining[invPile] = self.InventoryPeriodsRemaining[invPile] - 1
        #MZ updates inventory periods remaining. Decreases amount of days remaining (periods) by one, if there are more than zero days remaining before inventory can arrive. 
    
    def ProcessScheduleDemand(self, simDay=0, DPList=[], RootList=[], rootConsumptionVec=[]): # Process demand for an end node with a demand schedule
        todayDemand = int(self.demandSched[simDay]) # Demand for today
        todayInventoryAvailable = 0
        for invPile in range(len(self.InventoryDPID)): # Iterate through our inventory piles and count how much is presently available
            if self.InventoryPeriodsRemaining[invPile] == 0:
                todayInventoryAvailable = todayInventoryAvailable + self.InventoryLevel[invPile] #*MZ adds up how much inventory is available and can be distributed today. 
        # Take the minimum of what's available or the demand; balance demand becomes stockout
        inventoryToTake = min(todayDemand,todayInventoryAvailable)
        stockoutDemand = todayDemand - inventoryToTake #MZ how much more inventory do we need to satisfy all of today's demand.
        self.demandResults[1] += stockoutDemand # Increment dissatisfied customers
        currInvPile = 0
        while inventoryToTake > 0:
            if (self.InventoryPeriodsRemaining[currInvPile]==0) and (self.InventoryLevel[currInvPile]>0): # Inventory is available and can be distributed today
                pileAmountToTake = min(inventoryToTake,self.InventoryLevel[currInvPile]) # Can only take inventory if it's there #* MZ is it even possible fot inventoryToTake to be greater than inventorylevel?
                self.InventoryLevel[currInvPile] -= pileAmountToTake # Increment down the inventory pile
                inventoryToTake -= pileAmountToTake # Increment amount left to take
                self.demandResults[0] += pileAmountToTake # Increment satisfied customers
                currDPID = self.InventoryDPID[currInvPile] # Grab the DPID for obtaining statistics
                for indDP in DPList:
                    if indDP.id == currDPID: # Found the matching DPID
                        rootSupplier = indDP.nodeIDList[0] #root node will be first value of node id list
                        rootSupplierInd = RootList.index(rootSupplier) #finding which index in rootList corresponfs with this rootSupplier id
                        rootConsumptionVec[rootSupplierInd] += pileAmountToTake # Add pile amount to root nodes stat. MZ: The amount consumed.
                    
            currInvPile += 1 # Go to the next inventory pile
        ### END WHILE LOOP FOR INVENTORY
        return DPList, RootList, rootConsumptionVec  #*MZ: seems like only value changing from fxn would be rootConsumptionVec, why return all?
        #removes available inventory from stock pile for each individual drug with inventory available to distribute to demand. 
        #Also updated the amount of rootconsumptionvector and adds the amount of drug packets distributed to the root of the drug we are taking from.
        
    def MakeOrder(self, rootList, intermediateList, DPList): # Check if current inventory levels are leq the reorder point, and make orders according to the preferred supplier list
        if np.sum(self.InventoryLevel) <= self.r: # We place an order
            # Iterate through the prefered suppliers list until we hit a root node
            Rremaining = self.R # Initialize the desired reorder amount
            suppInd = 0 # Supplier number is distinct from where the supplier's LT is in the list
            for currSupplier in self.PreferenceList:
                currLT = self.PreferenceLTsList[suppInd] # LT from this supplier
                suppInd += 1
                # Check if it is a root node
                isRoot = False
                for root in rootList: #MZ: loop thru all roots in rootList, see if any root values match the currsupplier
                    if root==currSupplier: # We have a root
                        isRoot = True
                # Place order if supplier is a root; otherwise need to check if inventory is available
                if isRoot == True:
                    # If this node has a positive falsifier procurement probability,
                    # generate a random uniform and set the current supplier to 
                    # the falsifier node (node 1) if the uniform is below the threshold
                    if self.FalsifierProbability > random.uniform(0,1):  #MZ: use UnifRandom var to randomly see if the currsupplier is as falsifier
                        currSupplier = 1 #MZ: denotes currSupplier is a falsifier
                    
                    
                    newDP = DrugPacket(size=Rremaining, nodeIDList=[currSupplier,self.id]) 
                    self.InventoryLevel.append(newDP.size)
                    self.InventoryDPID.append(newDP.id)
                    self.InventoryPeriodsRemaining.append(currLT)
                    DPList.append(newDP)
                    #MZ: if the preferred supplier is at a root node, add new drug packet(from global supplier) w the reorder amount of drug packets to the inventory list. 
                    
                else: # Need to order as much as possible from this intermediate node
                    for intermediate in intermediateList:
                        if intermediate.id == currSupplier: # We have located our intermediate node
                            supplierInventoryAvailable = 0
                            for invPile in range(len(intermediate.InventoryDPID)): # Iterate through our inventory piles and count how much is presently available
                                if intermediate.InventoryPeriodsRemaining[invPile] == 0:
                                    supplierInventoryAvailable += intermediate.InventoryLevel[invPile] # Increment the available inventory from this intermediate node
                            # Take the minimum of what's available or the demand; balance demand becomes stockout
                            inventoryToTake = min(Rremaining,supplierInventoryAvailable)
                            Rremaining -= inventoryToTake
                            stockoutDemand = max(Rremaining,0) #******MZ: i feel  like  Rremaining is getting subtracted twice. Here the updated value of Rremaining is being used. # Record stockout demand if taken inventory is less than the remaining reorder amount.  Yes I was right!
                            intermediate.demandResults[1] += stockoutDemand # Increment dissatisfied customers
                            currInvPile = 0 # Start with intermediary's first inventory pile
                            while inventoryToTake > 0:
                                if (intermediate.InventoryPeriodsRemaining[currInvPile]==0) and (intermediate.InventoryLevel[currInvPile]>0): # Inventory is available
                                    pileAmountToTake = min(inventoryToTake,intermediate.InventoryLevel[currInvPile]) # Can only take inventory if it's there
                                    tempDPID = intermediate.InventoryDPID[currInvPile]
                                    # Need to find the current DP
                                    for DPobj in DPList:
                                        if DPobj.id==tempDPID: #Found it. store as new class characteristic (parent). initialize as empty list. what was dpid right b4 got here (at intermediate node). histogram. see how reorder policies at intermediate nodes effect histogram.
                                            oldDPList = DPobj.nodeIDList
                                            self.prev_id_int = DPobj.id #MZ added. Store prev_id_list of current DP 
                                    tempDPIDList = oldDPList + [self.id]
                                    newDP = DrugPacket(size=pileAmountToTake,nodeIDList=tempDPIDList)
                                    self.InventoryLevel.append(newDP.size) # Add pile to inventory list
                                    self.InventoryDPID.append(newDP.id)
                                    self.InventoryPeriodsRemaining.append(currLT)
                                    DPList.append(newDP)
                                    intermediate.InventoryLevel[currInvPile] -= pileAmountToTake # Increment down the inventory pile
                                    inventoryToTake -= pileAmountToTake # Increment amount left to take
                                    intermediate.demandResults[0] += pileAmountToTake # Increment satisfied customers
                                # END IF FOR AVAILABLE INVENTORY IN PILE
                                currInvPile += 1 # Go to the next inventory pile #*** MZ: if you never take all of inventoryToTake, the currInvPile index will be out o bounds for inventory level. 
                            ### END WHILE LOOP FOR INVENTORY
                                                   
        return intermediateList, DPList # We may have made changes to these lists
        #MZ: Finding preferred node to take inventory from. If its a root node, you can take as much as you demand. Add this new inventory as a new DP with the size you detemine
        #Otherwise, if preferred node is an inventory node, decrease inventory levels by the amount the preferred supplier has available. Decrease this amount from  Rremaining. Add the amount you can take from this intermediate node as a new drug packet. 
                            
                            
                            
                            
                            
                            
                            
        
        
        
        
        
        
        
        
        
        
        
        