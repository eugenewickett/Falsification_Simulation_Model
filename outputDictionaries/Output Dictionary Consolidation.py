# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""
'''
This script combines output dictionaries with identical inputs and stores them
as new dictionaries. You need to set the directory location
'''
import os # for directories
import pickle # for saving/loading objects in Python



#os.getcwd() Run this command to get the current working directory string
# SET DIRECTORY LOCATION HERE
directory = r'C:\Users\eugen\OneDrive\Documents\EAGER Project\Simulator\Sim Model Files TOO BIG FOR REPO\THIRTEENDEC_RUNS\Raw Files\temp' # Set directory    
os.chdir(directory)

# ONLY COMBINE DICTIONARIES IF THEY HAVE THE SAME LEADING STRING
# First grab all unique string beginnings         
nameList = []
for filename in os.listdir(directory):
    lastUnderInd = filename.rfind('_')
    currBegin = filename[:lastUnderInd]
    if not currBegin in nameList:
        # Add it to the unique names list
        nameList.append(currBegin)


# Now combine all files starting with the same name in nameList
for nameStr in nameList:
    # Generate new dictionary
    newOPDict = {}
    for filename in os.listdir(directory):
        if nameStr == filename[:len(nameStr)]:
            # Add replications to new dictionary
            dirLoc = directory + '\\' + str(filename)
            currOPdict = pickle.load(open(dirLoc, 'rb'))
            for repNum in currOPdict.keys():
                newOPDict[len(newOPDict)] = currOPdict[repNum]
    # Now we store the new dictionary
    pickle.dump(newOPDict,open(nameStr,'wb'))
        
len(nameList)
