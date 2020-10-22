# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""

import os # for directories
import pickle # for saving/loading objects in Python

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

import Falsification_Sim_Modules as simModules

directory = r'C:\Users\eugen\OneDrive\Documents\EAGER Project\Simulator'+\
         '\Sim Model Files TOO BIG FOR REPO\Consolidated Budgets & Diagnostics'
# BAD INTERMEDIATE NODE
OPFileNames = os.listdir(directory)
OPDicts = []
for item in OPFileNames:
    fileName1 = directory + '\\'+item
    outputDict1 = pickle.load(open(fileName1, 'rb'))
    OPDicts.append(outputDict1)


OPDicts1 = OPDicts[0:25]
OPDicts2 = OPDicts[25:50]
OPDicts3 = OPDicts[50:75]
OPDicts4 = OPDicts[75:]

simModules.SimSFEstimateOutput(OPDicts3)

'''
for i in range(5):
    print(OPDicts[2][i]['intSFTrueValues'])

len(OPDicts)
len(OPDicts[6])
OPDicts[8][3]['intFalseEstimates']
'''