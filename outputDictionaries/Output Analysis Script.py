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

# BAD INTERMEDIATE NODE
OPFileNames = ['OP_Static_200_0.6_0.99_0.99_0.0_33946773.29296132',\
               'OP_Static_200_0.9_0.99_0.99_0.0_6143606.200114578',\
               'OP_Static_200_1.2_0.99_0.99_0.0_6143606.200414138']
OPDicts = []
for item in OPFileNames:
    fileName1 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator'+\
    '\\Sim Model Files TOO BIG FOR REPO\\Sampling Budgets and Diagnostics\\'+item
    outputDict1 = pickle.load(open(fileName1, 'rb'))
    OPDicts.append(outputDict1)

simModules.SimSFEstimateOutput(OPDicts)

'''
for i in range(5):
    print(OPDicts[2][i]['falsePerc_LklhdSamples'])
'''