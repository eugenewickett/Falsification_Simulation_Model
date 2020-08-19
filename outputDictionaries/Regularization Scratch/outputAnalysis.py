# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""

import os # for directories
import pickle # for saving/loading objects in Python
import Falsification_Sim_Modules as simModules

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\outputDictionaries\\Regularization Scratch'
os.chdir(dirStr) # Set directory    

strVec = ['001','01','05','d1','d5','1','10','100','1000']

outputDicts = []
for strEnd in strVec:
    fileName = dirStr + '\\OP_' + strEnd  
    currDict = pickle.load(open(fileName, 'rb'))
    outputDicts.append(currDict)

simModules.SimSFEstimateOutput(outputDicts)











