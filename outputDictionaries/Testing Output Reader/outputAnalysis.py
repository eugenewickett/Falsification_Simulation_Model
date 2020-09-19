# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""
import os # for directories
import pickle # for saving/loading objects in Python

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\outputDictionaries\\Testing Output Reader 2'
dirStr2 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model'
os.chdir(dirStr2) # Set directory   
import Falsification_Sim_Modules as simModules
os.chdir(dirStr) # Set directory   

strVec = ['OP_Static','OP_NUTSthresh_100_0.3','OP_NUTSthresh_200_0.3']

outputDicts = []
for strEnd in strVec:
    fileName = dirStr + '\\' + strEnd  
    currDict = pickle.load(open(fileName, 'rb'))
    outputDicts.append(currDict)

simModules.SimSFEstimateOutput(outputDicts,strVec,0.30)










