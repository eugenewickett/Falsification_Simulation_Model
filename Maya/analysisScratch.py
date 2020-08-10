# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""

import numpy as np
import matplotlib.pyplot as plt
import os # for directories
import pickle # for saving/loading objects in Python
import Falsification_Sim_Modules as simModules

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
dirStr = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\Maya'
os.chdir(dirStr) # Set directory    

fileName = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\Maya\\Fals_Sim_debug_OUTPUT'
OPdict = pickle.load(open(fileName, 'rb'))

simModules.SimReplicationOutput(OPdict)
OPdict[4]






