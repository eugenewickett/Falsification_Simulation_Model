# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 19:15:37 2020

@author: eugen
"""
import numpy as np
import matplotlib.pyplot as plt
import random #for seeds
import sys
import csv # for csv manipulation
import time # for time tracking
import os # for directories
from tabulate import tabulate # for making outputs
import pickle # for saving/loading objects in Python
import winsound # for beeps

# Run supporting files
#os.getcwd() Run this command to get the current working directory string
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

import Falsification_Sim_Classes as simClasses # our class objects and methods
import Falsification_Sim_Modules as simModules

# Standard model
fileName1 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Falsification_Simulator_OUTPUT'
outputDict1 = pickle.load(open(fileName1, 'rb'))
simModules.SimReplicationOutput(outputDict1)

# One importer with 20% falsified batches
fileName2 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\SCEN1_OUTPUT'
outputDict2 = pickle.load(open(fileName2, 'rb'))
simModules.SimReplicationOutput(outputDict2)

# Heightened impatience of outlets
fileName3 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\SCEN2_OUTPUT'
outputDict3 = pickle.load(open(fileName3, 'rb'))
simModules.SimReplicationOutput(outputDict3)

# Heightened impatience of outlets
fileName4 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\SCEN3_OUTPUT'
outputDict4 = pickle.load(open(fileName4, 'rb'))
simModules.SimReplicationOutput(outputDict4)