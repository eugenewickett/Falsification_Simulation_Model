# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

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



# BAD INTERMEDIATE NODE
fileName1 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Falsification_Simulator_BADINT_OUTPUT'
outputDict1 = pickle.load(open(fileName1, 'rb'))

# BAD END NODES
fileName2 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Falsification_Simulator_BADENDS_OUTPUT'
outputDict2 = pickle.load(open(fileName2, 'rb'))

# BAD END NODES
fileName3 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Falsification_Simulator_BADENDS2_OUTPUT'
outputDict3 = pickle.load(open(fileName3, 'rb'))



simModules.SimReplicationOutput(outputDict1)
simModules.SimReplicationOutput(outputDict2)
simModules.SimReplicationOutput(outputDict3)

