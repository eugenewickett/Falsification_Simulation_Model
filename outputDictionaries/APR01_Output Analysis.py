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
fileName1 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\output dictionaries\\Fals_Sim_debug_OUTPUT_d'
outputDict1 = pickle.load(open(fileName1, 'rb'))

simModules.SimReplicationOutput(outputDict1)

