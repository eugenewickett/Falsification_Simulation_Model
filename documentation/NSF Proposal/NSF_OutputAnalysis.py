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
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\NSF Proposal') # Set directory     
 
import Falsification_Sim_Classes as simClasses # our class objects and methods 
import Falsification_Sim_Modules as simModules 
 
 
 
# Output A 
fileName1 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\NSF Proposal\\output dictionaries\\NSF_OUTPUTA' 
outputDict1 = pickle.load(open(fileName1, 'rb')) 
 
# Output B 
fileName2 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\NSF Proposal\\output dictionaries\\NSF_OUTPUTB' 
outputDict2 = pickle.load(open(fileName2, 'rb')) 
 
# Output C 
fileName3 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\NSF Proposal\\output dictionaries\\NSF_OUTPUTC' 
outputDict3 = pickle.load(open(fileName3, 'rb')) 
 
# Output D 
fileName4 = 'C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model\\NSF Proposal\\output dictionaries\\NSF_OUTPUTD' 
outputDict4 = pickle.load(open(fileName4, 'rb')) 
 
 
simModules.SimReplicationOutput(outputDict1) 
simModules.SimReplicationOutput(outputDict2) 
simModules.SimReplicationOutput(outputDict3) 
simModules.SimReplicationOutput(outputDict4)