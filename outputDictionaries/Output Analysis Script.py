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
         '\Sim Model Files TOO BIG FOR REPO'
# BAD INTERMEDIATE NODE
OPFileNames = os.listdir(directory)
OPDicts = []
for item in OPFileNames:
    fileName1 = directory + '\\'+item
    outputDict1 = pickle.load(open(fileName1, 'rb'))
    OPDicts.append(outputDict1)


OPDicts1 = OPDicts[0:20]
OPFileNames1 = OPFileNames[0:20]
OPDicts2 = OPDicts[20:40]
OPFileNames2 = OPFileNames[20:40]
OPDicts3 = OPDicts[40:60]
OPFileNames3 = OPFileNames[40:60]
OPDicts4 = OPDicts[60:80]
OPFileNames4 = OPFileNames[60:80]
OPDicts5 = OPDicts[80:]
OPFileNames5 = OPFileNames[80:]

simModules.SimSFEstimateOutput(OPDicts2,OPFileNames2)

OPDicts_PerfDiag_Names = ['OP_Static_200_0.6_0.99_0.99_0.0',\
                          'OP_Static_200_0.9_0.99_0.99_0.0',\
                          'OP_Static_200_1.2_0.99_0.99_0.0',\
                          'OP_Static_200_1.5_0.99_0.99_0.0',\
                          'OP_Static_200_1.8_0.99_0.99_0.0',\
                          'OP_Static_200_2.1_0.99_0.99_0.0',\
                          'OP_Static_200_2.4_0.99_0.99_0.0',\
                          'OP_Static_200_2.7_0.99_0.99_0.0',\
                          'OP_Static_200_3.0_0.99_0.99_0.0']
OPDicts_PerfDiag = []
for item in OPDicts_PerfDiag_Names:
    fileName1 = directory + '\\'+item
    outputDict1 = pickle.load(open(fileName1, 'rb'))
    OPDicts_PerfDiag.append(outputDict1)

simModules.SimSFEstimateOutput(OPDicts_PerfDiag,OPDicts_PerfDiag_Names)

'''
for i in range(5):
    print(OPDicts1[8][i+6]['intFalseEstimates'])

len(OPDicts)
len(OPDicts[6])
OPDicts[1][0]['intFalseEstimates_Plum']

currOutputLine = {'inputParameterDictionary':inputParameterDictionary,
                          'rootConsumption':List_RootConsumption,
                          'intDemandResults':List_demandResultsInt,
                          'endDemandResults':List_demandResultsEnd,
                          'testResults':TestReportTblToSend,
                          'dynTestResults':List_TestResults,
                          'intFalseEstimates':estIntFalsePercList,
                          'endFalseEstimates':estEndFalsePercList,
                          'intFalseEstimates_Bern':estIntFalsePercList_Bern,
                          'endFalseEstimates_Bern':estEndFalsePercList_Bern,
                          'intFalseEstimates_Plum':estIntFalsePercList_Plum,
                          'endFalseEstimates_Plum':estEndFalsePercList_Plum,
                          'intFalseEstimates_SampMLE':estIntFalsePercList_SampMLE,
                          'endFalseEstimates_SampMLE':estEndFalsePercList_SampMLE,
                          'falsePerc_LklhdSamples':estFalsePerc_LklhdSamples,
                          'intSFTrueValues':intSFVec,'endSFTrueValues':endSFVecCombo,
                          'simStartTime':startTime,
                          'simRunTime':totalRunTime
                          }
'''