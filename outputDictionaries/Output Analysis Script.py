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
         '\Sim Model Files TOO BIG FOR REPO\SIXNOV_RUNS\Consolidated Files'
# BAD INTERMEDIATE NODE
OPFileNames = os.listdir(directory)
'''
OPFileNames=['OP_600_10.0_0.99_0.85_0.0_0.0_0.5',
             'OP_600_10.0_0.99_0.95_0.0_0.0_0.5',
             'OP_600_10.0_0.99_0.9_0.0_0.0_0.5',
             'OP_600_10.0_0.99_0.99_0.0_0.0_0.5',
             'OP_600_10.0_0.99_0.85_80.0_0.0_0.5',
             'OP_600_10.0_0.99_0.95_80.0_0.0_0.5',
             'OP_600_10.0_0.99_0.9_80.0_0.0_0.5',
             'OP_600_10.0_0.99_0.99_80.0_0.0_0.5'
             ]
'''

OPDicts = []
for item in OPFileNames:
    if len(OPDicts) < 10:
        fileName1 = directory + '\\'+item
        outputDict1 = pickle.load(open(fileName1, 'rb'))
        OPDicts.append(outputDict1)


simModules.SimSFEstimateOutput(OPDicts,OPFileNames)


OPDicts[2][0]['inputParameterDictionary']



'''
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