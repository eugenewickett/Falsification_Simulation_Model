# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 00:42:07 2020

@author: eugen
"""

import os # for directories
import pickle # for saving/loading objects in Python
import matplotlib.pyplot as plt
import numpy as np
# Run supporting files
#os.getcwd() Run this command to get the current working directory string
os.chdir('C:\\Users\\eugen\\OneDrive\\Documents\\EAGER Project\\Simulator\\Falsification_Simulation_Model') # Set directory    

import SFP_Sim_Helpers as simHelpers

directory = r'C:\Users\eugen\OneDrive\Documents\EAGER Project\Simulator'+\
         '\Sim Model Files TOO BIG FOR REPO\THIRTEENDEC_RUNS\Consolidated Files'
         
OPFileNames = os.listdir(directory)
'''
GROUP 1
OPFileNames=['OP_600_0.25_0.7_0.95_0.0_0.0_0.2_0.0',
             'OP_600_0.25_0.7_0.95_0.0_0.0_0.2_1.0',
             'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_0.0',
             'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_1.0',
             'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_0.0',
             'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_1.0',
             'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_0.0',
             'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_1.0',
              'OP_600_1.0_0.9_0.99_0.0_0.0_0.2_0.0',
              'OP_600_1.0_0.9_0.99_0.0_0.0_0.2_1.0',
             ]
GROUP 2
OPFileNames=['OP_600_0.25_0.7_0.95_0.0_0.0_0.2_0.0',
             'OP_600_0.25_0.7_0.95_0.0_0.0_0.2_1.0',
             'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_0.0',
             'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_1.0',
             'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_0.0',
             'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_1.0',
             'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_0.0',
             'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_1.0',
             'OP_600_1.0_0.9_0.99_0.0_0.0_0.2_1.0',
             'OP_600_1.0_0.9_0.99_80.0_0.0_0.2_1.0'
             ]

 
 'OP_600_0.25_0.7_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.7_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.7_0.99_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.7_0.99_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.7_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.7_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.7_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.99_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.99_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.8_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.8_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.99_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.99_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.99_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.99_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.99_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.99_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.25_0.9_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.25_0.9_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.99_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.99_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.7_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.7_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.95_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.8_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.95_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.8_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.99_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.8_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.9_0.0_0.0_0.2_0.0',
 'OP_600_0.5_0.8_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.8_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.8_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.99_0.9_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.95_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.95_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.99_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.99_80.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.9_0.0_0.0_0.2_1.0',
 'OP_600_0.5_0.9_0.9_80.0_0.0_0.2_0.0',
 'OP_600_0.5_0.9_0.9_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.95_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.95_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.99_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.99_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.9_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.7_0.9_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.95_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.95_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.99_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.99_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.9_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.8_0.9_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.95_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.95_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.99_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.99_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.9_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.99_0.9_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.95_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.95_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.99_0.0_0.0_0.2_0.0',
 'OP_600_1.0_0.9_0.99_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.99_80.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.9_0.0_0.0_0.2_1.0',
 'OP_600_1.0_0.9_0.9_80.0_0.0_0.2_1.0']
'''

OPDicts = []
for item in OPFileNames:
    if len(OPDicts) < 10:
        fileName1 = directory + '\\'+item
        outputDict1 = pickle.load(open(fileName1, 'rb'))
        OPDicts.append(outputDict1)


simHelpers.SimSFPEstimateOutput(OPDicts,OPFileNames)


numSampsColleced = []
for i in range(len(OPDicts)):
    currDict = OPDicts[i]
    
    for j in range(len(currDict.keys())):
        numSampsColleced.append(np.sum([item[1] for item in currDict[j]['dynTestResults']]))
plt.hist([i for i in numSampsColleced if i > 501])
np.mean(numSampsColleced)
'''
timeVec = []
for i in range(len(OPDicts)):
    currDict = OPDicts[i]
    
    for j in range(len(currDict.keys())):
        timeVec.append(currDict[j]['simRunTime'])

plt.hist(timeVec)
np.mean(timeVec)/60


OPDicts[2][0]['inputParameterDictionary']
OPDicts[2][0]['intSFTrueValues']
np.sum([i[1] for i in OPDicts[0][5]['dynTestResults']])

for imp in range(10):
    print(np.mean(simModules.invlogit(OPDicts[2][0]['postSamps_Untracked'][:,imp])))
    print(np.mean(simModules.invlogit(OPDicts[2][0]['postSamps_Tracked'][:,imp])))
'''


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