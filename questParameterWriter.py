# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:42:01 2020

@author: eugen
"""
import csv

'''
sens = [0.99,0.95,0.9,0.85,0.8]
spec = [0.99,0.95,0.9,0.85]
sampleSize = [1,5,10]
simDays = [600]
globalDem = [0.,40.,80.]
sampPols = [0.]
regulWts = [0.5]
'''

sens = [0.85,0.95,0.9]
spec = [0.99,0.95,0.9]
sampleSize = [0.25,0.5,1]
simDays = [600]
globalDem = [0.,40.,80.,120]
sampPols = [0.]
regulWts = [0.1]



with open('parameterFile.csv', 'w') as f:
    for se in sens:
        for sp in spec:
            for samps in sampleSize:
                for days in simDays:
                    for dem in globalDem:
                        for pol in sampPols:
                            for wt in regulWts:
                                l = []
                                l.extend([se,sp,samps,days,dem,pol,wt])
                                l = [str(l[i]) for i in range(len(l))]
                                wrtr = csv.writer(f,lineterminator = '\n')
                                wrtr.writerow(l)
                                #f.write('\t'.join(l) + '\n')
                        

