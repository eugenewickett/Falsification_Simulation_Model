# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:42:01 2020

@author: eugen
"""
import csv
'''
sens = [0.75,0.85,0.95,0.99]
spec = [0.75,0.85,0.95,0.99]
sampleMult = [1,5,10]
simDays = [100,400,700]
globalDem = [0.,35.,70.]
'''
sens = [0.99,0.95,0.9,0.85,0.8,0.75]
spec = [0.99,0.95,0.9,0.85]
sampleMult = [3,5,7]
simDays = [200,400,600]
globalDem = [0.,35.,70.]

with open('parameterFile.csv', 'w') as f:
    for se in sens:
        for sp in spec:
            for mult in sampleMult:
                for days in simDays:
                    for dem in globalDem:
                        l = []
                        l.extend([se,sp,mult,days,dem])
                        l = [str(l[i]) for i in range(len(l))]
                        wrtr = csv.writer(f,lineterminator = '\n')
                        wrtr.writerow(l)
                        #f.write('\t'.join(l) + '\n')
                        

