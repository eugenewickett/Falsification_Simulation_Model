# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 19:42:01 2020

@author: eugen
"""


sens = [0.75,0.85,0.95,0.99]
spec = [0.75,0.85,0.95,0.99]
sampleMult = [1,5,10]
simDays = [100,400,700]
globalDem = [0.,35.,70.]

with open('parameterFile.txt', 'w') as f:
    for se in sens:
        for sp in spec:
            for mult in sampleMult:
                for days in simDays:
                    for dem in globalDem:
                        l = []
                        l.extend([se,sp,mult,days,dem])
                        l = [str(l[i]) for i in range(len(l))]
                        f.write('\t'.join(l) + '\n')
                        

