# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 18:41:32 2020

@author: eugen
"""
import matplotlib.pyplot as plt
import numpy as np

xList = list(range(1,100))
x = np.array([i/100 for i in xList])

y = x-3*x**2+2*x**3

plt.plot(x,y)


y = x-x**2

plt.plot(x,y)
