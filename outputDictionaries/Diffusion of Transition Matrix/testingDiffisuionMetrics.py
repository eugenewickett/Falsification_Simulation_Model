# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:20:01 2020

@author: eugen
"""

'''
We first generate some matrices to compare measurements:
    A0: no diffusion
    A1: diffuse very little
    A2: diffuse some
    A3: diffuse a lot
    A4: perfectly diffuse
'''
import numpy as np
import matplotlib.pyplot as plt

A0 = np.array([[1,0,0],
               [0,0,1],
               [0,1,0],
               [1,0,0],
               [0,1,0],
               [0,0,1]])
A1 = np.array([[0.9,0.1,0],
               [0.1,0,0.9],
               [0,0.9,0.1],
               [0.9,0,0.1],
               [0.1,0.9,0],
               [0,0.1,0.9]])
A2 = np.array([[0.7,0.2,0.1],
               [0.2,0.1,0.7],
               [0.1,0.7,0.2],
               [0.9,0.1,0.0],
               [0,0.9,0.1],
               [0,0.2,0.8]])
A3 = np.array([[0.5,0.3,0.2],
               [0.1,0.4,0.5],
               [0.2,0.6,0.2],
               [0.4,0.3,0.3],
               [0.1,0.5,0.4],
               [0.3,0.1,0.6]])
A4 = np.array([[0.33,0.33,0.34],
               [0.35,0.33,0.32],
               [0.32,0.36,0.32],
               [0.33,0.34,0.33],
               [0.34,0.35,0.31],
               [0.35,0.31,0.34]])

    
np.trace(A0.T @ A0)
np.trace(A0 @ A0.T)
np.linalg.det(A0.T @ A0)
np.linalg.det(A0 @ A0.T)

np.trace(A1.T @ A1)
np.trace(A1 @ A1.T)
np.linalg.det(A1.T @ A1)
np.linalg.det(A1 @ A1.T)

np.trace(A2.T @ A2)
np.trace(A2 @ A2.T)
np.linalg.det(A2.T @ A2)
np.linalg.det(A2 @ A2.T)

np.trace(A3.T @ A3)
np.trace(A3 @ A3.T)
np.linalg.det(A3.T @ A3)
np.linalg.det(A3 @ A3.T)

np.trace(A4.T @ A4)
np.trace(A4 @ A4.T)
np.linalg.det(A4.T @ A4)
np.linalg.det(A4 @ A4.T)

x = [np.trace(A0.T @ A0),np.trace(A1.T @ A1),np.trace(A2.T @ A2),\
     np.trace(A3.T @ A3),np.trace(A4.T @ A4)]

y = [np.linalg.det(A0.T @ A0),np.linalg.det(A1.T @ A1),\
     np.linalg.det(A2.T @ A2),np.linalg.det(A3.T @ A3),\
     np.linalg.det(A4.T @ A4)]

plt.plot(x,y)
plt.plot(x,np.log(y))
plt.plot(x,np.sqrt(y))




