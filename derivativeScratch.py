# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 19:48:23 2020

@author: eugen
"""

import sympy as sym
betaA = sym.Symbol('betaA')
s = sym.Symbol('s')
r = sym.Symbol('r')
betaB = sym.Symbol('betaB')
y = sym.Symbol('y')
nSam = sym.Symbol('n')
th = sym.exp(betaB)/(sym.exp(betaB)+1)
pi = sym.exp(betaA)/(sym.exp(betaA)+1)
z = pi+th-pi*th
pTilde = s*z + (1-r)*(1-z)

#from the outlet perspective
jac_betaA = (1-th) * (s+r-1) * sym.diff(pi,betaA)*\
            (y/pTilde - (nSam-y)/(1-pTilde))

hess_betaA_betaA = sym.diff(jac_betaA,betaA)


Sens,Spec = 0.95, 0.95

outletPart = 0
for j in range(m):
    bA, bB = -4.5,-4.5
    posits, numSamples = 5, 100
    repl = [(s,Sens),(r,Spec),(betaA,bA),(betaB,bB),(y,posits),(nSam,numSamples)]
    summand = hess_betaA_betaA.subs(repl)
    outletPart += summand



