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


### now for UNTRACKED
Q_ab = sym.Symbol('Q_ab')
z = pi+(1-pi)*th*Q_ab
pTilde = s*z + (1-r)*(1-z)
#from the outlet perspective
loglklhd = y*sym.log(pTilde)+(nSam-y)*sym.log(1-pTilde)
l_jac_betaA = sym.diff(loglklhd,betaA)

#for getting partials WRT betaB,betaA
import Falsification_Sim_Modules as simMod
betaA = sym.Symbol('betaA')
s = sym.Symbol('s')
r = sym.Symbol('r')
betaB = sym.Symbol('betaB')
y = sym.Symbol('y')
nSam = sym.Symbol('n')
th = sym.exp(betaB)/(sym.exp(betaB)+1)
pi = sym.exp(betaA)/(sym.exp(betaA)+1)
c1 = sym.Symbol('c1') #bunch of stuff indep. of betaA
c2 = sym.Symbol('c2') #summation
z = pi + (1-pi)*c2
pT = s*z + (1-r)*(1-z)

l_jac_betaB = c1*(1-pi)*(y/pT - (nSam-y)/(1-pT))

l_hess_betaB_betaA = sym.diff(l_jac_betaB,betaA)

bA,bB = beta0[3],beta0[0]
posits,numSamples = Y[0],N[0]
c1temp = Q[0,0]*(Sens+Spec-1)*exp(bB)/((exp(bB)+1)**2)
c2temp = np.sum([Q[0,j]*simMod.invlogit(beta0[j])  for j in range(m)])
repl = [(s,Sens),(r,Spec),(betaA,bA),(betaB,bB),(y,posits),(nSam,numSamples),\
        (c1,c1temp),(c2,c2temp)]

l_hess_betaB_betaA.subs(repl)

#for getting partials WRT betaB,betaB'
import Falsification_Sim_Modules as simMod
betaA = sym.Symbol('betaA')
s = sym.Symbol('s')
r = sym.Symbol('r')
betaB = sym.Symbol('betaB')
betaBpr = sym.Symbol('betaBpr')
y = sym.Symbol('y')
nSam = sym.Symbol('n')
th = sym.exp(betaB)/(sym.exp(betaB)+1)
thpr = sym.exp(betaBpr)/(sym.exp(betaBpr)+1)
pi = sym.exp(betaA)/(sym.exp(betaA)+1)
c1 = sym.Symbol('c1') #bunch of stuff indep. of betaBpr
c2 = sym.Symbol('c2')
c3 = sym.Symbol('c3')
z = c2+c3*thpr
pT = s*z + (1-r)*(1-z)

l_jac_betaB = c1*(y/pT - (nSam-y)/(1-pT))

l_hess_betaB_betaBpr = sym.diff(l_jac_betaB,betaBpr)

partial=0
for i in range(n):
    bA,bB, bBpr = beta0[i+m],beta0[0], beta0[1]
    posits,numSamples = Y[i],N[i]
    c1temp = Q[i,0]*(Sens+Spec-1)*(1-simMod.invlogit(bA))*exp(bB)/((exp(bB)+1)**2)
    c2temp = np.sum([Q[i,j]*simMod.invlogit(beta0[j])  for j in range(m)]) - Q[i,1]*simMod.invlogit(bBpr)
    c3temp = (1-simMod.invlogit(bA))*Q[i,1]
    repl = [(s,Sens),(r,Spec),(betaA,bA),(betaB,bB),(betaBpr,bBpr),(y,posits),(nSam,numSamples),\
            (c1,c1temp),(c2,c2temp),(c3,c3temp)]
    summand = l_hess_betaB_betaBpr.subs(repl)
    partial += summand
partial = partial*-1
print(partial)



























