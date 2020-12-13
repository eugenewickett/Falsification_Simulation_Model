# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 23:55:58 2020

@author: eugen
"""

'''
This script reads an CSV list of PMS testing results and returns a list of
estimation intervals for SFP rates for the Linear, Untracked, and Tracked
methods, as well as posterior distribution samples for the Untracked and
Tracked methods.
'''

import numpy as np
import DRA_EstimationMethods as estMethods

### 
fileName = ''

