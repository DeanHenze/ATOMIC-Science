# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:35:59 2022

@author: Dean

'Time series' plot of P-3 cloud module key heights:
    - mixed layer top height
    - LCL
    - trade inversion top and bottom
    - cloud top height
    
Also include plot of symbols denoting a few qualitative regimes - e.g. weak 
trade inversion, multiple inversions, etc.
"""


import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm

# Local code
import thermo


# Load key altitudes table:
keyalts = pd.read_csv("./cldmod_keyaltitudes.csv")



fig = plt.figure(figsize=(6.5, 4))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])


ax.scatter(
    keyalts['ncld'], keyalts['z_mltop'], 
    marker='^', color='green'
    )
ax.scatter(
    keyalts['ncld'], keyalts['z_lcl'], 
    marker='o', facecolor='none', edgecolor='k'
    )
ax.scatter(
    keyalts['ncld'], keyalts['z_tib'], 
    marker='s', facecolor='none', edgecolor='k'
    )
ax.scatter(
    keyalts['ncld'], keyalts['z_tit'], 
    marker='x', facecolor='r', edgecolor='k'
    )
ax.set_ylim(0, 3000)














