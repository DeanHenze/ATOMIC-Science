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


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Load key altitudes table:
keyalts = pd.read_csv("./cldmod_keyaltitudes.csv")


# Plot:
fig = plt.figure(figsize=(6.5, 4))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])

pltmarkers = dict(
    z_mltop=dict(marker='^', color='green'), 
    z_lcl=dict(marker='o', facecolor='none', edgecolor='k'), 
    z_tib=dict(marker='s', facecolor='none', edgecolor='k'),
    z_tit=dict(marker='x', facecolor='r', edgecolor='r'), 
    z_ct=dict(marker='P', facecolor='y', edgecolor='y')
    )
    
for varkey in pltmarkers.keys():
    ax.scatter(
        keyalts['ncld'], keyalts[varkey], 
        **pltmarkers[varkey]
        )

# Figure legend, limits, labels:
ax.set_ylim(0, 3300)
ax.set_yticks(np.arange(0, 3500, 500))
ax.set_ylabel("altitude (m)", fontsize=12)

ax.set_xticks(np.arange(0, 16.1, 1))
ax.set_xlabel("cloud module #", fontsize=12)

legendmarkers = dict(
    z_mltop=dict(marker='^', color='green', linestyle=''), 
    z_lcl=dict(marker='o', markerfacecolor='none', 
               markeredgecolor='k', linestyle=''), 
    z_tib=dict(marker='s', markerfacecolor='none', 
               markeredgecolor='k', linestyle=''),
    z_tit=dict(marker='x', markerfacecolor='r', 
               markeredgecolor='r', linestyle=''), 
    z_ct=dict(marker='P', markerfacecolor='y', 
              markeredgecolor='y', linestyle='')
    )
legend_elements = [
    Line2D([0], [0], label=varkey,**legendmarkers[varkey], markersize=9)
    for varkey in legendmarkers.keys()
    ]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, ncol=5)












