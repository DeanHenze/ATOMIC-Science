# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:15:49 2022

@author: Dean

Individual and mean profiles of vertical velocity skewness for P-3 cloud 
modules.
"""


import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import thermo
import rangescaler


# Cloud groups
ncld_g1 = [7, 9, 11, 10, 12, 6]
ncld_g2 = [15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4, 8]



## Load and merge data tables:
wskew = pd.read_csv("./wskew_levlegs/WP3_wskewness_levlegs.csv")
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
data = wskew.merge(
    keyalts_table, 
    left_on="ncld", right_on="ncld", how="left"
    )


# Restructure data into a dataframe where each column is a separate profile 
# and the index is scaled altitude:
def restruct_wksewprofiles(data, ncld_list):
    """
    """
    wskew_prfs = pd.DataFrame({})
    
    for ncld, data_cld in data.groupby(by='ncld'):
    
        if ncld in ncld_list:
            
            # Scale altitude.
            # 0=sea level, 1=LCL, 2=trade inversion bottom or cloud top.
            z1 = data_cld['z_lcl'].iloc[0]
            z2 = data_cld['z_tib'].iloc[0]
            data_cld['altleg_scaled'] = rangescaler.piecewise_linscale(
                    data_cld['altleg'].values, 
                    (0, z1, z2), 
                    (0,1,2)
                    )
            
            # Append:
            data_cld.set_index('altleg_scaled', inplace=True)
            data_cld = data_cld.rename(columns={"wskew": "wskew_cld%i"%ncld})
            wskew_prfs = pd.merge(
                    wskew_prfs, data_cld[["wskew_cld%i"%ncld]], 
                    #how='outer', left_on='altleg_scaled', right_on='altleg_scaled'
                    how='outer', left_index=True, right_index=True
                    )   
            
    return wskew_prfs
        


def plot_wskewprfs(wskew_prfs, ax, color):
    """
    """
    # Plot individual profiles:
    for k in wskew_prfs.columns:
        colplot = wskew_prfs[k].dropna()
        ax.plot(colplot, colplot.index, c=color, alpha=0.5)
        ax.scatter(colplot, colplot.index, c=color, s=10, alpha=0.5)  
        
    # Plot mean profile:
    altgrouped = np.round(wskew_prfs.index/0.5)*0.5 # vertical binning.
    wskew_grouped = wskew_prfs.groupby(altgrouped, axis=0, as_index=True)
    meanprf = wskew_grouped.mean().mean(axis=1)
    ax.plot(meanprf.values, meanprf.index, c=color, linewidth=5)
    ax.scatter(meanprf.values, meanprf.index, c=color, s=10)
    
    
    
# Get wskew profiles for each cloud group and plot:
fig = plt.figure(figsize=(4,8))
ax = fig.add_axes([0.3, 0.1, 0.65, 0.8])

pltcolors = ['grey', 'blue', 'red']
for ncld_list, c in zip([ncld_g1, ncld_g2, ncld_g3], pltcolors):
    wskew_prfs = restruct_wksewprofiles(data, ncld_list)      
    plot_wskewprfs(wskew_prfs, ax, c)

# Vertical line at wskew=0 for visual reference:
ax.plot([0,0], [0,4], 'k-', linewidth=1)


# Figure limits, labels, etc:
ax.set_xlim(-1,2.5)
ax.set_ylim(0,3.5)
yticks = [0,1,2,3]
ax.set_yticks(yticks)
yticklabels = [
    "0", r"$z_{ML}$", r"$z_{IB} = z_{ML} + \Delta z_{CL}$", 
    r"$z_{ML} + 2\Delta z_{CL}$"
    ]
ax.set_yticklabels(yticklabels)
ax.set_xlabel(r"<w'w'w'>", fontsize=12)


fig.savefig("./fig_wskewprofiles.png")


















