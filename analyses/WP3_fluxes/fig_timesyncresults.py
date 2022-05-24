# -*- coding: utf-8 -*-
"""
Created on Sun May 22 11:10:48 2022

@author: Dean

Summary plot of time synchronization results.
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Load and merge time sync and key altitudes tables:
tsyncresults = pd.read_csv("./other_data/timesync_wind-water_results.csv")
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
data = tsyncresults.merge(keyalts_table, on='ncld', how='left')


# Split data into 3 categories:
#       1) Trade inversion idenifiable and leg is below it.
#       2) Trade inversion idenifiable and leg is above it.
#       3) Trade inversion is not idenifiable.
below_invbottom = data['alt_leg'] < data['z_tib']
invbottom_exists = data['z_tib'].notnull()

data_c1 = data.loc[below_invbottom & invbottom_exists]

data_c2 = data.loc[~below_invbottom & invbottom_exists]

data_c3 = data.loc[~invbottom_exists]


# Scatter, colored by cross-correlation value:
fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.2, 0.15, 0.525, 0.8])
cax = fig.add_axes([0.825, 0.15, 0.05, 0.6])

sc = ax.scatter(
    data_c1['tshift'], data_c1['alt_leg'], 
    marker='o', s=14,
    c=data_c1['xcormax'], cmap='seismic', vmin=-0.2, vmax=0.4
    )
ax.scatter(
    data_c2['tshift'], data_c2['alt_leg'], 
    marker='*', s=14,
    c=data_c2['xcormax'], cmap='seismic', vmin=-0.2, vmax=0.4
    )
ax.scatter(
    data_c3['tshift'], data_c3['alt_leg'], 
    marker='^', s=14,
    c=data_c3['xcormax'], cmap='seismic', vmin=-0.2, vmax=0.4
    )
fig.colorbar(sc, cax=cax)


# Legend, labels, limits, etc:
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='below TIB',
           markerfacecolor='k', markersize=9),
    Line2D([0], [0], marker='*', color='w', label='above TIB',
           markerfacecolor='k', markersize=9),
    Line2D([0], [0], marker='^', color='w', label='TIB N/A',
           markerfacecolor='k', markersize=9),
    ]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    
ax.set_xlabel("time shift (s)", fontsize=12)
ax.set_ylabel("altitude (m)", fontsize=12)
cax.set_title("correlation", fontsize=12)


fig.savefig("./fig_timesyncresults.png")









