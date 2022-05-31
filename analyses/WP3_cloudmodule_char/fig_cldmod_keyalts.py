# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:35:59 2022

@author: Dean

'Time series' plot of P-3 cloud module key heights:
    - mixed layer top height
    - LCL
    - trade inversion top and bottom
    - cloud top height
    
Saves figure to working directory.

To do
------
Also include plot of symbols denoting a few qualitative regimes - e.g. weak 
trade inversion, multiple inversions, etc.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Local code
import oversampler



def plot_keyalts(keyalts, ax, xplot):
    """
    keyalts: pandas.Dataframe
        Key altitudes data set.
    """
    
    # Plot cloud tops as 50th - 95th percentile ranges:
    ax.errorbar(
        xplot, keyalts['z_ctmean_50p95p'],
        yerr=(
            abs(keyalts['z_ct50p']-keyalts['z_ctmean_50p95p']), 
            abs(keyalts['z_ct95p']-keyalts['z_ctmean_50p95p'])
            ), 
        linestyle='', marker='_', color='black',
        elinewidth=8, ecolor='blue', alpha=0.7
        )
    
    
    # Plot remaining heights as single markers:
    pltmarkers = dict(
        z_mltop=dict(marker='^', color='green'), 
        z_lcl=dict(marker='o', facecolor='none', edgecolor='k'), 
        z_tib=dict(marker='s', facecolor='none', edgecolor='k'),
        z_tit=dict(marker='x', facecolor='r', edgecolor='r'), 
        #z_ct=dict(marker='P', facecolor='y', edgecolor='y')
        )
        
    for varkey in pltmarkers.keys():
        ax.scatter(
            xplot, keyalts[varkey], 
            **pltmarkers[varkey], zorder=10
            )
    
    
    # Figure limits, labels:
    ax.set_xlim(np.min(xplot)-1, np.max(xplot)+1)
    ax.set_xticks(xplot)
    ax.set_xlabel("cloud module #", fontsize=12)
    ax.set_xticklabels([str(n).zfill(2) for n in keyalts['ncld']])
    
    ax.set_ylim(0, 3300)
    ax.set_yticks(np.arange(0, 3500, 500))
    ax.set_ylabel("altitude (m)", fontsize=12)
    
    
    # Figure legend:
    legendmarkers = dict(
        z_mltop=dict(marker='^', color='green', linestyle=''), 
        z_lcl=dict(marker='o', markerfacecolor='none', 
                   markeredgecolor='k', linestyle=''), 
        z_tib=dict(marker='s', markerfacecolor='none', 
                   markeredgecolor='k', linestyle=''),
        z_tit=dict(marker='x', markerfacecolor='r', 
                   markeredgecolor='r', linestyle=''), 
        #z_ct=dict(marker='P', markerfacecolor='y', 
        #          markeredgecolor='y', linestyle='')
        )
    legend_elements = [
        Line2D([0], [0], label=varkey, **legendmarkers[varkey], markersize=9)
        for varkey in legendmarkers.keys()
        ]
    legend_elements.append(Patch(facecolor='blue', edgecolor='blue', label='z_ct'))
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9, ncol=5)



def addthetacontours(dir_drpsnd, ncld, varkey, fig, ax):
    """
    """
    fnames_drpsnd = [
        f for f in os.listdir(dir_drpsnd) 
        if f.startswith("p3cld_dropsondes") and f.endswith(".nc") 
        ]
    
    
    # Get mean profiles of chosen variable for each cloud module:
    meanvarprfs = []
    altvals = []
    xvals = []
    xplot = ax.get_xticks()
    
    for i in range(len(ncld)):
        nstr = str(ncld[i]).zfill(2)
        f = [f for f in fnames_drpsnd if "_ncld%s" % nstr in f][0]
        drpsnds = xr.load_dataset(dir_drpsnd + f)
        
        prf = drpsnds[varkey].mean(dim='sounding').dropna(dim='alt')
        meanvarprfs.append(prf.values)
        altvals.append(prf['alt'].values)
        xvals.append([xplot[i] for e in prf])

    
    # Get smoothed map via oversampling:
        # Prep lists for use with oversampler:
    listflatten = lambda l: [item for sublist in l for item in sublist]
    meanvarprfs = listflatten(meanvarprfs)
    altvals = listflatten(altvals)
    xvals = listflatten(xvals)
        # oversampler:
    yplot = np.arange(0, 3001, 10)
    varmap = oversampler.oversampler(
        meanvarprfs, xvals, altvals, 
        xplot, yplot, 
        w_add=None, ffact=0.3, return_stdev='no'
        )
    
    
    # Contour plot:
    levels  = np.round(np.arange(297, 312, 2))
    cf = ax.contourf(
        xplot, yplot, varmap['mean'],
        levels=levels, cmap='Greys'
        )
    c = ax.contour(
        xplot, yplot, varmap['mean'],
        levels=levels[1:4], colors='black', linewidths=0
        )
    fmt = dict(zip(levels, [r"$\theta=%i K$" % l for l in levels]))
    manuallocs=[(11, 500), (13, 1300), (12, 2000)]
    ax.clabel(
        c, levels=levels[1:4], fmt=fmt, 
        inline=True, manual=manuallocs, fontsize=8
        )
    #fig.colorbar(c, ax=ax)



if __name__=="__main__":
    
    keyalts = pd.read_csv("./cldmod_keyaltitudes.csv") # key altitudes table
    keyalts = keyalts.sort_values('z_ctmean_50p95p')
    
    dir_drpsnd = "./cldmod_datafiles/"
    
    fig = plt.figure(figsize=(6.5, 4))
    ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
    xplot = np.arange(len(keyalts.index)) # x-axis values for plotting
    ax.set_xticks(xplot)
    addthetacontours(dir_drpsnd, keyalts['ncld'], 'theta', fig, ax)
    plot_keyalts(keyalts, ax, xplot)

    fig.savefig("./fig_cloudmod_keyalts.png")










