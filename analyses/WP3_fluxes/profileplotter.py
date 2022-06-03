# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 15:49:15 2022

@author: Dean

Plot vertical profiles with various options:
    - plot mean, median profile
    - plot profiles with scaled altitude
    - plt shaded regions for std's or min/max values at each altitude.
"""



import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import thermo
import rangescaler



def plotprf_singlevar(data, ax, pcolor='grey'):
    """
    Plot both individual profiles and mean profile for passed data.
    
    Inputs
    ------
    flux_prfs_dict: dictionary of pandas.DataFrame's.
        Each DataFrame contains profile data for one of the variables. The 
        index of the DataFrames is altitude and each column corresponds to one 
        of the profiles.
        
    varkeysplot: list of str's.
        Keys in flux_prfs_dict to plot.
        
    axset: list of matplotlib.pyplot.Axes.
        Same length as varkeysplot. Axes to plot on.
    """
    
    # Group by altitude bins with pandas:
    #altgrouped = np.round(data.index/0.25)*0.25 # vertical binning.
    altgrouped = np.round(data.index/0.33)*0.33 # vertical binning.
    fluxvar_grouped = data.groupby(altgrouped, axis=0, as_index=True)
    
    
    # For each group get mean, median, max, min:
    alt_bincenter = []
    meanprf = []
    medianprf = []
    minvals, maxvals = [], []
    for altbc, grp in fluxvar_grouped:
        alt_bincenter.append(altbc)
        grpvals_1d = grp.values.flatten()
        meanprf.append(np.nanmean(grpvals_1d))
        medianprf.append(np.nanmedian(grpvals_1d))
        minvals.append(np.nanmin(grpvals_1d))
        maxvals.append(np.nanmax(grpvals_1d))
    meanprf = np.array(meanprf)
    medianprf = np.array(medianprf)
    minvals = np.array(minvals)
    maxvals = np.array(maxvals)
    alt_bincenter = np.array(alt_bincenter)


    # Remove and levels with less than 2 data points:            
    lessthan2points = (fluxvar_grouped.count().sum(axis=1) < 2).values
    alt_bincenter = alt_bincenter[~lessthan2points]
    meanprf = meanprf[~lessthan2points]
    medianprf = medianprf[~lessthan2points]
    minvals = minvals[~lessthan2points]
    maxvals = maxvals[~lessthan2points]


    # Plot:
    ax.fill_betweenx(
        alt_bincenter, minvals, x2=maxvals, 
        color=pcolor, edgecolor='none', alpha=0.3
        )
    ax.plot(
        meanprf, alt_bincenter, 
        color=pcolor, linestyle='-', linewidth=3, zorder=10
        )
    ax.plot(
        medianprf, alt_bincenter, 
        color=pcolor, linestyle='--', linewidth=2, zorder=10
        )