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
from scipy import interpolate

# Local code
import thermo
import rangescaler



def plotprf_singlevar(data, ax, altbinwidth=0.33, npts_thresh=2, pcolor='grey'):
    """
    Plot a summary profile for a collection of individual profiles. Plot:
        - mean (solid line)
        - median (dashed)
        - min / max values at each altitude (shaded region)
    
    Inputs
    ------
    data: pandas.DataFrame.
        Profile data for one of the variables. The index is the altitude and 
        each column corresponds to one of the individual profiles.
                
    ax: matplotlib.pyplot.Axes.
        Axes to plot on.
        
    altbinwidth: scalar.
        Bin / average data by altitude bins of this width.
        
    npts_thresh: scalar.
        Minimum number of points in an altitude bin, otherwise data are 
        discarded.
        
    pcolor: str.
        Color to pass to matplotlib.
    """
    
    # Group by altitude bins with pandas:
    #altgrouped = np.round(data.index/0.25)*0.25 # vertical binning.
    altgrouped = np.round(data.index/altbinwidth)*altbinwidth # vertical binning.
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


    # Remove and levels with less than data point threshold:            
    ltnpoints = (fluxvar_grouped.count().sum(axis=1) < npts_thresh).values
    alt_bincenter = alt_bincenter[~ltnpoints]
    meanprf = meanprf[~ltnpoints]
    medianprf = medianprf[~ltnpoints]
    minvals = minvals[~ltnpoints]
    maxvals = maxvals[~ltnpoints]
    
    
    # Spline interpolation:
    alts_interp = np.arange(min(alt_bincenter), max(alt_bincenter), altbinwidth/1.5)
    alts_interp = np.append(alts_interp, max(alt_bincenter))
    meanprf_interp = interpolate.interp1d(alt_bincenter, meanprf, kind='cubic')(alts_interp)
    medianprf_interp = interpolate.interp1d(alt_bincenter, medianprf, kind='cubic')(alts_interp)
    minvals_interp = interpolate.interp1d(alt_bincenter, minvals, kind='cubic')(alts_interp)
    maxvals_interp = interpolate.interp1d(alt_bincenter, maxvals, kind='cubic')(alts_interp)


    # Plot:
    #ax.fill_betweenx(
    #    alt_bincenter, minvals, x2=maxvals, 
    #    color=pcolor, alpha=0.3, 
    #    )
    #ax.plot(
    #    meanprf, alt_bincenter, 
    #    color=pcolor, linestyle='-', linewidth=3, zorder=10
    #    )
    #ax.plot(
    #    medianprf, alt_bincenter, 
    #    color=pcolor, linestyle='--', linewidth=2, zorder=10
    #    )
    
    ax.fill_betweenx(
        alts_interp, minvals_interp, x2=maxvals_interp, 
        color=pcolor, alpha=0.3, 
        )
    ax.plot(
        meanprf_interp, alts_interp, 
        color=pcolor, linestyle='-', linewidth=3, zorder=10
        )
    ax.plot(
        medianprf_interp, alts_interp, 
        color=pcolor, linestyle='--', linewidth=2, zorder=10
        )
 
    