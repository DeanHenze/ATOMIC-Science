# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:47:26 2022

@author: Dean

Figure for vertical profiles of <w'w'>, sensible and latent heat fluxes, 
and dD of the flux. 
"""



import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import thermo
import rangescaler



## Data directories and filenames.
##_____________________________________________________________________________    
dir_flux = "./fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]


dir_prf = "../../data/WP3/cloud_modules/"
fnames_prf = [f for f in os.listdir(dir_prf) if f.endswith('.nc')]
##_____________________________________________________________________________
## Data directories and filenames.



def get_fluxprofiles(ncld_list, keyalts_table, dir_flux, 
                     scale_altkeys=[]):
    """
    Inputs
    ------
    keyalts_table: pandas.DataFrame.
        If scaling altitude by a subset of columns in this table, those 
        columns cannot have NAN's.
        
    scale_altkeys: iterable of str's.
        Optional. Keys in keyalts_table to scale altitude by. If not empty,
        altitudes will be mapped as:
            0 <-- 0
            1 <-- 1st key in scale_alt
            2 <-- 2nd key in scale_alt
            3 ...
    """

    keyalts_table = keyalts_table.set_index('ncld')


    # All flux filenames:
    fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]


    # Collect flux profiles into a dictionary of pandas.DataFrame's.
    # Each DataFrame will contain all profiles for a specific var.
    fluxvarkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE", 
                   "flux_sh", "flux_lh", "flux_b", "dD_flux"
                   ]
    fluxprfs_dict = {}
    for k in fluxvarkeys: fluxprfs_dict[k] = pd.DataFrame({})
    
    
    for n in ncld_list:
        
        print("Working on cloud module %i." % n)
        
        # Key altitude info for this module:
        keyalts = keyalts_table.loc[n]
        
        # Load flux data:
        fname_flux = [f for f in fnames_flux if "_cld%s_" % str(n).zfill(2) in f]
        if len(fname_flux) == 0: continue
        fname_flux = fname_flux[0]        
        fluxes = xr.load_dataset(dir_flux + fname_flux)
    
        # Optional scale altitude.
        # 0=sea level.
        if len(scale_altkeys) != 0:
            alt_scalepoints = [0] + [keyalts[k] for k in scale_altkeys]
            #alt = rangescaler.piecewise_linscale(
            #    fluxes['alt'].values, 
            #    (0, keyalts['z_lcl'], keyalts['z_ctmean_50p95p']), 
            #    (0,1,2)
            #    )
            alt = rangescaler.piecewise_linscale(
                fluxes['alt'].values, 
                alt_scalepoints, np.arange(len(alt_scalepoints))
                )
        else:
            alt = fluxes['alt'].values
            
        # Merge flux profiles of current module into resp. DataFrame's:
        for k in fluxvarkeys:
            
            fluxvar_pd = pd.DataFrame(  # easier to work in pandas.
                {k+"_cld%i"%n: fluxes[k].values}, 
                index=alt
                )
                
            fluxprfs_dict[k] = pd.merge(
                fluxprfs_dict[k], fluxvar_pd, 
                how='outer', left_index=True, right_index=True
                )   
            
            
    return fluxprfs_dict



def plot_fluxprofiles(fluxprfs_dict, varkeysplot, axset, 
                      plotmeans=False, pcolor='grey'):
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
        
    plotmeans: bool.
        If set to True, plot the average proifiles as well as individuals.
    """
    for varkey, ax in zip(varkeysplot, axset):

        # Individual profiles:
        for k in fluxprfs_dict[varkey].columns:
            colplot = fluxprfs_dict[varkey][k].dropna()
            ax.plot(colplot, colplot.index, c=pcolor, alpha=0.5)
            ax.scatter(colplot, colplot.index, c=pcolor, s=10, alpha=0.5)  
        
        # Optional compute and plot mean profile:
        if plotmeans:
            #altgrouped = np.round(fluxprfs_dict[varkey].index/0.5)*0.5 # vertical binning.
            altgrouped = fluxprfs_dict[varkey].index.values//0.25
            altgrouped[altgrouped % 2 == 0] += 1
            altgrouped = altgrouped*0.25
            fluxvar_grouped = fluxprfs_dict[varkey].groupby(altgrouped, axis=0, as_index=True)
            meanprf = fluxvar_grouped.mean().mean(axis=1)
            ax.plot(
                meanprf.values, meanprf.index, 
                color=pcolor, linestyle='-', linewidth=4
                )
        


## Plot P-3 flux profiles for wind components and TKE.
##_____________________________________________________________________________
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
#ncld_group1 = [1,6,7,9,10,11]
#ncld_group2 = [4,5,8,13]
ncld_group1 = [1, 7, 5, 9, 4, 11, 10, 12, 6]
ncld_group2 = [8, 15, 3, 2, 13, 16, 14]
scale_altkeys = ["z_lcl", "z_ctmean_50p95p"]
fluxprfs_g1 = get_fluxprofiles(
    ncld_group1, keyalts_table, dir_flux, scale_altkeys=scale_altkeys)
fluxprfs_g2 = get_fluxprofiles(
    ncld_group2, keyalts_table, dir_flux, scale_altkeys=scale_altkeys)

fig_wind, axes_wind = plt.subplots(1, 4, figsize=(10, 5))
windkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE"]
plot_fluxprofiles(fluxprfs_g1, windkeys, axes_wind, plotmeans=True, pcolor='grey')
plot_fluxprofiles(
    fluxprfs_g2, windkeys, axes_wind, plotmeans=True, pcolor='blue')
#axes_wind[0].legend()
##_____________________________________________________________________________
## Plot P-3 flux profiles for wind components and TKE.
    
   
    
## Plot P-3 flux profiles for SHF, LHF, and dD of flux.
##_____________________________________________________________________________
fig_scalar, axes_scalar = plt.subplots(1, 4, figsize=(10, 5))
scalarfluxkeys = ["flux_sh", "flux_lh", "flux_b", "dD_flux"]
plot_fluxprofiles(fluxprfs_g1, scalarfluxkeys, axes_scalar, plotmeans=True, pcolor='grey')
plot_fluxprofiles(fluxprfs_g2, scalarfluxkeys, axes_scalar, plotmeans=True, pcolor='blue')
##_____________________________________________________________________________
## Plot P-3 flux profiles for SHF, LHF, and dD of flux.
   

    
## Add surface flux mean +/- std for RHB measurements for the time period 
## 2020-01-17 to 2020-02-12.
##_____________________________________________________________________________
dir_rhbflux = "../../data/RHB/metflux/"
fname = "EUREC4A_ATOMIC_RonBrown_10min_nav_met_sea_flux_20200109-20200212_v1.3.nc"

rhb_metflux = xr.load_dataset(dir_rhbflux+fname)
dts_wp3start = np.datetime64("2020-01-17")
    # Remove measurements for dates before the P-3 started sampling:
rhb_metflux = rhb_metflux.where(rhb_metflux['time'] > dts_wp3start, drop=True)

rhbfluxkeys = ["hs_bulk", "hl_bulk", "ustar", "gust"]
flux_P3iopmean = rhb_metflux[rhbfluxkeys].mean(dim='obs')
flux_P3iopstd = rhb_metflux[rhbfluxkeys].std(dim='obs')

rhb_samplingalt = 20/500 # height of RHB instruments, roughly normalized by z_ML.
axes_scalar[0].errorbar(
    -1*flux_P3iopmean['hs_bulk'], rhb_samplingalt, 
    xerr = flux_P3iopstd['hs_bulk'], 
    marker='o', markersize=5, markerfacecolor='red', ecolor='red'
    )
axes_scalar[1].errorbar(
    -1*flux_P3iopmean['hl_bulk'], rhb_samplingalt, 
    xerr = flux_P3iopstd['hl_bulk'], 
    marker='o', markersize=5, markerfacecolor='red', ecolor='red'
    )
##_____________________________________________________________________________
  ## Add surface flux mean +/- std for RHB measurements.
    


## Axes limits, labels, etc for figure 1:    
for ax in axes_wind[0:3]: ax.set_xlim(0, 0.8)
yticks = [0,1,2,3,4]
for ax in axes_wind: ax.set_yticks(yticks)
yticklabels = [
    "0", r"$z_{ML}$", r"$z_{IB} = z_{ML} + \Delta z_{CL}$", 
    r"$z_{ML} + 2\Delta z_{CL}$", r"$z_{ML} + 3\Delta z_{CL}$"
    ]
axes_wind[0].set_yticklabels(yticklabels)
for ax in axes_wind[1:]: ax.set_yticklabels(["" for t in yticks])

axes_wind[0].set_xlabel(r"u'u' (m/s)$^2$", fontsize=14)
axes_wind[1].set_xlabel(r"v'v' (m/s)$^2$", fontsize=14)
axes_wind[2].set_xlabel(r"w'w' (m/s)$^2$", fontsize=14)
axes_wind[3].set_xlabel(r"TKE (m/s)$^2$", fontsize=14)



## Axes limits, labels, etc for figure 2:    
axes_scalar[1].set_xlim(-100, 650)
axes_scalar[3].set_xlim(-160, -30)
for ax in axes_scalar: ax.set_yticks(yticks)
axes_scalar[0].set_yticklabels(yticklabels)
for ax in axes_scalar[1:]: ax.set_yticklabels(["" for t in yticks])
axes_scalar[0].set_xlabel(r"$F_{sh}$ (W/m$^2$)", fontsize=14)
axes_scalar[1].set_xlabel(r"$F_{lh}$ (W/m$^2$)", fontsize=14)
axes_scalar[2].set_xlabel(r"$F_{b}$ (W/m$^2$)", fontsize=14)
axes_scalar[3].set_xlabel(r"$\delta D_{flux}$ (permil)", fontsize=14)



## Save figure:
fig_wind.savefig("./fig_fluxprofiles_windvars.png")
fig_scalar.savefig("./fig_fluxprofiles_scalarvars.png")


"""
## Night flights:
## Plot P-3 flux profiles for wind components and TKE.
##_____________________________________________________________________________
#ncld_list = [2,3,16]
ncld_list = [2,3,14,15]
flux_prfs = get_fluxprofiles(ncld_list, keyalts_table, dir_flux)

fig3, axes3 = plt.subplots(1, 4, figsize=(10, 5))
windkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE"]
plot_fluxprofiles(flux_prfs, windkeys, axes3)
##_____________________________________________________________________________
## Plot P-3 flux profiles for wind components and TKE.
    
   
    
## Plot P-3 flux profiles for SHF, LHF, and dD of flux.
##_____________________________________________________________________________
fig4, axes4 = plt.subplots(1, 3, figsize=(10, 5))
scalarfluxkeys = ["flux_sh", "flux_lh", "dD_flux"]
plot_fluxprofiles(flux_prfs, scalarfluxkeys, axes4)
##_____________________________________________________________________________
## Plot P-3 flux profiles for SHF, LHF, and dD of flux.
"""


