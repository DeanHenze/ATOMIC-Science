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
path_keyaltstable = "../WP3_cloudmodule_char/cldmod_keyaltitudes.csv"

dir_flux = "./fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]

dir_prf = "../../data/WP3/cloud_modules/"
fnames_prf = [f for f in os.listdir(dir_prf) if f.endswith('.nc')]

path_rhbflux = ("../../data/RHB/metflux/EUREC4A_ATOMIC_RonBrown_10min"
                "_nav_met_sea_flux_20200109-20200212_v1.3.nc"
                )
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

        """
        # Individual profiles:
        for k in fluxprfs_dict[varkey].columns:
            colplot = fluxprfs_dict[varkey][k].dropna()
            ax.plot(colplot, colplot.index, c=pcolor, alpha=0.5)
            ax.scatter(colplot, colplot.index, c=pcolor, s=10, alpha=0.3)  
        """
        
        # Optional compute and plot mean profile:
        if plotmeans:
            # Group by altitude bins with pandas:
            #altgrouped = np.round(fluxprfs_dict[varkey].index/0.25)*0.25 # vertical binning.
            altgrouped = np.round(fluxprfs_dict[varkey].index/0.33)*0.33 # vertical binning.
            fluxvar_grouped = fluxprfs_dict[varkey].groupby(altgrouped, axis=0, as_index=True)
            
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
        


def plot_RHBmeanfluxes(path_rhbflux, ax_sh, ax_lh):
    """
    Plot mean, std of RHB surface fluxes during IOP portion with P-3 flights. 
    Plot SHF and LHF on ax_sh and ax_lh repspectively.
    """    
    rhb_metflux = xr.load_dataset(path_rhbflux)
    
    # Remove measurements for dates before the P-3 started sampling:
    dts_wp3start = np.datetime64("2020-01-17")
    rhb_metflux = rhb_metflux.where(rhb_metflux['time'] > dts_wp3start, drop=True)

    # Means and std's:
    rhbfluxkeys = ["hs_bulk", "hl_bulk", "ustar", "gust"]
    flux_P3iopmean = rhb_metflux[rhbfluxkeys].mean(dim='obs')
    flux_P3iopstd = rhb_metflux[rhbfluxkeys].std(dim='obs')
    
    # Error bar plot:
    rhb_altsamp = 20/500 # height of RHB instruments, roughly normalized by z_ML.
    ax_sh.errorbar(
        -1*flux_P3iopmean['hs_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hs_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black'
        )
    ax_lh.errorbar(
        -1*flux_P3iopmean['hl_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hl_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black'
        )



def xaxes_cleanup(axset_wind, axset_scalar):
    """
    Set x-axis labels, limits, etc. for sets of wind and flux axes.
    """

    ## Limits, labels, etc for wind / turbulence axes:    
    for ax in axset_wind[0:3]: ax.set_xlim(0, 0.8)    
    axset_wind[0].set_xlabel(r"u'u' (m/s)$^2$", fontsize=12)
    axset_wind[1].set_xlabel(r"v'v' (m/s)$^2$", fontsize=12)
    axset_wind[2].set_xlabel(r"w'w' (m/s)$^2$", fontsize=12)
    axset_wind[3].set_xlabel(r"TKE (m/s)$^2$", fontsize=12)

    ## Limits, labels, etc for flux axes:    
    axset_scalar[0].set_xlim(-30, 35)
    axset_scalar[1].set_xlim(-100, 650)
    axset_scalar[2].set_xlim(-5e-4, 1.7e-3)
    axset_scalar[3].set_xlim(-160, -30)
    axset_scalar[0].set_xlabel(r"$F_{sh}$ (W/m$^2$)", fontsize=12)
    axset_scalar[1].set_xlabel(r"$F_{lh}$ (W/m$^2$)", fontsize=12)
    axset_scalar[2].set_xlabel(r"$F_{b}$ ($m^2/s^3$)", fontsize=12)
    axset_scalar[3].set_xlabel(r"$\delta D_{flux}$ (permil)", fontsize=12)



def yaxes_cleanup(axset_wind, axset_scalar, yticks, yticklabels):
    """
    Set y-axis labels, limits, etc. for sets of wind and flux axes.
    """
    ## Limits, labels, etc for wind / turbulence axes:    
    for ax in axset_wind:
        ax.set_ylim(min(yticks), max(yticks))
        ax.set_yticks(yticks)
    axset_wind[0].set_yticklabels(yticklabels)
    for ax in axset_wind[1:]: ax.set_yticklabels(["" for t in yticks])
    
    ## Limits, labels, etc for flux axes:    
    for ax in axset_scalar:
        ax.set_ylim(min(yticks), max(yticks))
        ax.set_yticks(yticks)
    axset_scalar[0].set_yticklabels(yticklabels)
    for ax in axset_scalar[1:]: ax.set_yticklabels(["" for t in yticks])



def analyzeplotgroupings(ncld_groups, colors, scale_altkeys, keyalts_table, 
                         axset_wind, axset_scalar):
    """
    """
    # Get flux profiles grouped; groups are elements of a list:
    fluxprfs_grouped = []
    for cg in ncld_groups:
        fluxprfs_grouped.append(
            get_fluxprofiles(       # Returns flux profiles as a dictionary of 
                cg, keyalts_table,  # pandas.DataFrame's. 
                dir_flux, scale_altkeys=scale_altkeys
                )
            )
        
    # Plot wind / turbulence profiles with means overlain:
    windkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE"] # vars to plot.   
    for fpgroup, c in zip(fluxprfs_grouped, colors):
        plot_fluxprofiles(
            fpgroup, windkeys, axset_wind, 
            plotmeans=True, pcolor=c
            )

    # Plot scalar flux profiles with means overlain:
    scalarfluxkeys = ["flux_sh", "flux_lh", "flux_b", "dD_flux"]
    for fpgroup, c in zip(fluxprfs_grouped, colors):
        plot_fluxprofiles(
            fpgroup, scalarfluxkeys, axset_scalar, 
            plotmeans=True, pcolor=c
            )


def fig_LCLCTscaling():
    """
    Create and save figures for turbulence and flux profiles where altitude 
    is scaled by LCL and cloud top height.
    """
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    #ncld_g1 = [1, 7, 5, 9, 4, 11, 10, 12, 6]
    #ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
    ncld_g1 = [7, 9, 11, 10, 12, 6]
    ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
    ncld_g3 = [1, 5, 4]
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    
    scale_altkeys = ["z_lcl", "z_ctmean_50p95p"] # scale altitude by these quantities.

    fig_wind, axset_wind = plt.subplots(1, 4, figsize=(10, 5))
    fig_scalar, axset_scalar = plt.subplots(1, 4, figsize=(10, 5))
    
    analyzeplotgroupings(
        ncld_groups, ['grey', 'blue', 'red'], 
        scale_altkeys, keyalts_table, 
        axset_wind, axset_scalar
        )
    
    
    ## Add RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1])

    
    ## Axes limits, labels, etc.
    xaxes_cleanup(axset_wind, axset_scalar)

    yticks = [0,1,2,3]
    yticklabels = ["0", r"$z_{LCL}$", r"$z_{CT}$", 
                   r"$z_{LCL}$+2$\Delta z_{(CT-LCL)}$"
                   ]  
    yaxes_cleanup(axset_wind, axset_scalar, yticks, yticklabels)       
      

    ## Save figure:
    fig_wind.savefig("./fig_wind+turb_profiles_LCLCTscaling.png")
    fig_scalar.savefig("./fig_scalarflux_profiles_LCLTCTscaling.png")
    
    
    
def fig_LCLTIBscaling():
    """
    Create and save figures for turbulence and flux profiles where altitude 
    is scaled by LCL and trade inversion bottom.
    """
    
    ## P-3 flux profiles:
    keyalts_table = pd.read_csv(path_keyaltstable)
    
    ncld_group1 = [1, 7, 5, 9, 4, 11, 10, 12, 6]
    ncld_group2 = [8, 15, 3, 2, 13, 16, 14]
    scale_altkeys = ["z_lcl", "z_tib"]
    
    fluxprfs_g1 = get_fluxprofiles(
        ncld_group1, keyalts_table, dir_flux, scale_altkeys=scale_altkeys)
    fluxprfs_g2 = get_fluxprofiles(
        ncld_group2, keyalts_table, dir_flux, scale_altkeys=scale_altkeys)
    

    ## Plot profiles for wind components and TKE.
    fig_wind, axset_wind = plt.subplots(1, 4, figsize=(10, 5))
    windkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE"]    
    plot_fluxprofiles(fluxprfs_g1, windkeys, axset_wind, plotmeans=True, pcolor='grey')
    plot_fluxprofiles(
        fluxprfs_g2, windkeys, axset_wind, plotmeans=True, pcolor='blue')
        
        
    ## Plot profiles for SHF, LHF, buoyancy flux, and dD of flux.
    fig_scalar, axset_scalar = plt.subplots(1, 4, figsize=(10, 5))
    scalarfluxkeys = ["flux_sh", "flux_lh", "flux_b", "dD_flux"]
    plot_fluxprofiles(fluxprfs_g1, scalarfluxkeys, axset_scalar, plotmeans=True, pcolor='grey')
    plot_fluxprofiles(fluxprfs_g2, scalarfluxkeys, axset_scalar, plotmeans=True, pcolor='blue')

    
    ## Add RHB surface flux means, stds for the P-3 sampling time period:
    path_rhbflux = ("../../data/RHB/metflux/EUREC4A_ATOMIC_RonBrown_10min"
                    "_nav_met_sea_flux_20200109-20200212_v1.3.nc"
                    )
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1])

    
    ## Axes limits, labels, etc.
    xaxes_cleanup(axset_wind, axset_scalar)

    yticks = [0,1,2,3,4]
    yticklabels = [
        "0", r"$z_{LCL}$", r"$z_{IB}$", r"$z_{LCL} + 2\Delta z_{(IB-LCL)}$", 
        r"$z_{LCL} + 3\Delta z_{(IB-LCL)}$"
        ]  
    yaxes_cleanup(axset_wind, axset_scalar, yticks, yticklabels)       
      

    ## Save figure:
    fig_wind.savefig("./fig_wind+turb_profiles_LCLTIBscaling.png")
    fig_scalar.savefig("./fig_scalarflux_profiles_LCLTIBscaling.png")



if __name__=="__main__":
    fig_LCLCTscaling()
    #fig_LCLTIBscaling()