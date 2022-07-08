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
import profileplotter
import rangescaler
import prfrestructure
import thermo_deSzoeke



## Data directories and filenames.
##_____________________________________________________________________________    
path_keyaltstable = "../WP3_cloudmodule_char/cldmod_keyaltitudes.csv"

dir_flux = "./fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]

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
    fluxvarkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE", "TKE_h",
                   "flux_sh", "flux_lh", "flux_b", "dD_flux", 
                   "q'w'_bar", "T'w'_bar"
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



def plot_fluxprofiles(fluxprfs_dict, varkeysplot, axset, pcolor='grey'):
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
    for varkey, ax in zip(varkeysplot, axset):        
        #profileplotter.plotprf_singlevar(fluxprfs_dict[varkey], ax, pcolor=pcolor)
        profileplotter.plotprf_singlevar(
            fluxprfs_dict[varkey], ax, pcolor=pcolor, 
            altbinwidth=0.5, npts_thresh=2, 
            cubic_interp=False
            )



def plot_RHBmeanfluxes(path_rhbflux, ax_sh, ax_lh, ax_bflux):
    """
    Plot mean, std of RHB surface fluxes during IOP portion with P-3 flights. 
    Plot SHF and LHF on ax_sh and ax_lh repspectively.
    """    
    rhb_metflux = xr.load_dataset(path_rhbflux)
    
    # Remove measurements for dates before the P-3 started sampling:
    dts_wp3start = np.datetime64("2020-01-17")
    rhb_metflux = rhb_metflux.where(rhb_metflux['time'] > dts_wp3start, drop=True)

    # Compute buoyancy flux:
    flux_b = thermo_deSzoeke.buoyancy_flux(
        rhb_metflux['hs_bulk'], rhb_metflux['hl_bulk'], 
        T=rhb_metflux['tair'], p=rhb_metflux["psealevel_ship"]*100, 
        q=rhb_metflux["qair"]/1000
        )
    rhb_metflux['bflux_bulk'] = flux_b[0]

    # Means and std's:
    rhbfluxkeys = ["hs_bulk", "hl_bulk", "bflux_bulk", "ustar", "gust"]
    flux_P3iopmean = rhb_metflux[rhbfluxkeys].mean(dim='obs')
    flux_P3iopstd = rhb_metflux[rhbfluxkeys].std(dim='obs')
    
    # Error bar plot:
    rhb_altsamp = 20/500 # height of RHB instruments, roughly normalized by z_ML.
    ax_sh.errorbar(
        -1*flux_P3iopmean['hs_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hs_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )
    ax_lh.errorbar(
        -1*flux_P3iopmean['hl_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hl_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )
    ax_bflux.errorbar(
        -1*flux_P3iopmean['bflux_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['bflux_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )



def xaxes_cleanup(axset_wind, axset_scalar):
    """
    Set x-axis labels, limits, etc. for sets of wind and flux axes.
    """

    ## Limits, labels, etc for wind / turbulence axes:    
    for ax in axset_wind[0:3]: 
        ax.set_xlim(0, 0.8)
        ax.set_xticks(np.arange(0, 1, 0.2))
    axset_wind[3].set_xlim(0, 1.2)  
    axset_wind[3].set_xticks(np.arange(0, 1.2, 0.2))
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
    axset_wind[0].set_ylabel('altitude (m)', fontsize=12)
    
    ## Limits, labels, etc for flux axes:    
    for ax in axset_scalar:
        ax.set_ylim(min(yticks), max(yticks))
        ax.set_yticks(yticks)
    axset_scalar[0].set_yticklabels(yticklabels)
    for ax in axset_scalar[1:]: ax.set_yticklabels(["" for t in yticks])
    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)



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
        plot_fluxprofiles(fpgroup, windkeys, axset_wind, pcolor=c)

    # Plot scalar flux profiles with means overlain:
    scalarfluxkeys = ["flux_sh", "flux_lh", "flux_b", "dD_flux"]
    for fpgroup, c in zip(fluxprfs_grouped, colors):
        plot_fluxprofiles(fpgroup, scalarfluxkeys, axset_scalar, pcolor=c)
        
        
        
def analyzeplotgroupings_nonscaled(ncld_groups, colors, keyalts_table,
                         axset_wind, axset_scalar):
    """
    """
    # Get flux profiles grouped; groups are elements of a list:
    fluxprfs_grouped = []
    for cg in ncld_groups:
        fluxprfs_grouped.append(
            get_fluxprofiles(       # Returns flux profiles as a dictionary of 
                cg, keyalts_table,  # pandas.DataFrame's. 
                dir_flux, scale_altkeys=[]
                )
            )
        
    
    def plotprfs(fluxprfs_dict, varkeysplot, axset, pcolor='grey'):
        for varkey, ax in zip(varkeysplot, axset):        
            profileplotter.plotprf_singlevar(
                fluxprfs_dict[varkey], ax, 
                altbinwidth=400, npts_thresh=1, cubic_interp=False, 
                pcolor=pcolor
                )
        
        
    # Plot wind / turbulence profiles with means overlain:
    windkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE"] # vars to plot.   
    for fpgroup, c in zip(fluxprfs_grouped, colors):
        plotprfs(fpgroup, windkeys, axset_wind, pcolor=c)

    # Plot scalar flux profiles with means overlain:
    scalarfluxkeys = ["flux_sh", "flux_lh", "flux_b", "dD_flux"]
    for fpgroup, c in zip(fluxprfs_grouped, colors):
        plotprfs(fpgroup, scalarfluxkeys, axset_scalar, pcolor=c)
        
        
    xaxes_cleanup(axset_wind, axset_scalar)
    yticks = np.arange(0, 3001, 500)
    yaxes_cleanup(
        axset_wind, axset_scalar, 
        yticks, yticks.astype(str)
        )
    for ax in axset_wind.flatten(): ax.set_ylim(-100, 3200)
    for ax in axset_scalar.flatten(): ax.set_ylim(-100, 3200)
        


def fig_noscaling():
    """
    Create and save figures for turbulence and flux profiles with altitude 
    in meters.
    """
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    ncld_g1 = [7, 9, 11, 10, 12, 6]
    ncld_g2 = [15, 3, 2, 13, 16, 14]
    ncld_g3 = [1, 5, 4, 8]
    
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    

    # Plot:
    fig_wind, axset_wind = plt.subplots(1, 4, figsize=(10, 5))
    fig_scalar, axset_scalar = plt.subplots(1, 4, figsize=(10, 5))
        # Profiles:
    analyzeplotgroupings_nonscaled(
        ncld_groups, ['grey', 'blue', 'red'], 
        keyalts_table, 
        axset_wind, axset_scalar
        )
        # RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1])
        # Axes limits, labels, etc.
    xaxes_cleanup(axset_wind, axset_scalar)


    ## Save figures:
    fig_wind.savefig("./fig_wind+turb_profiles.png")
    fig_scalar.savefig("./fig_scalarflux_profiles.png")



def fig_LCLCTscaling():
    """
    Create and save figures for turbulence and flux profiles where altitude 
    is scaled by LCL and cloud top height.
    """
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    ncld_g1 = [7, 9, 11, 10, 12, 6]
    ncld_g2 = [15, 3, 2, 13, 16, 14]
    #ncld_g1 = [7, 9, 11, 10, 12, 6]
    #ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
    ncld_g3 = [1, 5, 4, 8]
    
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    
    scale_altkeys = ["z_lcl", "z_ctmean_50p95p"] # scale altitude by these quantities.


    # Plot:
    fig_wind, axset_wind = plt.subplots(1, 4, figsize=(10, 5))
    fig_scalar, axset_scalar = plt.subplots(1, 4, figsize=(10, 5))
        # Profiles:
    analyzeplotgroupings(
        ncld_groups, ['grey', 'blue', 'red'], 
        scale_altkeys, keyalts_table, 
        axset_wind, axset_scalar
        )
        # RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1])
        # Axes limits, labels, etc.
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
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    ncld_g1 = [7, 9, 11, 10, 12, 6]
    ncld_g2 = [15, 3, 2, 13, 16, 14]
    #ncld_g1 = [7, 9, 11, 10, 12, 6]
    #ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
    ncld_g3 = [1, 5, 4, 8]
    
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    
    #scale_altkeys = ["z_lcl", "z_madev"] # scale altitude by these quantities.
    scale_altkeys = ["z_lcl"] # scale altitude by these quantities.


    # Plot:
    fig_wind, axset_wind = plt.subplots(1, 4, figsize=(10, 5))
    fig_scalar, axset_scalar = plt.subplots(1, 4, figsize=(10, 5))
        
        # Profiles:
    analyzeplotgroupings(
        ncld_groups, ['grey', 'blue', 'red'], 
        scale_altkeys, keyalts_table, 
        axset_wind, axset_scalar
        )
        # RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1])
        # Axes limits, labels, etc.
    xaxes_cleanup(axset_wind, axset_scalar)
    yticks = [0,1,2,3]
    yticklabels = ["0", r"$z_{LCL}$", r"$z_{IB}$", 
                   r"$z_{LCL} + 2\Delta z_{(IB-LCL)}$", 
                   ]  
    yaxes_cleanup(axset_wind, axset_scalar, yticks, yticklabels)       
      

    ## Save figure:
    fig_wind.savefig("./fig_wind+turb_profiles_LCLTIBscaling.png")
    fig_scalar.savefig("./fig_scalarflux_profiles_LCLTIBscaling.png")
    
    

def fig_scatter():
    """
    Create and save figures for turbulence and flux profiles where altitude 
    is scaled by LCL and trade inversion bottom.
    """
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    #ncld_g1 = [1, 5, 4, 8]
    #ncld_g2 = [7, 9, 11, 10, 12, 6]
    #ncld_g3 = [15, 3, 2, 13, 16, 14]  
    ncld_g1 = [1, 5, 4]
    ncld_g2 = [8, 7, 9, 11, 10, 6]
    ncld_g3 = [15, 3, 2, 13, 16, 14, 12]  
    
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    
    #scale_altkeys = ["z_lcl"] # scale altitude by these quantities.
    #altbinwidth=0.5
    scale_altkeys = [] # scale altitude by these quantities.
    altbinwidth = 250

    # Get flux profiles grouped; groups are elements of a list:
    fluxprfs_grouped = []
    for cg in ncld_groups:
        fluxprfs_grouped.append(
            get_fluxprofiles(       # Returns flux profiles as a dictionary of 
                cg, keyalts_table,  # pandas.DataFrame's. 
                dir_flux, scale_altkeys=scale_altkeys
                )
            )
        
        
    ## Plot turbulence and flux profiles
    ##_________________________________________________________________________
    def scatterwithmean(data_prfs, ax, color, altbinwidth):
        """
        """
        for key_prf in data_prfs.columns:
            prf = data_prfs[key_prf]
            ax.scatter(prf, data_prfs.index, color=c, s=8)
        
        altgrouped = np.round(data_prfs.index/altbinwidth)*altbinwidth # vertical binning.
        fluxvar_grouped = data_prfs.groupby(altgrouped, axis=0, as_index=True)
        alt_bincenter = []
        meanprf = []
        for altbc, grp in fluxvar_grouped:
            alt_bincenter.append(altbc)
            grpvals_1d = grp.values.flatten()
            meanprf.append(np.nanmean(grpvals_1d))
        
        df = pd.DataFrame({'data_prfs':meanprf}, index=alt_bincenter)
        refalt = pd.DataFrame(
            index=pd.Index(np.arange(0, np.nanmax(data_prfs), altbinwidth)))
        df = df.merge(refalt, 
                      left_index=True, right_index=True, how='outer')
        #df = data_prfs.stack().reset_index()[['level_0', 0]]
        test = df.rolling(
            window=4, min_periods=1, win_type='hamming', 
            center=True, closed='neither'
            )
        test = test.mean()
        ax.plot(test.iloc[1:], test.index[1:], color=c, linewidth=2)
        #ax.plot(test[0].iloc[1:], test['level_0'][1:], color=c, linewidth=2)
        

    fig_scalar = plt.figure(figsize=(6.5, 3))
    #axset_scalar = (
    #    fig_scalar.add_axes([0.1, 0.2, 0.25, 0.75]),
    #    fig_scalar.add_axes([0.4, 0.2, 0.25, 0.75]),
    #    fig_scalar.add_axes([0.7, 0.2, 0.25, 0.75]),
    #    )
    axset_scalar = (
        fig_scalar.add_axes([0.1, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.325, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.55, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.775, 0.2, 0.2, 0.75]),
        )
    fig_wind = plt.figure(figsize=(6.5, 3))
    axset_wind = (
        fig_wind.add_axes([0.1, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.325, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.55, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.775, 0.2, 0.2, 0.75]),
        )
    
    
    colors = ['red', 'grey', 'blue']
    for fpgroup, c in zip(fluxprfs_grouped, colors):

        scatterwithmean(fpgroup["TKE_h"], axset_wind[0], c, altbinwidth)                  
        scatterwithmean(fpgroup["w'w'_bar"], axset_wind[1], c, altbinwidth)                  
        scatterwithmean(fpgroup["TKE"], axset_wind[2], c, altbinwidth)                  

        scatterwithmean(fpgroup["flux_sh"], axset_scalar[0], c, altbinwidth)        
        scatterwithmean(fpgroup["flux_lh"], axset_scalar[1], c, altbinwidth)        
        scatterwithmean(fpgroup["flux_b"], axset_scalar[2], c, altbinwidth)        
        scatterwithmean(fpgroup["dD_flux"], axset_scalar[3], c, altbinwidth)        

        
    # RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1], axset_scalar[2])
    
    
    # Ref line at x=0 for buoyancy flux:
    axset_scalar[2].vlines(
        0, -100, 3300, 
        color='black', linewidth=1.5, linestyle='dashed', zorder=100
        )
    ##_________________________________________________________________________    
    
    
    ## Plot w' skewness
    ##_________________________________________________________________________
    # Load and merge data tables:
    wmom = pd.read_csv("./WP3_wmoments_levlegs.csv")
    keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
    wmomdata = wmom.merge(
        keyalts_table, 
        left_on="ncld", right_on="ncld", how="left"
        )
    
    
    wmomprfs = []
    for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
        wmomprfs.append(
            prfrestructure.restruct_profiles(
                wmomdata, 'wskew', ncld_list, 
                keyalts_table, scale_altkeys=[]
                )
            )
        
    # Plot skewness profiles:
    for prfs, c in zip(wmomprfs, colors):  
            scatterwithmean(prfs, axset_wind[3], c, altbinwidth)
            
    # Ref line at x=0:
    axset_wind[3].vlines(
        0, -100, 3300, 
        color='black', linewidth=1.5, linestyle='dashed', zorder=100
        )
    ##_________________________________________________________________________

    
    axset_wind[0].set_yticks(axset_wind[0].get_yticks())
    axset_wind[0].set_yticklabels(['' for t in axset_wind[1].get_yticks()])
    axset_wind[0].set_ylim(-100, 3300)
    axset_wind[0].set_xlabel(r"TKE$_h$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[0].set_xlim(0, 1.2)
    axset_wind[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
    axset_wind[0].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

    axset_wind[1].set_yticks(axset_wind[0].get_yticks())
    axset_wind[1].set_ylim(-100, 3300) 
    axset_wind[1].set_xlabel(r"$\bar{w'w'}$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[1].set_xlim(0, 0.84)
    axset_wind[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    axset_wind[1].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)
    
    axset_wind[2].set_yticks(axset_wind[0].get_yticks())
    axset_wind[2].set_yticklabels(['' for t in axset_wind[1].get_yticks()])
    axset_wind[2].set_ylim(-100, 3300)
    axset_wind[2].set_xlabel(r"TKE (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[2].set_xlim(0, 1.3)
    axset_wind[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
    axset_wind[2].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

    axset_wind[3].set_yticks(axset_wind[0].get_yticks())
    axset_wind[3].set_ylim(-100, 3300)
    axset_wind[3].set_xlabel(r"$\bar{w'w'w'}$", fontsize=12)
    axset_wind[3].set_xlim(-1.3, 3.1)
    axset_wind[3].set_xticks([-1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3])
    axset_wind[3].set_xticklabels(['-1.2', '', '0', '', '1.2', '', '2.4', ''], fontsize=9)

    axset_scalar[0].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[0].set_xlim(-25, 25) 
    axset_scalar[0].set_ylim(-100, 3300) 
    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)
    axset_scalar[0].set_xlabel(r"SHF (W m$^{-2}$)", fontsize=12)

    axset_scalar[1].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[1].set_xlim(-60, 490) 
    axset_scalar[1].set_ylim(-100, 3300) 
    axset_scalar[1].set_xlabel(r"LHF (W m$^{-2}$)", fontsize=12)

    axset_scalar[2].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[2].set_xlim(-0.0005, 0.002) 
    axset_scalar[2].set_ylim(-100, 3300) 
    axset_scalar[2].set_xlabel(r"F$_b$ ($m^2/s^3$)", fontsize=12)
    
    axset_scalar[3].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[3].set_ylim(-100, 3300) 
    axset_scalar[3].set_xlim(-170, 50) 
    axset_scalar[3].set_xlabel(r"$\delta D_{flux}$"+u'(\u2030)', fontsize=12)
    
    yticklabels = axset_wind[0].get_yticks().astype(int).astype(str)
    axset_wind[0].set_yticklabels(yticklabels, fontsize=9)
    axset_wind[0].set_ylabel('altitude (m)', fontsize=12)
    for ax in axset_wind[1:]:
        ax.set_yticklabels(['' for t in ax.get_yticks()])

    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)
    axset_scalar[0].set_yticklabels(yticklabels, fontsize=9)
    for ax in axset_scalar[1:]:
        ax.set_yticklabels(['' for t in ax.get_yticks()])

    
    fig_wind.savefig("./fig_wind+turb_profiles.png")
    fig_scalar.savefig("./fig_scalarflux_profiles.png")



if __name__=="__main__":
    #fig_LCLCTscaling()
    #fig_LCLTIBscaling()
    #fig_noscaling()
    fig_scatter()