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


dir_flux = "./fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]


dir_prf = "../../data/WP3/cloud_modules/"
fnames_prf = [f for f in os.listdir(dir_prf) if f.endswith('.nc')]


ncld = np.arange(1, 17, 1).astype(str)


def plot_fluxprf(ncld):
    """
    """
    fname_flux = [f for f in fnames_flux if "_cld%s_" % n in f]
    if len(fname_flux) == 0: return
    fname_flux = fname_flux[0]
    fname_prf = [f for f in fnames_prf if "_%s.nc" % n.zfill(2) in f]
    fname_prf = fname_prf[0]
    
    
    fluxes = xr.load_dataset(dir_flux + fname_flux)
    prf = xr.load_dataset(dir_prf + fname_prf)
    prfmean = prf.groupby(np.round(prf['alt']/50)*50).mean()

    
    fig, axes = plt.subplots(1, 3, figsize=(8, 3))
    for k, ax in zip(["w'w'_bar", "flux_sh", "flux_lh"], axes):
        ax.scatter(fluxes[k], fluxes['alt'])
    
    
    axes[2].twiny().plot(prfmean['mr'], prfmean['alt'])
    theta = thermo.ptemp(prfmean['Ta'], prfmean['press'])
    axes[0].twiny().plot(theta, prfmean['alt'])
    
    
    axes[0].set_ylabel(fname_flux[4:18], fontsize=12)
    for ax in axes:
        ax.set_ylim(0, 3250)
    axes[0].set_xlim(0, np.max(fluxes["w'w'_bar"]+0.02))
    axes[2].set_xlim(0, np.max(fluxes["flux_lh"]+20))


#for n in ncld:
#    plot_fluxprf(ncld)    
    

# ## Profiles with scaled altitude

keyalts = pd.read_csv("./cldmod_keyaltitudes.csv")
g1 = (keyalts['ncld']>=4) & (keyalts['ncld']<=11)
keyalts_g1 = keyalts.loc[g1]
keyalts_g2 = keyalts.loc[~g1 & (keyalts['z_mltop']!=(-1))]

"""
fig, axes = plt.subplots(1, 4, figsize=(10, 5))

for n in keyalts_g1['ncld']:
    
    n = str(n)
    keyalt_info = keyalts_g1.loc[keyalts_g1['ncld']==int(n)]
    
    fname_flux = [f for f in fnames_flux if "_cld%s_" % n in f]
    if len(fname_flux) == 0: continue
    fname_flux = fname_flux[0]
    fname_prf = [f for f in fnames_prf if "_%s.nc" % n.zfill(2) in f]
    fname_prf = fname_prf[0]
    
    
    fluxes = xr.load_dataset(dir_flux + fname_flux)
    prf = xr.load_dataset(dir_prf + fname_prf)
    prfmean = prf.groupby(np.round(prf['alt']/50)*50).mean()

    zflux_scaled = rangescaler.piecewise_linscale(
        fluxes['alt'].values, (0, keyalt_info['z_mltop'].values, keyalt_info['z_tradeinversion'].values), 
        (0,1,2)
        )
    
    for k, ax in zip(["w'w'_bar", "flux_sh", "flux_lh", "dD_flux"], axes):
        ax.scatter(fluxes[k], zflux_scaled, s=3, c='grey')
        ax.plot(fluxes[k], zflux_scaled, c='grey', alpha=0.3)
        
    
axes[0].set_xlim(0, 0.65)
axes[2].set_xlim(-20, 650)
axes[3].set_xlim(-120, -40)
yticks = [0,1,2,3,4]
for ax in axes:
    ax.set_yticks(yticks)
axes[0].set_yticklabels([
    "0", r"$z_{ML}$", r"$z_{IB} = z_{ML} + \Delta z_{CL}$", 
    r"$z_{ML} + 2\Delta z_{CL}$", r"$z_{ML} + 3\Delta z_{CL}$"
    ])
for ax in axes[1:]: ax.set_yticklabels(["" for t in yticks])
"""


def get_fluxprofiles(ncld_list, keyalts_table, dir_flux):
    """
    keyalts_table: pandas.DataFrame.
        Cannot have NAN's for the mixed layer top or trade inversion bottom.
    """

    keyalts_table = keyalts_table.set_index('ncld')


    # All flux filenames:
    fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]


    # Collect flux profiles into a dictionary of pandas.DataFrame's.
    # Each DataFrame will contain all profiles for a specific var.
    fluxvarkeys = ["w'w'_bar", "flux_sh", "flux_lh", "dD_flux"]
    fluxprfs_dict = {}
    for k in fluxvarkeys: fluxprfs_dict[k] = pd.DataFrame({})
    
    
    for n in ncld_list:
        
        # Key altitude info for this module:
        keyalts = keyalts_table.loc[n]
        
        # Load flux data:
        fname_flux = [f for f in fnames_flux if "_cld%i_" % n in f]
        if len(fname_flux) == 0: continue
        fname_flux = fname_flux[0]        
        fluxes = xr.load_dataset(dir_flux + fname_flux)
    
        # Scale altitude.
        # 0=sea level, 1=mixed layer top, 2=trade inversion bottom.
        alt_scaled = rangescaler.piecewise_linscale(
            fluxes['alt'].values, 
            (0, keyalts['z_mltop'], keyalts['z_tradeinversion']), 
            (0,1,2)
            )
        
        # Merge flux profiles of current module into resp. DataFrame's:
        for k in fluxvarkeys:
            
            fluxvar_pd = pd.DataFrame(  # easier to work in pandas.
                {k+"_cld%i"%n: fluxes[k].values}, 
                index=alt_scaled
                )
                
            fluxprfs_dict[k] = pd.merge(
                fluxprfs_dict[k], fluxvar_pd, 
                how='outer', left_index=True, right_index=True
                )   
            
            
    return fluxprfs_dict


keyalts_table = pd.read_csv("./cldmod_keyaltitudes.csv")
ncld_list = np.arange(4, 12, 1)
flux_prfs = get_fluxprofiles(ncld_list, keyalts_table, dir_flux)

fig, axes = plt.subplots(1, 4, figsize=(10, 5))

for varkey, ax in zip(flux_prfs.keys(), axes):

    for k in flux_prfs[varkey].columns:
        colplot = flux_prfs[varkey][k].dropna()
        ax.plot(colplot, colplot.index, c='grey', alpha=0.5)
        ax.scatter(colplot, colplot.index, c='grey', s=10, alpha=0.5)  
    
    #altgrouped = np.round(flux_prfs[varkey].index/0.25)*0.25
    altgrouped = np.round(flux_prfs[varkey].index/0.5)*0.5
    #altgrouped = flux_prfs[varkey].index.values//0.25
    #altgrouped[altgrouped % 2 == 0] += 1
    #altgrouped = altgrouped*0.25
    fluxvar_grouped = flux_prfs[varkey].groupby(altgrouped, axis=0, as_index=True)
    
    meanprf = fluxvar_grouped.mean().mean(axis=1)
    ax.plot(meanprf.values, meanprf.index, 'b-', linewidth=5)
    ax.scatter(meanprf.values, meanprf.index, c='b', s=10)
      
    
axes[0].set_xlim(0, 0.65)
axes[2].set_xlim(-20, 650)
axes[3].set_xlim(-120, -40)
yticks = [0,1,2,3,4]
for ax in axes:
    ax.set_yticks(yticks)
axes[0].set_yticklabels([
    "0", r"$z_{ML}$", r"$z_{IB} = z_{ML} + \Delta z_{CL}$", 
    r"$z_{ML} + 2\Delta z_{CL}$", r"$z_{ML} + 3\Delta z_{CL}$"
    ])
for ax in axes[1:]: ax.set_yticklabels(["" for t in yticks])

axes[0].set_xlabel(r"w'w' (m/s)$^2$", fontsize=14)
axes[1].set_xlabel("SH flux (W/m$^2$)", fontsize=14)
axes[2].set_xlabel("LH flux (W/m$^2$)", fontsize=14)
axes[3].set_xlabel(r"$\delta D_{flux}$ (permil)", fontsize=14)



