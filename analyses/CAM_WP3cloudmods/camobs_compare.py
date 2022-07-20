# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 09:54:59 2022

@author: Dean
"""


# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



## I/O paths:
path_camdir = "./CAM_extracted/"
path_p3drpsnd_dir = "../WP3_cloudmodule_char/cldmod_datafiles/"
path_p3prfflux_dir = "../WP3_fluxes/mean_profiles/"
path_p3prfthermo_dir = "../WP3_cloudmodule_char/mean_profiles/"



ncld_g1 = [1, 5, 4]
ncld_g2 = [8, 7, 9, 11, 10, 6]
ncld_g3 = [15, 3, 2, 13, 16, 14, 12]  

ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()
ncld_groups = [ncld_g1, ncld_g2, ncld_g3]



def plotthermo_singlecldmod(ncld, color, axset, varkeys): # input str with zfill=2.

    """
    # Load dropsondedata:
    fnames_drpsnds = [f for f in os.listdir(path_p3drpsnd_dir) if "_dropsondes_" in f]
    fname_drpsnds = [f for f in fnames_drpsnds if "_ncld%s" % ncld in f]
    fname_drpsnds = fname_drpsnds[0]
    #fname_drpsnds = "p3cld_dropsondes_20200117_ncld01.nc"
    drpsnds = xr.load_dataset(os.path.join(path_p3drpsnd_dir, fname_drpsnds))
    """
    # Load CAM data:
    fnames_cam = os.listdir(path_camdir)
    fname_cam = [f for f in fnames_cam if "_cld%s" % ncld in f]
    fname_cam = fname_cam[0]
    #fname_cam = "cam_cldextract_20200117_cld01.nc"
    cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
    
    
    # Plot CAM vertical profiles:
    #fig, axset = plt.subplots(1, 2, figsize=(10, 4))
    for ax, vk in zip(axset, varkeys):
        ax.plot(cam[vk], cam['P'], color=color)
        ax.set_ylim(100100, 70000)
    
    
    # Plot dropsonde vertical profiles:
    def plot_alldrpsnds(data, varkey, ax):
        for snd in data['sounding']:
            data_snd = data.sel(sounding=snd)
            ax.plot(data_snd[varkey], data_snd['p'], color='black')
    """        
    plot_alldrpsnds(drpsnds, 'theta', axset[0])
    plot_alldrpsnds(drpsnds, 'q', axset[1])
    """
    
    #axset[0].set_ylabel('altitude (m)', fontsize=12)
    #axset[0].set_xlabel('theta (K)', fontsize=12)
    #axset[1].set_xlabel('q (kg/kg)', fontsize=12)



def plotflux_singlecldmod(ncld, color, axset, varkeys): # input str with zfill=2.

    # Load CAM data:
    fnames_cam = os.listdir(path_camdir)
    fname_cam = [f for f in fnames_cam if "_cld%s" % ncld in f]
    fname_cam = fname_cam[0]
    #fname_cam = "cam_cldextract_20200117_cld01.nc"
    cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
    
    
    # Plot CAM vertical profiles:
    #fig, axset = plt.subplots(1, 2, figsize=(10, 4))
    for ax, vk in zip(axset, varkeys):
        ax.plot(cam[vk][1:], cam['P'], color=color)
        ax.set_ylim(100100, 70000)
        
        

def max_prf(camdata, varkey):
    """
    Return vertical profile at CAM grid point (gp) with the maximum value at 
    any one of its levels across the passed domain.
    """
    varmax = camdata[varkey].max()
    camcoords_gpmax = camdata.where(camdata[varkey]==varmax, drop=True).coords
    return camdata.sel(lon=camcoords_gpmax['lon'], lat=camcoords_gpmax['lat'])
    
    

def camstats(ncld_list, varkeys):
    """
    Returns follow for CAM output at timestamps for each of the input cloud 
    module numbers:
        (1) mean and standard deviation of CAM data over the ATOMIC study 
            region and the set of cloud modules.
        (2) mean and standard deviation of gridpoint profiles with maximum 
            vertical velocity variance for the set of cloud modules.
    """
    cam_list = []
    # Load and append cam data for each cloud module:
    for n in ncld_list:
        n_str = str(n).zfill(2)
        fnames_cam = [f for f in os.listdir(path_camdir) if "_ATOMICextract" in f]
        fname_cam = [f for f in fnames_cam if "_cld%s" % n_str in f]
        fname_cam = fname_cam[0]
        cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
        cam_list.append(cam[varkeys + ['P', 'P_ilevs', 'WP2_CLUBB']])
        
    # All data over P-3 region for all cloud modules:
    cam_all = xr.concat(cam_list, dim='ncld')
        
    # Gridpoint profiles with maximum vertical velocity variance for each 
    # cloud module:
    cam_maxprfs = [max_prf(camdata, 'WP2_CLUBB') for camdata in cam_list]
    cam_maxprfs = xr.concat(cam_maxprfs, dim='ncld')
    

    return (
        (cam_all.mean(dim=['ncld','lat','lon']), cam_all.std(dim=['ncld','lat','lon'])), 
        (cam_maxprfs.mean(dim=['ncld','lat','lon']), cam_maxprfs.std(dim=['ncld','lat','lon']))
        )



def camstats_v2(ncld_list, varkeys):
    """
    Returns follow for CAM output at timestamps for each of the input cloud 
    module numbers:
        (1) mean and standard deviation of CAM data over the ATOMIC study 
            region and the set of cloud modules.
        (2) mean and standard deviation of gridpoint profiles with maximum 
            vertical velocity variance for the set of cloud modules.
    """
    cam_list = []
    # Load and append cam data for each cloud module:
    for n in ncld_list:
        n_str = str(n).zfill(2)
        fnames_cam = [f for f in os.listdir(path_camdir) if "_ATOMICextract" in f]
        fname_cam = [f for f in fnames_cam if "_cld%s" % n_str in f]
        fname_cam = fname_cam[0]
        cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
        cam_list.append(cam[varkeys + ['P', 'P_ilevs', 'WP2_CLUBB']])
        
    # All data over P-3 region for all cloud modules:
    cam_all = xr.concat(cam_list, dim='ncld')
    
    # Top/bottom quantiles of column integrated vertical velocity variance:
    cam_all['WP2_CLUBB_column'] = cam_all['WP2_CLUBB'].sum(dim='ilev')
    q_15p = np.quantile(cam_all['WP2_CLUBB_column'], 0.15)
    q_85p = np.quantile(cam_all['WP2_CLUBB_column'], 0.85)
    cam_wp15p = cam_all.where(cam_all['WP2_CLUBB_column'] < q_15p, drop=True)
    cam_wp85p = cam_all.where(cam_all['WP2_CLUBB_column'] > q_85p, drop=True)
    cam_wp15p = cam_wp15p.mean(dim=['lon','lat','ncld'])
    cam_wp85p = cam_wp85p.mean(dim=['lon','lat','ncld'])

    cam_mean = cam_all.mean(dim=['ncld','lat','lon'])

    return (cam_mean, cam_wp15p, cam_wp85p)



def collect_prfs(ncld_list, varkeys):
    """
    Return mean profiles and standard deviation on the mean.
    """
    prfs_list = []
    # Load and append cam data for each cloud module:
    for n in ncld_list:
        n_str = str(n).zfill(2)
        fnames_cam = [f for f in os.listdir(path_camdir) if "_cldextract" in f]
        fname_cam = [f for f in fnames_cam if "_cld%s" % n_str in f]
        fname_cam = fname_cam[0]
        cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
        prfs_list.append(cam[varkeys + ['P']])
    
    # Return mean, std:
    prfs_all = xr.concat(prfs_list, dim='profile')
    meanprfs = prfs_all.mean(dim='profile')
    stdprfs = prfs_all.std(dim='profile')/len(ncld_list)**0.5
    return meanprfs, stdprfs




def plot_prfs(ncld_list, color, axset, varkeys, alpha, 
              levs='mid', prftype='region'):
    """
    For CAM quantities on interface-levels.
    """
    #meanprfs, stdprfs = collect_prfs(ncld_list, varkeys)
    regionstats, maxstats = camstats(ncld_list, varkeys)
    if prftype=='region':
        mean, std = regionstats
    elif prftype=='max':
        mean, std = maxstats
    
    if levs=='mid':
        z = (mean['P'][-1]-mean['P'])/10 # Rough altitude estimation, needs later revision.
    elif levs=='int':
        z = (mean['P_ilevs'][-1]-mean['P_ilevs'])/10 # Rough altitude estimation, needs later revision.
    
    for ax, vk in zip(axset, varkeys):
        #ax.plot(meanprfs[vk][1:], meanprfs['P']/100, color=color)
        ax.plot(mean[vk], z, color=color, linewidth=2)
        ax.fill_betweenx(
            #meanprfs['P']/100, 
            z, 
            mean[vk] - std[vk], 
            mean[vk] + std[vk], 
            color=color, alpha=alpha
            )



def plot_p3prfs(path_p3prfdir, varkey, ax, 
                colors=['red', 'grey', 'blue'], 
                binning=False):
    """
    Plot P-3 mean profiles of a quantity on the specified axes. Plot  one 
    variable, 3 profiles (one for each cloud group).
    """
    for n, c in zip([1, 2, 3], colors):
        fname = "meanprf_%s_WP3.csv" % varkey
        p3data = pd.read_csv(path_p3prfdir + fname)
        if binning:
            p3data = p3data.groupby(150*np.round(p3data['alt']/150)).mean()
        ax.plot(
            p3data['cg%i' % n], p3data['alt'], 
            color=c, linestyle='dashed', linewidth=1.
            )    



## Thermodynamic quantity profiles
##_____________________________________________________________________________
varkeys_thermo = ['theta', 'H216OV', 'RH', 'dD']
#fig, axset = plt.subplots(1, 4, figsize=(10, 4))
fig_thermo = plt.figure(figsize=(6.5, 2.5))
axset = [
    fig_thermo.add_axes([0.1, 0.2, 0.2, 0.75]),
    fig_thermo.add_axes([0.325, 0.2, 0.2, 0.75]),
    fig_thermo.add_axes([0.55, 0.2, 0.2, 0.75]),
    fig_thermo.add_axes([0.775, 0.2, 0.2, 0.75]),
    ]

# Plot CAM:
#for n in ncld_g1: plotthermo_singlecldmod(str(n).zfill(2), 'red', axset, varkeys_thermo)
#for n in ncld_g2: plotthermo_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys_thermo)
#for n in ncld_g3: plotthermo_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys_thermo)
plot_prfs(ncld_g3, 'blue', axset, varkeys_thermo, 0.1, levs='mid')
plot_prfs(ncld_g2, 'grey', axset, varkeys_thermo, 0.2, levs='mid')
plot_prfs(ncld_g1, 'red', axset, varkeys_thermo, 0.2, levs='mid')

# Plot P-3:
p3_prfvarkeys = ["theta", "q", "RH", "dD"]
for ax, vk in zip(axset, p3_prfvarkeys):
    plot_p3prfs(path_p3prfthermo_dir, vk, ax, binning=True)    
    

axset[0].set_xlim(294, 315)
axset[0].set_xticks(np.arange(295, 315, 5))
axset[0].set_xticklabels(axset[0].get_xticks().astype(str), fontsize=9)
axset[0].set_xlabel(r'$\theta$ (K)', fontsize=12)

axset[1].set_xlim(-0.0006, 0.02)
axset[1].set_xticks(np.arange(0, 0.02, 0.005))
axset[1].set_xticklabels(
    [(1000*t).astype(int).astype(str) for t in axset[1].get_xticks()],
    fontsize=9
    ) # factor of 1000 to convert to g/kg.
axset[1].set_xlabel(r'$q$ (g/kg)', fontsize=12)

axset[2].set_xlim(-0.02, 1.05)
axset[2].set_xticks(np.arange(0, 1.01, 0.25))
axset[2].set_xticklabels(
    [(100*t).astype(int).astype(str) for t in axset[2].get_xticks()],
    fontsize=9
    ) # factor of 100 to convert to %.
axset[2].set_xlabel(r'$RH$ (%)', fontsize=12)

axset[3].set_xlim(-285, -40)
axset[3].set_xticks(np.arange(-250, -49, 50))
axset[3].set_xticklabels(axset[3].get_xticks().astype(str), fontsize=9)
axset[3].set_xlabel(r'$\delta D$ '+u'(\u2030)', fontsize=12)

for ax in axset:
    ax.set_ylim(0, 3200)
    ax.set_yticks(np.arange(0, 3200, 500))
axset[0].set_yticklabels(ax.get_yticks().astype(str), fontsize=9)
axset[0].set_ylabel('altitude (m)', fontsize=12)
for ax in axset[1:]: 
    ax.set_yticklabels(['' for t in ax.get_yticks()], fontsize=9)
    

legend_lines = [Line2D([0], [0], color='grey', lw=1.5, linestyle='dashed'),
                Line2D([0], [0], color='grey', lw=2, linestyle='solid')]
axset[0].legend(legend_lines, ['P-3', 'CAM'], loc='lower right', 
                handletextpad=0.3, frameon=False)    
    
    
fig_thermo.savefig("./fig_thermoprofiles_CAM.png")
##_____________________________________________________________________________
## Thermodynamic quantity profiles



def turbulenceflux_figure(prftype='region'):
    """
    """
    fig_turb = plt.figure(figsize=(6.5, 5))
    axset_toprow = [
        fig_turb.add_axes([0.1, 0.6, 0.2, 0.375]),
        fig_turb.add_axes([0.325, 0.6, 0.2, 0.375]),
        fig_turb.add_axes([0.55, 0.6, 0.2, 0.375]),
        fig_turb.add_axes([0.775, 0.6, 0.2, 0.375]),
        ]
    axset_bottomrow = [
        fig_turb.add_axes([0.2, 0.1, 0.2, 0.375]),
        fig_turb.add_axes([0.45, 0.1, 0.2, 0.375]),
        fig_turb.add_axes([0.7, 0.1, 0.2, 0.375]),
        ]


    # Plot CAM top row:
    camvarkeys_toprow = ['TKE_h', 'WP2_CLUBB', 'TKE', 'anisotropy_ratio']
    plot_prfs(
        ncld_g3, 'blue', axset_toprow, camvarkeys_toprow, 0.1, 
        levs='int', prftype=prftype
        )
    plot_prfs(
        ncld_g2, 'grey', axset_toprow, camvarkeys_toprow, 0.2, 
        levs='int', prftype=prftype
        )
    plot_prfs(
        ncld_g1, 'red', axset_toprow, camvarkeys_toprow, 0.2, 
        levs='int', prftype=prftype
        )
    
    
    # Plot P-3 top row:
    p3varkeys_toprow = ["TKE_h", "wp2_bar", "TKE", "anisotropy_ratio"]
    #for ax, vk in zip(axset[0, :], p3varkeys_toprow):
    for ax, vk in zip(axset_toprow, p3varkeys_toprow):
        plot_p3prfs(path_p3prfflux_dir, vk, ax)
        
        
    # Plot CAM bottom row:
    camvarkeys_bottomrow = ['WPTHLP_CLUBB', 'WPRTP_CLUBB', 'Fb_m2s3']
    plot_prfs(
        ncld_g3, 'blue', axset_bottomrow[0:2], camvarkeys_bottomrow[0:2], 0.1, 
        levs='int', prftype=prftype
        )
    plot_prfs(
        ncld_g2, 'grey', axset_bottomrow[0:2], camvarkeys_bottomrow[0:2], 0.2, 
        levs='int', prftype=prftype
        )
    plot_prfs(
        ncld_g1, 'red', axset_bottomrow[0:2], camvarkeys_bottomrow[0:2], 0.2, 
        levs='int', prftype=prftype
        )
    plot_prfs(
        ncld_g3, 'blue', [axset_bottomrow[2]], [camvarkeys_bottomrow[2]], 0.1, 
        levs='mid', prftype=prftype
        )
    plot_prfs(
        ncld_g2, 'grey', [axset_bottomrow[2]], [camvarkeys_bottomrow[2]], 0.2, 
        levs='mid', prftype=prftype
        )
    plot_prfs(
        ncld_g1, 'red', [axset_bottomrow[2]], [camvarkeys_bottomrow[2]], 0.2, 
        levs='mid', prftype=prftype
        )
    
    
    # Plot P-3 bottom row:
    p3varkeys_bottomrow = ["flux_sh", "flux_lh", "flux_b"]
    #for ax, vk in zip(axset[1, 0:3], p3varkeys_bottomrow):
    for ax, vk in zip(axset_bottomrow, p3varkeys_bottomrow):
        plot_p3prfs(path_p3prfflux_dir, vk, ax)


    for ax in axset_toprow + axset_bottomrow:
        ax.set_ylim(0, 3200)
        ax.set_yticks(np.arange(0, 3200, 500))
    axset_toprow[0].set_yticklabels(axset_toprow[0].get_yticks().astype(str), fontsize=9)
    axset_toprow[0].set_ylabel('altitude (m)', fontsize=12)
    axset_bottomrow[0].set_yticklabels(axset_bottomrow[0].get_yticks().astype(str), fontsize=9)
    axset_bottomrow[0].set_ylabel('altitude (m)', fontsize=12)
    for ax in axset_toprow[1:] + axset_bottomrow[1:]: 
        ax.set_yticklabels(['' for t in ax.get_yticks()], fontsize=9)


    axset_toprow[0].set_xlabel(r"TKE$_h$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_toprow[0].set_xlim(0, 1.2)
    axset_toprow[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
    axset_toprow[0].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)
 
    axset_toprow[1].set_xlabel(r"$\bar{w'w'}$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_toprow[1].set_xlim(0, 0.84)
    axset_toprow[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    axset_toprow[1].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)
    
    axset_toprow[2].set_xlabel(r"TKE (m$^2$ s$^{-2}$)", fontsize=12)
    axset_toprow[2].set_xlim(0, 1.3)
    axset_toprow[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
    axset_toprow[2].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

    axset_toprow[3].set_xlabel("anisotropy \nratio", fontsize=12)
    axset_toprow[3].set_xlim(0, 2.5)
    axset_toprow[3].set_xticks([0, 0.5, 1, 1.5, 2])
    axset_toprow[3].set_xticklabels(['0', '0.5', '1.0', '1.5', '2.0'], fontsize=9)


    axset_bottomrow[0].set_ylabel('altitude (m)', fontsize=12)
    axset_bottomrow[0].set_xlabel(r"SHF (W m$^{-2}$)", fontsize=12)

    axset_bottomrow[1].set_xlim(-60, 490) 
    axset_bottomrow[1].set_ylim(-100, 3300) 
    axset_bottomrow[1].set_xlabel(r"LHF (W m$^{-2}$)", fontsize=12)

    axset_bottomrow[2].set_xlim(-0.00025, 0.0015) 
    axset_bottomrow[2].set_xticks([0, 0.0005, 0.001, 0.0015]) 
    axset_bottomrow[2].set_xticklabels(['0', '0.5', '1.0', '1.5']) 
    axset_bottomrow[2].set_ylim(-100, 3300) 
    axset_bottomrow[2].set_xlabel(r"F$_b$ ($m^2/s^3$ $10^{-3}$)", fontsize=12)
    
    legend_lines = [Line2D([0], [0], color='grey', lw=1.5, linestyle='dashed'),
                    Line2D([0], [0], color='grey', lw=2, linestyle='solid')]
    axset_toprow[0].legend(legend_lines, ['P-3', 'CAM'], loc='upper right', 
                    handletextpad=0.3, frameon=False)      


    return fig_turb


fig_turbregion = turbulenceflux_figure(prftype='region')
fig_turbregion.savefig("./fig_turbulence+flux_profiles_CAM_region.png")

fig_turbmax = turbulenceflux_figure(prftype='max')
fig_turbmax.savefig("./fig_turbulence+flux_profiles_CAM_max.png")


"""
## Mean turbulence/flux profiles for ATOMIC region
##_____________________________________________________________________________
varkeys_turb = ['TKE_h', 'WP2_CLUBB', 'TKE']
fig, axset_turb = plt.subplots(1, 4, figsize=(10, 4))

# Plot CAM:
plot_prfs(ncld_g3, 'blue', axset_turb[0:3], varkeys_turb, 0.2, levs='int')
plot_prfs(ncld_g2, 'grey', axset_turb[0:3], varkeys_turb, 0.4, levs='int')
plot_prfs(ncld_g1, 'red', axset_turb[0:3], varkeys_turb, 0.4, levs='int')

# Plot P-3:
#varkeys = ["wp2_bar", "Sw"]
varkeys = ["TKE_h", "wp2_bar", "TKE"]
#for ax, vk in zip(axset_turb[[1,3]], varkeys):
for ax, vk in zip(axset_turb[0:3], varkeys):
    plot_p3prfs(path_p3prfflux_dir, vk, ax)    


#axset_turb[0].set_yticks(axset_turb[0].get_yticks())
#axset_turb[0].set_yticklabels(['' for t in axset_turb[1].get_yticks()])
#axset_turb[0].set_ylim(-100, 3300)
#axset_turb[0].set_xlabel(r"TKE$_h$ (m$^2$ s$^{-2}$)", fontsize=12)
#axset_turb[0].set_xlim(0, 1.2)
#axset_turb[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#axset_turb[0].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

#axset_turb[1].set_yticks(axset_turb[0].get_yticks())
#axset_turb[1].set_ylim(-100, 3300) 
#axset_turb[1].set_xlabel(r"$\bar{w'w'}$ (m$^2$ s$^{-2}$)", fontsize=12)
#axset_turb[1].set_xlim(0, 0.84)
#axset_turb[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
#axset_turb[1].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)

#axset_turb[2].set_yticks(axset_turb[0].get_yticks())
#axset_turb[2].set_yticklabels(['' for t in axset_turb[1].get_yticks()])
#axset_turb[2].set_ylim(-100, 3300)
#axset_turb[2].set_xlabel(r"TKE (m$^2$ s$^{-2}$)", fontsize=12)
#axset_turb[2].set_xlim(0, 1.3)
#axset_turb[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#axset_turb[2].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

#axset_turb[3].set_yticks(axset_turb[0].get_yticks())
#axset_turb[3].set_ylim(-100, 3300)
##axset_turb[3].set_xlabel(r"$\bar{w'w'w'}$", fontsize=12)
#axset_turb[3].set_xlabel(r"$S_w$", fontsize=12)
#axset_turb[3].set_xlim(-1.3, 3.1)
#axset_turb[3].set_xticks([-1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3])
#axset_turb[3].set_xticklabels(['-1.2', '', '0', '', '1.2', '', '2.4', ''], fontsize=9)

for ax in axset_turb:
    ax.set_ylim(0, 3200)
    ax.set_yticks(np.arange(0, 3200, 500))
axset_turb[0].set_yticklabels(ax.get_yticks().astype(str), fontsize=9)
axset_turb[0].set_ylabel('altitude (m)', fontsize=12)
for ax in axset_turb[1:]: 
    ax.set_yticklabels(['' for t in ax.get_yticks()], fontsize=9)
    
fig.savefig("./fig_wind+turb_profiles_CAM.png")
##_____________________________________________________________________________
## Mean turbulence/flux profiles for ATOMIC region
"""
    
""" 
## Flux profiles
##_____________________________________________________________________________
varkeys_flux = ['WPTHLP_CLUBB', 'WPRTP_CLUBB', 'WPTHVP_CLUBB']
fig, axset_flux = plt.subplots(1, 4, figsize=(10, 4))

# Plot CAM:
#for n in ncld_g1: plotflux_singlecldmod(str(n).zfill(2), 'red', axset, varkeys_flux)
#for n in ncld_g2: plotflux_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys_flux)
#for n in ncld_g3: plotflux_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys_flux)
plot_prfs(ncld_g3, 'blue', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.2, levs='int')
plot_prfs(ncld_g2, 'grey', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.4, levs='int')
plot_prfs(ncld_g1, 'red', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.4, levs='int')
plot_prfs(ncld_g3, 'blue', [axset_flux[2]], [varkeys_flux[2]], 0.2, levs='mid')
plot_prfs(ncld_g2, 'grey', [axset_flux[2]], [varkeys_flux[2]], 0.4, levs='mid')
plot_prfs(ncld_g1, 'red', [axset_flux[2]], [varkeys_flux[2]], 0.4, levs='mid')

# Plot P-3:
varkeys = ["flux_sh", "flux_lh"]
for ax, vk in zip(axset_flux[[0,1]], varkeys):
    plot_p3prfs(path_p3prfflux_dir, vk, ax)      
    
for ax in axset_flux:
    ax.set_ylim(0, 3200)
    ax.set_yticks(np.arange(0, 3200, 500))
axset_flux[0].set_yticklabels(ax.get_yticks().astype(str), fontsize=9)
axset_flux[0].set_ylabel('altitude (m)', fontsize=12)
for ax in axset_flux[1:]: 
    ax.set_yticklabels(['' for t in ax.get_yticks()], fontsize=9)
    
axset_flux[0].set_yticks(axset_flux[0].get_yticks())
axset_flux[0].set_xlim(-25, 25) 
axset_flux[0].set_ylim(-100, 3300) 
axset_flux[0].set_ylabel('altitude (m)', fontsize=12)
axset_flux[0].set_xlabel(r"SHF (W m$^{-2}$)", fontsize=12)

axset_flux[1].set_yticks(axset_flux[0].get_yticks())
axset_flux[1].set_xlim(-60, 490) 
axset_flux[1].set_ylim(-100, 3300) 
axset_flux[1].set_xlabel(r"LHF (W m$^{-2}$)", fontsize=12)

#axset_flux[2].set_yticks(axset_flux[0].get_yticks())
#axset_flux[2].set_xlim(-0.0005, 0.002) 
#axset_flux[2].set_ylim(-100, 3300) 
axset_flux[2].set_xlabel(r"F$_b$ (W/$m^2$)", fontsize=12)

#axset_flux[3].set_yticks(axset_flux[0].get_yticks())
#axset_flux[3].set_ylim(-100, 3300) 
#axset_flux[3].set_xlim(-170, 50) 
#axset_flux[3].set_xlabel(r"$\delta D_{flux}$"+u'(\u2030)', fontsize=12)

fig.savefig("./fig_scalarflux_profiles.png")
##_____________________________________________________________________________
## Flux profiles
"""

fig_turb = plt.figure(figsize=(6.5, 5))
axset_toprow = [
    fig_turb.add_axes([0.1, 0.6, 0.2, 0.375]),
    fig_turb.add_axes([0.325, 0.6, 0.2, 0.375]),
    fig_turb.add_axes([0.55, 0.6, 0.2, 0.375]),
    fig_turb.add_axes([0.775, 0.6, 0.2, 0.375]),
    ]
axset_bottomrow = [
    fig_turb.add_axes([0.2, 0.1, 0.2, 0.375]),
    fig_turb.add_axes([0.45, 0.1, 0.2, 0.375]),
    fig_turb.add_axes([0.7, 0.1, 0.2, 0.375]),
    ]

varkeys = ['TKE_h', 'WP2_CLUBB', 'TKE', 'anisotropy_ratio',
           'WPTHLP_CLUBB', 'WPRTP_CLUBB', 'WPTHVP_CLUBB']
levtype = ['int', 'int', 'int', 'int', 'int', 'int', 'mid']
ncld_list = list(range(1,17))
results = camstats_v2(ncld_list, varkeys)


for ds in results:
    ds['z_mid'] = (ds['P'][-1]-ds['P'])/10 # Rough altitude estimation, needs later revision.
    ds['z_int'] = (ds['P_ilevs'][-1]-ds['P_ilevs'])/10 # Rough altitude estimation, needs later revision.
    

for vk, lvt, ax in zip(varkeys, levtype, axset_toprow + axset_bottomrow):
    
    ax.plot(results[0][vk], results[0]['z_'+lvt], color='black')
    ax.plot(results[1][vk], results[1]['z_'+lvt], color='grey')
    ax.plot(results[2][vk], results[2]['z_'+lvt], color='red')

    ax.set_ylim(0, 3500)



