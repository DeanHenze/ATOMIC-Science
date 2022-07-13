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
        
        
        
def camstats(ncld_list, varkeys):
    """
    Return mean profiles and standard deviation on the mean.
    """
    prfs_list = []
    # Load and append cam data for each cloud module:
    for n in ncld_list:
        n_str = str(n).zfill(2)
        fnames_cam = os.listdir(path_camdir)
        fname_cam = [f for f in fnames_cam if "_cld%s" % n_str in f]
        fname_cam = fname_cam[0]
        cam = xr.load_dataset(os.path.join(path_camdir, fname_cam))
        prfs_list.append(cam[varkeys + ['P']])
    
    # Return mean, std:
    prfs_all = xr.concat(prfs_list, dim='profile')
    meanprfs = prfs_all.mean(dim='profile')
    stdprfs = prfs_all.std(dim='profile')/len(ncld_list)**0.5
    return meanprfs, stdprfs



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



def plot_meanprf_mlevs(ncld_list, color, axset, varkeys, alpha):
    """
    For CAM quantities on mid-levels.
    """
    meanprfs, stdprfs = collect_prfs(ncld_list, varkeys)
    
    z = (meanprfs['P'][-1]-meanprfs['P'])/10 # Rough altitude estimation, needs later revision.

    for ax, vk in zip(axset, varkeys):
        #ax.plot(meanprfs[vk], meanprfs['P']/100, color=color)
        ax.plot(meanprfs[vk], z, color=color, linewidth=3)
        ax.fill_betweenx(
            #meanprfs['P']/100, 
            z, 
            meanprfs[vk] - stdprfs[vk], meanprfs[vk] + stdprfs[vk], 
            color=color, alpha=alpha
            )

        
        
def plot_meanprf_ilevs(ncld_list, color, axset, varkeys, alpha):
    """
    For CAM quantities on interface-levels.
    """
    meanprfs, stdprfs = collect_prfs(ncld_list, varkeys)
    
    z = (meanprfs['P'][-1]-meanprfs['P'])/10 # Rough altitude estimation, needs later revision.
    
    for ax, vk in zip(axset, varkeys):
        #ax.plot(meanprfs[vk][1:], meanprfs['P']/100, color=color)
        ax.plot(meanprfs[vk][1:], z, color=color, linewidth=3)
        ax.fill_betweenx(
            #meanprfs['P']/100, 
            z, 
            meanprfs[vk][1:] - stdprfs[vk][1:], 
            meanprfs[vk][1:] + stdprfs[vk][1:], 
            color=color, alpha=alpha
            )



def plot_p3prfs(path_p3prfdir, varkey, ax, colors=['red', 'grey', 'blue']):
    """
    Plot P-3 mean profiles of a quantity on the specified axes. Plot  one 
    variable, 3 profiles (one for each cloud group).
    """
    for n, c in zip([1, 2, 3], colors):
        fname = "meanprf_%s_WP3.csv" % varkey
        p3data = pd.read_csv(path_p3prfdir + fname)
        ax.plot(
            p3data['cg%i' % n], p3data['alt'], 
            color=c, linestyle='dashed', linewidth=1.
            )    



## Thermodynamic quantity profiles
##_____________________________________________________________________________
varkeys_thermo = ['theta', 'H216OV', 'RH', 'dD']
fig, axset = plt.subplots(1, 4, figsize=(10, 4))

# Plot CAM:
#for n in ncld_g1: plotthermo_singlecldmod(str(n).zfill(2), 'red', axset, varkeys_thermo)
#for n in ncld_g2: plotthermo_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys_thermo)
#for n in ncld_g3: plotthermo_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys_thermo)
plot_meanprf_mlevs(ncld_g3, 'blue', axset, varkeys_thermo, 0.2)
plot_meanprf_mlevs(ncld_g2, 'grey', axset, varkeys_thermo, 0.4)
plot_meanprf_mlevs(ncld_g1, 'red', axset, varkeys_thermo, 0.4)

# Plot P-3:
p3_prfvarkeys = ["theta", "q", "RH", "dD"]
for ax, vk in zip(axset, p3_prfvarkeys):
    plot_p3prfs(path_p3prfthermo_dir, vk, ax)    
    

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
    
fig.savefig("./fig_thermoprofiles_CAM.png")
##_____________________________________________________________________________
## Thermodynamic quantity profiles



## Turbulence profiles
##_____________________________________________________________________________
varkeys_turb = ['WP2_CLUBB', 'Sw']
fig, axset_turb = plt.subplots(1, 4, figsize=(10, 4))

# Plot CAM:
#for n in ncld_g1: plotflux_singlecldmod(str(n).zfill(2), 'red', axset, varkeys_turb)
#for n in ncld_g2: plotflux_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys_turb)
#for n in ncld_g3: plotflux_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys_turb)
plot_meanprf_ilevs(ncld_g3, 'blue', axset_turb[[1,3]], varkeys_turb, 0.2)
plot_meanprf_ilevs(ncld_g2, 'grey', axset_turb[[1,3]], varkeys_turb, 0.4)
plot_meanprf_ilevs(ncld_g1, 'red', axset_turb[[1,3]], varkeys_turb, 0.4)

# Plot P-3:
varkeys = ["wp2_bar", "Sw"]
for ax, vk in zip(axset_turb[[1,3]], varkeys):
    plot_p3prfs(path_p3prfflux_dir, vk, ax)    


#axset_turb[0].set_yticks(axset_turb[0].get_yticks())
#axset_turb[0].set_yticklabels(['' for t in axset_turb[1].get_yticks()])
#axset_turb[0].set_ylim(-100, 3300)
#axset_turb[0].set_xlabel(r"TKE$_h$ (m$^2$ s$^{-2}$)", fontsize=12)
#axset_turb[0].set_xlim(0, 1.2)
#axset_turb[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#axset_turb[0].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

axset_turb[1].set_yticks(axset_turb[0].get_yticks())
axset_turb[1].set_ylim(-100, 3300) 
axset_turb[1].set_xlabel(r"$\bar{w'w'}$ (m$^2$ s$^{-2}$)", fontsize=12)
axset_turb[1].set_xlim(0, 0.84)
axset_turb[1].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
axset_turb[1].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)

#axset_turb[2].set_yticks(axset_turb[0].get_yticks())
#axset_turb[2].set_yticklabels(['' for t in axset_turb[1].get_yticks()])
#axset_turb[2].set_ylim(-100, 3300)
#axset_turb[2].set_xlabel(r"TKE (m$^2$ s$^{-2}$)", fontsize=12)
#axset_turb[2].set_xlim(0, 1.3)
#axset_turb[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1., 1.2])
#axset_turb[2].set_xticklabels(['0', '', '0.4', '', '0.8', '', '1.2'], fontsize=9)

axset_turb[3].set_yticks(axset_turb[0].get_yticks())
axset_turb[3].set_ylim(-100, 3300)
#axset_turb[3].set_xlabel(r"$\bar{w'w'w'}$", fontsize=12)
axset_turb[3].set_xlabel(r"$S_w$", fontsize=12)
axset_turb[3].set_xlim(-1.3, 3.1)
axset_turb[3].set_xticks([-1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3])
axset_turb[3].set_xticklabels(['-1.2', '', '0', '', '1.2', '', '2.4', ''], fontsize=9)

for ax in axset_turb:
    ax.set_ylim(0, 3200)
    ax.set_yticks(np.arange(0, 3200, 500))
axset_turb[0].set_yticklabels(ax.get_yticks().astype(str), fontsize=9)
axset_turb[0].set_ylabel('altitude (m)', fontsize=12)
for ax in axset_turb[1:]: 
    ax.set_yticklabels(['' for t in ax.get_yticks()], fontsize=9)
    
fig.savefig("./fig_wind+turb_profiles_CAM.png")
##_____________________________________________________________________________
## Turbulence profiles

    
 
## Flux profiles
##_____________________________________________________________________________
varkeys_flux = ['WPTHLP_CLUBB', 'WPRTP_CLUBB', 'WPTHVP_CLUBB']
fig, axset_flux = plt.subplots(1, 4, figsize=(10, 4))

# Plot CAM:
#for n in ncld_g1: plotflux_singlecldmod(str(n).zfill(2), 'red', axset, varkeys_flux)
#for n in ncld_g2: plotflux_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys_flux)
#for n in ncld_g3: plotflux_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys_flux)
plot_meanprf_ilevs(ncld_g3, 'blue', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.2)
plot_meanprf_ilevs(ncld_g2, 'grey', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.4)
plot_meanprf_ilevs(ncld_g1, 'red', axset_flux[[0,1]], varkeys_flux[slice(0,2)], 0.4)
plot_meanprf_mlevs(ncld_g3, 'blue', [axset_flux[2]], [varkeys_flux[2]], 0.2)
plot_meanprf_mlevs(ncld_g2, 'grey', [axset_flux[2]], [varkeys_flux[2]], 0.4)
plot_meanprf_mlevs(ncld_g1, 'red', [axset_flux[2]], [varkeys_flux[2]], 0.4)

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








