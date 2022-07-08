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
path_camdir = "./CAM_WP3cldmod_locs/"
path_p3drpsnd_dir = "../WP3_cloudmodule_char/cldmod_datafiles/"



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
    


#ncld = [str(n).zfill(2) for n in [1,2,3,4,5,6,7,8,9,10,11,12,13]]
#for n in ncld: plot_singlecldmod(n)

varkeys = ['theta', 'H216OV']
fig, axset = plt.subplots(1, 2, figsize=(10, 4))
for n in ncld_g1: plotthermo_singlecldmod(str(n).zfill(2), 'red', axset, varkeys)
for n in ncld_g2: plotthermo_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys)
for n in ncld_g3: plotthermo_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys)



varkeys = ['WP2_CLUBB', 'WP3_CLUBB', 'WPRTP_CLUBB', 'WPTHLP_CLUBB']
fig, axset = plt.subplots(1, 4, figsize=(10, 4))
for n in ncld_g1: plotflux_singlecldmod(str(n).zfill(2), 'red', axset, varkeys)
for n in ncld_g2: plotflux_singlecldmod(str(n).zfill(2), 'grey', axset, varkeys)
for n in ncld_g3: plotflux_singlecldmod(str(n).zfill(2), 'blue', axset, varkeys)















