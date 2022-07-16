# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:28:09 2022

@author: Dean

Compare surface pressure of ATOMIC study region during the low pressure drop 
over Jan 21-25 with the remainder of the observation period.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.text as text
from matplotlib.lines import Line2D

# Local code
import thermo



## I/O paths and filenames
##_____________________________________________________________________________
path_era5dir = r"../../data/reanalysis/ECMWF_ERA5/"
fnames_era5plevs = [f for f in os.listdir(path_era5dir) 
                    if "_plevs_hourly" in f and ".nc" in f]
fnames_era5plevs.sort()
fnames_era5singlelev = [f for f in os.listdir(path_era5dir) 
                        if "_single-lev_hourly" in f and ".nc" in f]
fnames_era5singlelev.sort()

path_cldmodtab = r"../WP3_cloudmodule_char/p3_cloudmodules.csv"
##_____________________________________________________________________________
## I/O paths and filenames



## Load cloud module table:
cldmods = pd.read_csv(path_cldmodtab)



def getstats(data):
    """
    Returns median, standard deviation, inner quartiles for a DataFrame or 
    Dataset.
    """
    return (
        data.median().item(),
        data.std().item(),
        np.quantile(data, 0.25),
        np.quantile(data, 0.75)
        )



## Get surface pressure data
## Lat/lon region to get data over
latbnds = (25, 10)
#latbnds = (16.5, 10)
lonbnds = (-60, -48)

dates_lowp = pd.date_range(start='1/21/2020', end='1/25/2020', freq='D').astype(str)
dates1_lowp = pd.date_range(start='1/22/2020', end='1/22/2020', freq='D').astype(str)
dates2_lowp = pd.date_range(start='1/23/2020', end='1/23/2020', freq='D').astype(str)
dates3_lowp = pd.date_range(start='1/24/2020', end='1/24/2020', freq='D').astype(str)

psfc_lowp = [] # Will hold data for low pressure time interval
psfc_lowp1 = [] # Will hold data for low pressure time interval
psfc_lowp2 = [] # Will hold data for low pressure time interval
psfc_lowp3 = [] # Will hold data for low pressure time interval
psfc = [] # Will hold data for all other dates

for f in fnames_era5singlelev:
    
    # Load era5 data:
    era5 = xr.load_dataset(os.path.join(path_era5dir+f))
    
    # data trimmed for ATOMIC region:
    era5 = era5.sel(
        latitude = slice(latbnds[0], latbnds[1]), 
        longitude = slice(lonbnds[0], lonbnds[1]), 
        )
    
    date_f = "%s-%s-%s" % tuple([f[-11:-7], f[-7:-5], f[-5:-3]])
    #if date_f in dates_lowp:
    #    psfc_lowp.append(era5['sp']/100)
    #else:
    #    psfc.append(era5['sp']/100)
    psfc.append(era5['sp']/100)        
    if date_f in dates_lowp:
        psfc_lowp.append(era5['sp']/100)
    if date_f in dates1_lowp:
        psfc_lowp1.append(era5['sp']/100)
    if date_f in dates2_lowp:
        psfc_lowp2.append(era5['sp']/100)
    if date_f in dates3_lowp:
        psfc_lowp3.append(era5['sp']/100)
        
        
## Mean surface pressure and pressure differences:
psfc_mean_lowp = xr.concat(psfc_lowp, dim='time').mean(dim='time')
psfc_mean_lowp1 = xr.concat(psfc_lowp1, dim='time').mean(dim='time')
psfc_mean_lowp2 = xr.concat(psfc_lowp2, dim='time').mean(dim='time')
psfc_mean_lowp3 = xr.concat(psfc_lowp3, dim='time').mean(dim='time')
psfc_mean = xr.concat(psfc, dim='time').mean(dim='time')
psfc_diff = psfc_mean_lowp - psfc_mean
psfc_diff1 = psfc_mean_lowp1 - psfc_mean
psfc_diff2 = psfc_mean_lowp2 - psfc_mean
psfc_diff3 = psfc_mean_lowp3 - psfc_mean
 
        
## Contour maps
fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(6.5, 3))        
cntr_lowp = ax1.contour(
    psfc_mean_lowp['longitude'], psfc_mean_lowp['latitude'],
    psfc_mean_lowp
    )
fig.colorbar(cntr_lowp, ax=ax1)

cntr_p = ax2.contour(
    psfc_mean['longitude'], psfc_mean['latitude'],
    psfc_mean
    )
fig.colorbar(cntr_p, ax=ax2)


fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(6.5, 3))        
cntr_lowp = ax1.contour(
    psfc_diff1['longitude'], psfc_diff1['latitude'],
    psfc_diff1, 
    vmin=-8, vmax=0
    )
#fig.colorbar(cntr_lowp, ax=ax1)
cntr_lowp = ax2.contour(
    psfc_diff2['longitude'], psfc_diff2['latitude'],
    psfc_diff2, 
    vmin=-8, vmax=0
    )
cntr_lowp = ax3.contour(
    psfc_diff3['longitude'], psfc_diff3['latitude'],
    psfc_diff3, 
    vmin=-8, vmax=0
    )
fig.colorbar(cntr_lowp, ax=ax3)



# Save surface pressure maps:
path_savedir = "./era5_Psfc/"
if not os.path.isdir(path_savedir): os.mkdir(path_savedir)
psfc_mean.to_netcdf(path_savedir + "ECMWF_ERA5_psfcmap_ATOMIC-region_meanIOP.nc")
psfc_diff1.to_netcdf(path_savedir + "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200122.nc")
psfc_diff2.to_netcdf(path_savedir + "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200123.nc")
psfc_diff3.to_netcdf(path_savedir + "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200124.nc")    





















