# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:28:09 2022

@author: Dean

Synoptic time series of 600 mb vertical velocity, surface windspeed, 
wind shear in the lower 3 km, surface pressure, and lower tropospheric stability. 

Other options:
    - windsheer over mixed layer and trade cumulus layer.
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
fnames_era5 = [f for f in os.listdir(path_era5dir) 
               if "_plevs_hourly" in f and ".nc" in f]
fnames_era5.sort()

path_cldmodtab = r"../WP3_cloudmodule_char/p3_cloudmodules.csv"
##_____________________________________________________________________________
## I/O paths and filenames


## Load cloud module table:
cldmods = pd.read_csv(path_cldmodtab)


windshear = lambda u1, u2, v1, v2, dz: (1/dz)*((u2-u1)**2 + (v2-v1)**2)**0.5
wspd = lambda u, v: (u**2 + v**2)**0.5


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


def extra_era5vars(era5):
    """
    LTS, surface windspeed, windshear, 600 hPa vertical velocity.
    """
    era5['usfc'] = wspd(era5['u'].sel(level=1000), era5['v'].sel(level=1000))
    era5['theta'] = thermo.theta(era5['t'], era5['level']*100)
    era5['w_600hPa'] = era5['w'].sel(level=600)
        # Windshear between surface and ~3km units of s^-1:
    era5['windshear_sfc700hpa'] = windshear(
        era5['u'].sel(level=1000), era5['u'].sel(level=700), 
        era5['v'].sel(level=1000), era5['v'].sel(level=700), 
        (1000-700)*10 # 1000 m / 100 hPa
        )    
    era5['LTS'] = (era5['theta'].sel(level=650) - era5['theta'].sel(level=1000))   
    return era5



## Get timeseries of statistics for large scale forcing from ERA5
##_____________________________________________________________________________
## Region to get forcing values for:
latbnds = (16.5, 10)
lonbnds = (-60, -48)

# List storage for date and large scale forcing variables:
date = []
usfc_stats = []
lts_stats = []
omega_stats = []
dudz_stats = []
psfc_stats = []

for f in fnames_era5:
    
    # Load era5 data and add some variables:
    era5 = xr.load_dataset(os.path.join(path_era5dir+f))
    era5 = extra_era5vars(era5)
    
    # data trimmed for ATOMIC region:
    era5_atomic = era5.sel(
        latitude = slice(latbnds[0], latbnds[1]), 
        longitude = slice(lonbnds[0], lonbnds[1]), 
        )
    
    # Compute and append date and forcing statistics:
        # date:
    dtime0 = pd.Timestamp(era5['time'].values[0])
    date.append(str(dtime0.date()))
        # Surface wind speed:
    usfc_stats.append(getstats(era5_atomic['usfc']))
        # Lower tropo stability:
    lts_stats.append(getstats(era5_atomic['LTS']))
        # Upper level vertical velocity:
    omega_stats.append(getstats(era5_atomic['w_600hPa']))
        # Surface pressure:
    #psfc_stats.append(getstats(era5_atomic['w_600hPa']))
        # Wind shear:
    dudz_stats.append(getstats(era5_atomic['windshear_sfc700hpa']))
    
    era5.close()
##_____________________________________________________________________________
## Get timeseries of statistics for large scale forcing from ERA5



## ERA5 data nearest the each cloud module time/lat/lon
##_____________________________________________________________________________
# List storage for date and large scale forcing variables:
date_p3 = []
usfc_p3 = []
lts_p3 = []
omega_p3 = []
dudz_p3 = []
psfc_p3 = []

for i, row in cldmods.iterrows():
    
    # ERA5 data on P-3 flight date with extra forcing variables: 
    fname_date = [f for f in fnames_era5 if "_%i" % row['flight_date'] in f]
    fname_date = fname_date[0]
    era5 = xr.load_dataset(os.path.join(path_era5dir+fname_date))
    era5 = extra_era5vars(era5)

    # ERA5 point nearest P-3 sampling time/location:
    era5_nearp3 = era5.sel(
        time=row['start_datetime'],
        longitude=row['lon_mean'], latitude=row['lat_mean'],
        method='nearest'
        )
    
    # Append date and forcing statistics:
    dtime0 = pd.Timestamp(era5_nearp3['time'].values)
    date_p3.append(str(dtime0.date()))
    usfc_p3.append(era5_nearp3['usfc'].item())
    lts_p3.append(era5_nearp3['LTS'].item())
    omega_p3.append(era5_nearp3['w_600hPa'].item())
    #psfc_p3.append(era5_nearp3['w_600hPa'].item())
    dudz_p3.append(era5_nearp3['windshear_sfc700hpa'].item())
    
date_p3 = np.array([pd.Timestamp(d) for d in date_p3])
usfc_p3 = np.array(usfc_p3)
lts_p3 = np.array(lts_p3)
omega_p3 = np.array(omega_p3)
dudz_p3 = np.array(dudz_p3)
psfc_p3 = np.array(psfc_p3)
##_____________________________________________________________________________
## ERA5 data nearest the each cloud module time/lat/lon

    
    
## Shaded regions for ERA5 time series stats    
##_____________________________________________________________________________
fig, axset = plt.subplots(4, 1, figsize=(6.5, 6))

date = [pd.Timestamp(d) for d in date]
omega_stats = np.array(omega_stats)
axset[0].fill_between(
    date, omega_stats[:,2], omega_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[0].plot(date, omega_stats[:,0], color='black')
    # Reference line at y=0:
axset[0].hlines(
    0, min(axset[0].get_xticks()), max(axset[0].get_xticks()), 
    colors='black', linestyles='solid', linewidth=1, zorder=10
    )

lts_stats = np.array(lts_stats)
axset[1].fill_between(
    date, lts_stats[:,2], lts_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[1].plot(date, lts_stats[:,0], color='black')

usfc_stats = np.array(usfc_stats)
axset[2].fill_between(
    date, usfc_stats[:,2], usfc_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[2].plot(date, usfc_stats[:,0], color='black')

dudz_stats = np.array(dudz_stats)
axset[3].fill_between(
    date, dudz_stats[:,2], dudz_stats[:,3], 
    color='orange', alpha=0.3
    )
p = axset[3].plot(date, dudz_stats[:,0], color='black')
    # Reference line at y=0.002:
axset[3].hlines(
    0.002, min(axset[3].get_xticks()), max(axset[3].get_xticks()), 
    colors='black', linestyles='solid', linewidth=1, zorder=10
    )
##_____________________________________________________________________________
## Shaded regions for ERA5 time series stats  



# Axes labels, limits, legend
##_____________________________________________________________________________
for ax in axset[0:3]: 
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(['' for t in ax.get_xticks()])
    
for ax in axset:
    ax.set_xlim(pd.Timestamp("2020-01-15"), pd.Timestamp("2020-02-13"))

    
axset[0].set_ylabel(r'$\omega_{600}$ (Pa * s$^{-1}$)', fontsize=12)

axset[1].set_ylim(14, 23)
axset[1].set_yticks(np.arange(15, 22, 2))
axset[1].set_ylabel(r'LTS (K)', fontsize=12)

axset[2].set_ylim(2, 13)
axset[2].set_yticks(np.arange(3, 13, 3))
axset[2].set_ylabel(r'$|U|_{sfc}$ (m/s)', fontsize=12)

axset[3].set_ylim(0.0007, 0.0037)
axset[3].set_yticks([0.001, 0.002, 0.003])
axset[3].set_ylabel(r'$(dU/dz)_{h}$ (s$^{-1}$)', fontsize=12)

xticklabels = [
    text.Text(*ticklab.get_position(), ticklab.get_text()[5:])
    for ticklab in axset[3].get_xticklabels()
    ]
axset[3].set_xticks(axset[3].get_xticks())
axset[3].set_xticklabels(xticklabels, rotation=45, rotation_mode=None, 
                         horizontalalignment='right')
#axset[3].tick_params(axis='x', labelrotation=30, labelright=True)
axset[3].set_xlabel('date', fontsize=12)
##_____________________________________________________________________________
# Axes labels, limits, legend



## Annotations for P-3 cloud modules
##_____________________________________________________________________________
# Annotations at top of graph:
cld_annot = [
    '01', '02\n03', '04\n05', '06\n07\n08', 
    '09\n10\n11', '12\n13', '14\n15', '16'
    ]
ha_locs = np.array([pd.Timestamp(str(d)) for d in cldmods['flight_date'].unique()])
ha_offsets = np.array([
    pd.Timedelta(0), pd.Timedelta(0), pd.Timedelta(0), pd.Timedelta(0), 
    pd.Timedelta(days=0.5), -pd.Timedelta(days=0.5), pd.Timedelta(0), 
    pd.Timedelta(days=0.5)
    ])
for xpos, x_off, txt in zip(ha_locs, ha_offsets, cld_annot):
    axset[0].text(
        xpos+x_off, 0.19, txt, fontsize=8, 
        ha='center', va='bottom'
        )
    
# Vertical lines at each day with cloud modules:
for ax in axset:
    ylims = ax.get_ylim()
    ax.vlines(
        ha_locs + ha_offsets, 
        np.ones(len(ha_locs))*ylims[0], np.ones(len(ha_locs))*ylims[1], 
        colors='red', linewidths=0.5
        )
##_____________________________________________________________________________
## Annotations for P-3 cloud modules



fig.savefig("./fig_synoptic_timeseries.png")








