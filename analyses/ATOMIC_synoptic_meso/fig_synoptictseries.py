# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:28:09 2022

@author: Dean

Synoptic time series of 600 mb vertical velocity, surface windspeed, 
wind shear in the lower 3 km, surface pressure, and lower tropospheric 
stability. 
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
    era5['gph_700hPa'] = era5['z'].sel(level=700)/9.8 # divide by gravitational accel.
    era5['gph_500hPa'] = era5['z'].sel(level=500)/9.8
    return era5



## Get timeseries of statistics for large scale forcing from ERA5
##_____________________________________________________________________________
## Region to get forcing values for:
latbnds = (16.5, 10)
lonbnds = (-60, -48)

# List storage for date and large scale forcing variables:
date = []
date_2dvars = [] # To make sure it matches 'date' for 3d vars.
usfc_stats = []
lts_stats = []
omega_stats = []
dudz_stats = []
psfc_stats = []
gph_stats = []

# ERA5 3D vars:
for f in fnames_era5plevs:
    
    # Load era5 3D data then add some variables:
    era5_3d = xr.load_dataset(os.path.join(path_era5dir+f))
    era5_3d = extra_era5vars(era5_3d)
    
    # data trimmed for ATOMIC region:
    era5_3d_atomic = era5_3d.sel(
        latitude = slice(latbnds[0], latbnds[1]), 
        longitude = slice(lonbnds[0], lonbnds[1]), 
        )
    
    # Compute and append date and forcing statistics:
        # date:
    dtime0 = pd.Timestamp(era5_3d['time'].values[0])
    date.append(str(dtime0.date()))
        # Surface wind speed:
    usfc_stats.append(getstats(era5_3d_atomic['usfc']))
        # Lower tropo stability:
    lts_stats.append(getstats(era5_3d_atomic['LTS']))
        # Upper level vertical velocity:
    omega_stats.append(getstats(era5_3d_atomic['w_600hPa']))
        # Wind shear:
    dudz_stats.append(getstats(era5_3d_atomic['windshear_sfc700hpa']))
        # 500 hPa geopotential height:
    gph_stats.append(getstats(era5_3d_atomic['gph_700hPa']))
    
    era5_3d.close()
    
    
# ERA5 2D vars:
for f in fnames_era5singlelev:
    
    # Load era5 3D data then add some variables:
    era5_2d = xr.load_dataset(os.path.join(path_era5dir+f))
    
    # data trimmed for ATOMIC region:
    era5_2d_atomic = era5_2d.sel(
        latitude = slice(latbnds[0], latbnds[1]), 
        longitude = slice(lonbnds[0], lonbnds[1]), 
        )
    
    # Compute and append date and forcing statistics:
        # date:
    dtime0 = pd.Timestamp(era5_2d['time'].values[0])
    date_2dvars.append(str(dtime0.date()))
        # Surface pressure:
    psfc_stats.append(getstats(era5_2d_atomic['sp']/100)) # 1/100 to convert to hPa.
        
    era5_2d.close()
    

# Verify that date lists for 3D and 2D vars match:
if date != date_2dvars:
    print("Warning: date lists from the 2D and 3D ERA5 vars do not match. "
          "This could result in an incorrect surface pressure time series.")
##_____________________________________________________________________________
## Get timeseries of statistics for large scale forcing from ERA5



## ERA5 data nearest the each cloud module time/lat/lon
##_____________________________________________________________________________
"""
# List storage for date and large scale forcing variables:
date_p3 = []
usfc_p3 = []
lts_p3 = []
omega_p3 = []
dudz_p3 = []
psfc_p3 = []

for i, row in cldmods.iterrows():
    
    # ERA5 pressure level data on P-3 flight date with extra forcing variables: 
    f_plev = [f for f in fnames_era5plevs if "_%i" % row['flight_date'] in f]
    f_plev = f_plev[0]
    era5_plev = xr.load_dataset(os.path.join(path_era5dir+f_plev))
    era5_plev = extra_era5vars(era5_plev)
    
    # Get ERA5 surface pressure from single level data on P-3 flight date: 
    f_singlelev = [f for f in fnames_era5singlelev if "_%i" % row['flight_date'] in f]
    f_singlelev = f_singlelev[0]
    era5_psfc = xr.load_dataset(os.path.join(path_era5dir+f_singlelev))
    

    # ERA5 point nearest P-3 sampling time/location:
    era5_nearp3 = era5_plev.sel(
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
"""
##_____________________________________________________________________________
## ERA5 data nearest the each cloud module time/lat/lon

    
    
## Shaded regions for ERA5 time series stats    
##_____________________________________________________________________________
fig, axset = plt.subplots(4, 1, figsize=(6.5, 6))

date_ts = [pd.Timestamp(d) for d in date]

# Vertical velocity:
omega_stats = np.array(omega_stats)
axset[0].fill_between(
    date_ts, omega_stats[:,2], omega_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[0].plot(date_ts, omega_stats[:,0], color='black')
    # Reference line at y=0:
axset[0].hlines(
    0, min(axset[0].get_xticks()), max(axset[0].get_xticks()), 
    colors='black', linestyles='solid', linewidth=1, zorder=10
    )

# LTS:
lts_stats = np.array(lts_stats)
axset[1].fill_between(
    date_ts, lts_stats[:,2], lts_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[1].plot(date_ts, lts_stats[:,0], color='black')

# Surface windspeed:
usfc_stats = np.array(usfc_stats)
axset[2].fill_between(
    date_ts, usfc_stats[:,2], usfc_stats[:,3], 
    color='orange', alpha=0.3
    )
axset[2].plot(date_ts, usfc_stats[:,0], color='black')

# Wind shear:
psfc_stats = np.array(psfc_stats)
axset[3].fill_between(
    date_ts, psfc_stats[:,2], psfc_stats[:,3], 
    color='orange', alpha=0.3
    )
p = axset[3].plot(date_ts, psfc_stats[:,0], color='black')

# Wind shear:
#dudz_stats = np.array(dudz_stats)
#axset[4].fill_between(
#    date_ts, dudz_stats[:,2], dudz_stats[:,3], 
#    color='orange', alpha=0.3
#    )
#p = axset[4].plot(date_ts, dudz_stats[:,0], color='black')





#gph_stats = np.array(gph_stats)
#axset[4].fill_between(
#    date_ts, gph_stats[:,2], gph_stats[:,3], 
#    color='orange', alpha=0.3
#    )
#p = axset[4].plot(date_ts, gph_stats[:,0], color='black')
##_____________________________________________________________________________
## Shaded regions for ERA5 time series stats  



# Axes labels, limits, legend
##_____________________________________________________________________________
for ax in axset[0:3]: 
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(['' for t in ax.get_xticks()])
    
for ax in axset:
    ax.set_xlim(pd.Timestamp("2020-01-15"), pd.Timestamp("2020-02-13"))

    
axset[0].set_ylabel(r'$\omega_{600}$ (Pa s$^{-1}$)', fontsize=11, labelpad=8)

axset[1].set_ylim(14, 23)
axset[1].set_yticks(np.arange(15, 22, 2))
axset[1].set_ylabel(r'LTS (K)', fontsize=11, labelpad=20)

axset[2].set_ylim(2, 13)
axset[2].set_yticks(np.arange(3, 13, 3))
axset[2].set_ylabel(r'$|U|_{sfc}$ (m/s)', fontsize=11, labelpad=20)

axset[3].set_ylim(1011.5, 1019.5)
axset[3].set_yticks([1012, 1015, 1018])
axset[3].set_ylabel(r'$P_{sfc}$ (hPa)', fontsize=11, labelpad=8)

#axset[4].set_ylim(0.0007, 0.0045)
#axset[4].set_yticks([0.001, 0.002, 0.003, 0.004])
#axset[4].set_ylabel(r'dU/dz (s$^{-1}$)', fontsize=11)

#xticklabels = [
#    text.Text(*ticklab.get_position(), ticklab.get_text()[5:])
#    for ticklab in axset[4].get_xticklabels()
#    ]
#axset[4].set_xticks(axset[4].get_xticks())
#axset[4].set_xticklabels(xticklabels, rotation=45, rotation_mode=None, 
#                         horizontalalignment='right')
##axset[4].tick_params(axis='x', labelrotation=30, labelright=True)
#axset[4].set_xlabel('2022 date', fontsize=11)

xticklabels = [
    text.Text(*ticklab.get_position(), ticklab.get_text()[5:])
    for ticklab in axset[3].get_xticklabels()
    ]
axset[3].set_xticks(axset[3].get_xticks())
axset[3].set_xticklabels(xticklabels, rotation=45, rotation_mode=None, 
                         horizontalalignment='right')
#axset[4].tick_params(axis='x', labelrotation=30, labelright=True)
axset[3].set_xlabel('2022 date', fontsize=11)
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
    
# Vertical lines at each day with cloud modules, colored by cloud group:
linecolors = ['red', 'blue', 'red', 'grey', 'grey', 'blue', 'blue', 'blue']
for ax in axset:
    ylims = ax.get_ylim()
    ax.vlines(
        ha_locs + ha_offsets, 
        np.ones(len(ha_locs))*ylims[0], np.ones(len(ha_locs))*ylims[1], 
        colors=linecolors, linewidths=1.
        )
##_____________________________________________________________________________
## Annotations for P-3 cloud modules



fig.savefig("./fig_synoptic_timeseries.png")



## Correlations between large scale variables
##_____________________________________________________________________________
df_means = pd.DataFrame({
    'omega': omega_stats[:,0], 
    
    })
##_____________________________________________________________________________





