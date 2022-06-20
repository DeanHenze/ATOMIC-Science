# -*- coding: utf-8 -*-
"""
Created on Sun May 29 12:41:53 2022

@author: Dean

Figure of mososcale and synoptic characteristics for each P-3 cloud module.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Local code
import thermo



path_era5dir = "./era5_cldmods/"
path_p3insitu = "../WP3_cloudmodule_char/cldmod_datafiles/"
fnames_p3insitu = [f for f in os.listdir(path_p3insitu) 
                   if '_insitu+remote_' in f]


varkeys = ['LTS', 'wspd_sfc', 'wdir_sfc', 
           'wdiv_950hPa', 'wdiv_850hPa', 'wdiv_750hPa'
           ]


# Wind related functions:
windspeed = lambda u, v: (u**2 + v**2)**0.5
winddirection = lambda u, v: np.arctan(u/v)*(180/np.pi)
windshear = lambda u1, u2, v1, v2, dz: (1/dz)*((u2-u1)**2 + (v2-v1)**2)**0.5


ncld_g1 = [1, 5, 4, 8]
ncld_g2 = [7, 9, 11, 10, 12, 6]
ncld_g3= [15, 3, 2, 13, 16, 14]


def div_ts(ncld, pltcolor, axset):
    """
    """
    # Load data:
    ncld_str = str(ncld).zfill(2)
    fname_era5 = "p3cld_ECMWF_ERA5_plevs_hourly_ncld%s.nc" % ncld_str
    era5 = xr.load_dataset(path_era5dir + fname_era5)
    fname_p3 = [f for f in fnames_p3insitu if "ncld%s" % ncld_str in f][0]
    p3 = xr.load_dataset(path_p3insitu + fname_p3)
    
    # Compute additional vars for ERA5:
    era5['theta'] = thermo.theta(era5['t'], era5['level']*100)
    era5_ns = era5.sel(level=1000) # near surface
    era5['wdiv_950hPa'] = era5['d'].sel(level=950)
    era5['wdiv_850hPa'] = era5['d'].sel(level=850)
    era5['wdiv_750hPa'] = era5['d'].sel(level=750)
    era5['LTS'] = era5['theta'].sel(level=700) - era5_ns['theta']
    era5['qlapse'] = era5['q'].sel(level=700) - era5_ns['q']
    era5['wspd_sfc'] = windspeed(era5_ns['u'], era5_ns['v'])
    era5['wdir_sfc'] = winddirection(era5_ns['u'], era5_ns['v'])
        # Windshear between surface and ~3km:
    era5['windshear_sfc700hpa'] = windshear(
        era5['u'].sel(level=1000), era5['u'].sel(level=700), 
        era5['v'].sel(level=1000), era5['v'].sel(level=700), 
        (1000-700)
        )
    
    
    # ERA5 grid point nearest the mean lat/lon of P-3 cloud module:
    #p3mean = p3.mean(dim='time')
    #era5_nearp3 = era5.sel(
    #    longitude=p3mean['lon'], latitude=p3mean['lat'],
    #    method='nearest'
    #    )
    #era5_nearp3 = era5_nearp3.mean(dim='time')
    

    #era5std = era5.std(dim=['longitude','latitude','time'])    
    era5mean_ts = era5.mean(dim=['longitude','latitude'])
    dt = np.flip(era5mean_ts['time'][-1] - era5mean_ts['time'])/(60*10**9)
    dt.values = dt.values.astype(float)
    era5mean_ts = era5mean_ts.where(dt<=120, drop=True)
    dt = dt[dt<=120]
    
    axset[0].plot(dt, era5mean_ts['wdiv_950hPa'], color=pltcolor)
    axset[1].plot(dt, era5mean_ts['wspd_sfc'], color=pltcolor)
    axset[2].plot(dt, era5mean_ts['windshear_sfc700hpa'], color=pltcolor)
    axset[3].plot(dt, era5mean_ts['LTS'], color=pltcolor)
    #axset[4].plot(dt, era5mean_ts['qlapse'], color=pltcolor)
    

## Create plot
##_____________________________________________________________________________
# Plot:
fig, axset = plt.subplots(4, 1, figsize=(6.5, 7))
ncldgroups = [ncld_g1, ncld_g2, ncld_g3]
pltcolors = ['red', 'grey', 'blue']
for ngroup, c in zip(ncldgroups, pltcolors):
    for n in ngroup:
        div_ts(n, c, axset)


# Axes labels, limits, legend:
for ax in axset: 
    ax.invert_xaxis()
    ax.set_xticklabels((ax.get_xticks()*-1).astype(int).astype(str))

axset[-1].set_xlabel('hours before P-3 sampling', fontsize=12)

axset[0].set_ylabel(r'$div_{950}$ (s$^{-1}$ * 10$^{-5}$)', fontsize=12)
axset[0].set_yticks(10**(-5)*np.array([-1, 0, 1]))
axset[0].set_yticklabels((axset[0].get_yticks()*10**5).astype(int).astype(str))
axset[1].set_ylabel(r'$|U|_{sfc}$ (m/s)', fontsize=12)
axset[2].set_ylabel(r'$(dU/dz)_{h}$ (s$^{-1}$)', fontsize=12)
axset[3].set_ylabel(r'LTS (K)', fontsize=12)

legend_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='blue', lw=4)]
axset[0].legend(
    legend_lines, ['cg1', 'cg2', 'cg3'], 
    loc='lower left', bbox_to_anchor= (0.5, 1.01), ncol=3,
    borderaxespad=0, frameon=False
    )


fig.savefig("./fig_mesoscale_timeseries.png")
##_____________________________________________________________________________
## Create plot














