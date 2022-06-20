# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 13:28:09 2022

@author: Dean

Synoptic time series of low level divergence, surface windspeed & direction, 
and EIS.

Other options:
    - windsheer over mixed layer and trade cumulus layer.
"""



# Built in
import os

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt



path_era5dir = r"../../data/reanalysis/ECMWF_ERA5/"
fnames_era5 = [f for f in os.listdir(path_era5dir) 
               if "_plevs_hourly" in f and ".nc" in f]
fnames_era5.sort()


latbnds = (10, 16.5)
lonbnds = (-60, -48)


date = []

div_mean = []
div_std = []
div_q25p = []
div_q75p = []

usfc_mean = []
usfc_std = []
usfc_q25p = []
usfc_q75p = []

phi_mean = []
phi_std = []
phi_q25p = []
phi_q75p = []


for f in fnames_era5:
    era5 = xr.load_dataset(os.path.join(path_era5dir+f))
    
    # data trimmed for ATOMIC region:
    era5_cldmod = era5.sel(
        latitude = slice(latbnds[0], latbnds[1]), 
        longitude = slice(lonbnds[0], lonbnds[1]), 
        )

    # date:
    dtime0 = pd.Timestamp(era5['time'].values[0])
    date.append(str(dtime0.date()))
    
    # Low level divergence:
    div_950hpa = era5['d'].sel(level=950)
    div_mean.append(div_950hpa.mean().item())
    div_std.append(div_950hpa.std().item())
    div_q25p.append(np.quantile(era5['d'], 0.25))
    div_q75p.append(np.quantile(era5['d'], 0.75))

    # Surface wind speed:
    wspd = lambda u, v: (u**2 + v**2)**0.5
    usfc = wspd(era5['u'].sel(level=1000), era5['v'].sel(level=1000))
    usfc_mean.append(usfc.mean().item())
    usfc_std.append(usfc.std().item())
    usfc_q25p.append(np.quantile(usfc, 0.25))
    usfc_q75p.append(np.quantile(usfc, 0.75))

    # Surface wind direction:    
    wdir = lambda u, v: np.arctan(u/v)*(180/np.pi)
    phi = wdir(era5['u'].sel(level=1000), era5['v'].sel(level=1000))
    phi_mean.append(phi.mean().item())
    phi_std.append(phi.std().item())
    phi_q25p.append(np.quantile(phi, 0.25))
    phi_q75p.append(np.quantile(phi, 0.75))
    
    
fig, axset = plt.subplots(3, 1, figsize=(6.5, 4))

axset[0].fill_between(date, div_q25p, div_q75p, color='grey')
axset[0].plot(date, div_mean, color='darkgrey')

axset[1].fill_between(date, usfc_q25p, usfc_q75p, color='grey')
axset[1].plot(date, usfc_mean, color='darkgrey')

axset[2].fill_between(date, phi_q25p, phi_q75p, color='grey')
axset[2].plot(date, phi_mean, color='darkgrey')


fig.savefig("./fig_synoptic_timeseries.png")








