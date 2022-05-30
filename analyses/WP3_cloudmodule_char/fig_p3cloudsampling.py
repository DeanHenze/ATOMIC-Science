# -*- coding: utf-8 -*-
"""
Created on Mon May 30 07:46:23 2022

@author: Dean

Figure for P-3 cloud sampling locations.
"""



# Built in
import os

# Third party
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs



# Set basemap:
fig = plt.figure(figsize=(8,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_xlim()
ax.set_extent([-61, -48, 7.5, 18], crs=None)
ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)


# Plot P-3 cloud module locations:
path_p3insitu = "./cldmod_datafiles/"
fnames_p3insitu = [f for f in os.listdir(path_p3insitu) 
                   if '_insitu+remote_' in f]

ncld_list = [str(i).zfill(2) for i in np.arange(1, 17, 1)]
for ncld in ncld_list:
    
    # Load p3 cloud module data:
    fname_p3 = [f for f in fnames_p3insitu if "ncld%s" % ncld in f][0]
    p3 = xr.load_dataset(path_p3insitu + fname_p3)

    # Plot mean location:
    p3mean = p3.mean(dim='time')
    ax.scatter(p3mean['lon'], p3mean['lat'], marker='^', color='black')
    
    
"""   
# Plot mean surface winds during P-3 sampling period:
path_era5dir = "../../data/reanalysis/ECMWF_ERA5/"

p3daterange = list(pd.date_range(start='1/17/2020', end='2/12/2020'))
p3daterange = [dt.strftime('%Y%m%d') for dt in p3daterange]
for date in p3daterange:
    fname_era5 = "ECMWF_ERA5_plevs_hourly_ATOMIC-region_%i.nc" % date
    era5 = xr.load_dataset(path_era5dir + fname_era5)
"""


fig.savefig("./fig_p3samplingmap.png")









