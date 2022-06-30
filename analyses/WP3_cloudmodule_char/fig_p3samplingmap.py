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



## Set basemap
##_____________________________________________________________________________
fig = plt.figure(figsize=(8,4))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=ccrs.PlateCarree())
ax.coastlines()
ax.set_xlim()
ax.set_extent([-61, -48, 7.5, 18], crs=None)
gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.xlines=False
gl.ylines=False
##_____________________________________________________________________________
## Set basemap


## Plot P-3 cloud module locations
##_____________________________________________________________________________
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
##_____________________________________________________________________________    
## Plot P-3 cloud module locations
    

## Plot mean surface winds during P-3 sampling period
##_____________________________________________________________________________
path_era5dir = "../../data/reanalysis/ECMWF_ERA5/"

p3daterange = list(pd.date_range(start='1/17/2020', end='2/12/2020'))
p3daterange = [dt.strftime('%Y%m%d') for dt in p3daterange]


def thinmap1d(data, di): # Thin a 1d array for use with plt.quiver()
    n = len(data)
    ithin = [i for i in range(n) if i%di==0]
    return data[ithin]


def thinmap2d(data, di): # Thin a 2d array for use with plt.quiver()
    nx, ny = np.shape(data)
    ixthin = [i for i in range(nx) if i%di==0]
    iythin = [i for i in range(ny) if i%di==0]
    datathin = data[ixthin]
    datathin = datathin[:, iythin]    
    return datathin


def sfcwinds_iopmean(path_era5dir, p3daterange):
    """
    Returns lat, lon, mean u-component, mean v-component.
    """
    
    # Daily average u, v components near surface:
    usfc, vsfc = [], []
    for date in p3daterange:
        fname_era5 = "ECMWF_ERA5_plevs_hourly_ATOMIC-region_%s.nc" % date
        era5 = xr.load_dataset(path_era5dir + fname_era5)
        usfc.append(era5['u'].sel(level=1000).mean(dim='time'))
        vsfc.append(era5['v'].sel(level=1000).mean(dim='time'))
        era5.close()
            
    # Average over days:
    usfc_iop = np.mean(usfc, axis=0)
    vsfc_iop = np.mean(vsfc, axis=0)
    
    return (
        usfc[0]['longitude'].values, usfc[0]['latitude'].values, 
        usfc_iop, vsfc_iop
        )


uvsfc_iop = sfcwinds_iopmean(path_era5dir, p3daterange)
uvsfc_iopthin = []
di = 4

qvr = ax.quiver(
    thinmap1d(uvsfc_iop[0], di), thinmap1d(uvsfc_iop[1], di), 
    thinmap2d(uvsfc_iop[2], di), thinmap2d(uvsfc_iop[3], di),
    color='grey'
    #scale=0.05, scale_units='width'
    )
#ax.quiverkey(qvr, -49, 9, 0.1, label='m/s')
##_____________________________________________________________________________
## Plot mean surface winds during P-3 sampling period


fig.savefig("./fig_p3samplingmap.png")









