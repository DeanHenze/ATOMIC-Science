# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:49:40 2022

@author: Dean
"""


import numpy as np
import pandas as pd
import xarray as xr



def getwindow(ds, dimkeys, dimcenters, dimhalfwidths):
    """
    """
    dimcenters, dimhalfwidths = np.array(dimcenters), np.array(dimhalfwidths)
    windows = [slice(xc - dx, xc + dx) for xc, dx in zip(dimcenters, dimhalfwidths)] 
    slicedict = dict(zip(dimkeys, windows))
    return ds.sel(slicedict)    
    


def era5_cldmodsingle(era5, latcen, loncen, dlat, dlon):
    """
    Get ERA5 data in the region around specified lat / lon center.
    
    era5: xarray.Dataset.
        Latitude and longitude must be two of the dimensions.
    """
    return era5.sel(
        lat=slice(latcen-dlat, latcen+dlat), 
        lon=slice(loncen-dlon, loncen+dlon)
        )



path_era5dir = "../../data/reanalysis/ECMWF_ERA5/"

cldmod_tab = pd.read_csv("./p3_cloudmodules.csv")
cldmod_tab['start_datetime'] = pd.to_datetime(cldmod_tab['start_datetime'])
cldmod_tab['end_datetime'] = pd.to_datetime(cldmod_tab['end_datetime'])

for i, row in cldmod_tab.iterrows():
        
    # mean time of P-3 cloud module:
    dthalf = (row['end_datetime'] - row['start_datetime'])/2
    tmean = row['start_datetime'] + dthalf
    
    # Load ERA5 data for flight date:
    fname_era5 = "ECMWF_ERA5_plevs_hourly_ATOMIC-region_%i.nc" % row['flight_date']
    era5 = xr.load_dataset(path_era5dir + fname_era5)

    # data in the spatial / temporal region near cloud module:
    dlat = 2                          # degrees
    dlon = 2
    dt = pd.Timedelta(hours=3.1)      # hours
    era5_cldmod = era5.sel(
        latitude = slice(row['lat_mean']+dlat, row['lat_mean']-dlat), 
        longitude = slice(row['lon_mean']-dlon, row['lon_mean']+dlon), 
        time = slice(tmean-dt, tmean+dt)
        )








