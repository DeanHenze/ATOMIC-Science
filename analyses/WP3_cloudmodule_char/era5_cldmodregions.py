# -*- coding: utf-8 -*-
"""
Created on Wed May 25 08:49:40 2022

@author: Dean


Get ERA5 data for latitude, longitude, and time windows centered on P-3 
mean location / time of cloud module sampling.

A call to main will generate datafiles for all P-3 cloud modules. 
Directory / file paths to load and save data are specified in the call to main.
"""



import pandas as pd
import xarray as xr



def createfiles(path_era5dir, path_cldmodtab, path_savedir, 
                dlat=2, dlon=2, dt=3.1):
    """
    Get ERA5 data for latitude, longitude, and time windows centered on P-3 
    mean location / time of cloud module sampling.
    
    Assumes ERA5 filename format: 
        'ECMWF_ERA5_plevs_hourly_ATOMIC-region_*date*.nc"\'
        
    Inputs
    ------
    path_era5dir, path_cldmodtab, path_savedir: str's.
        Paths to the ERA5 datafile directory, P-3 cloud module table, and 
        directory to save results to.
        
    dlat, dlon, dt: scalars.
        Window half width of latitude (degress), longitude (degrees), and 
        time (hours) to get data for.
    """
    
    cldmod_tab = pd.read_csv(path_cldmodtab)
    cldmod_tab['start_datetime'] = pd.to_datetime(cldmod_tab['start_datetime'])
    cldmod_tab['end_datetime'] = pd.to_datetime(cldmod_tab['end_datetime'])
    
    # Convert hour interval to pandas.Timedelta:
    dt = pd.Timedelta(hours=dt)

    for i, row in cldmod_tab.iterrows():
        
        print(row['flight_date'])
            
        # mean time of P-3 cloud module:
        dt_cldmod = row['end_datetime'] - row['start_datetime']
        tmean = row['start_datetime'] + dt_cldmod/2
        
        # Load ERA5 data for flight date:
        fname_era5 = "ECMWF_ERA5_plevs_hourly_ATOMIC-region_%i.nc" % row['flight_date']
        era5 = xr.load_dataset(path_era5dir + fname_era5)
    
        # data in the spatial / temporal region near cloud module:
        era5_cldmod = era5.sel(
            latitude = slice(row['lat_mean']+dlat, row['lat_mean']-dlat), 
            longitude = slice(row['lon_mean']-dlon, row['lon_mean']+dlon), 
            time = slice(tmean-dt, tmean+dt)
            )
        
        # Save new file:
        ncld = str(int(row['num_cld_iop'])).zfill(2)
        fnamesave = "p3cld_ECMWF_ERA5_plevs_hourly_ncld%s.nc" % ncld
        era5_cldmod.to_netcdf(path_savedir + fnamesave)
        
        
        
if __name__=="__main__":
    
    path_era5dir = "../../data/reanalysis/ECMWF_ERA5/"
    path_cldmodtab = "./p3_cloudmodules.csv"
    path_savedir = "./cldmod_datafiles/"
    
    createfiles(
        path_era5dir, path_cldmodtab, path_savedir,
        dlat=2, dlon=2, dt=1.5
        )






