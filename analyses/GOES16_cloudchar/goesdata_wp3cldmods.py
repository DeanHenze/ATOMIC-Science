# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:55:39 2021

@author: Dean


Get GOES-16 image data each P-3 cloud module. 

For each cloud module, the GOES image taken is the nearest 20 minute timestamp 
to the mean time of the module, and in a 2 deg x 2 deg lat, lon region about 
the P-3's mean location.

Reflectance, IR temperature, and cloud top data are taken. Also, an IR 
threshold mask is appended to the dataset.
"""


# Built in
import os

# Third party
import pandas as pd
import xarray as xr

# Local code
import datetimeconverter as dtconv
import threshmask



## I/O paths
##_____________________________________________________________________________
path_cldmodtab = "../WP3_cloudmodule_char/p3_cloudmodules.csv"
path_goesdir = "../../data/GOES-16/GOES16_for_cloudmodules/"
path_savedir = "./goes16_WP3cldmods/"
##_____________________________________________________________________________
## I/O paths



def goesfile_neartime(path_goesdir, t1, t2, varsubset=None):
    """
    Get a single GOES .nc data file for the mean time of t1 and t2 rounded to 
    the nearest 20 minutes.
    
    Inputs
    ------
    path_goesdir: str or path.
        Path to directory containing GOES data files.
    
    t1, t2: pd.Timestamp objects.
        t2 should be greater than t1.
        
    varsubset: iterables of str's.
        Optional to return only these variables keys. 

    Returns
    -------
    xarray.dataset
    """    
    t1 = pd.Timestamp(t1)
    t2 = pd.Timestamp(t2)
        
    # Mean datetime of cloud module rounded to nearest 20 min:
    dtime = (t1 + (t2 - t1)/2).round('min')
    dtime += pd.Timedelta(10, 'min') # Add 10 min so next line floors instead of rounds.
    dtime -= pd.Timedelta(dtime.minute % 20, 'min') # nearest 20 min.    
    
    # Load corresponding GOES-16 file and isolate a few variables:
    goesdata, fname = goesfile(path_goesdir, dtime)
    if varsubset is not None:
        goesdata = goesdata[varsubset]

    return goesdata, fname



def goesfile(path_goesdir, dtime):
    """
    Return GOES16 data file for specified date and time.
    
    dtime: str or pandas Timestamp object.
        If str, it should be in the format 'yyyymmddHHMM' (y=year, m=month, 
        d=day, H=hour, M = minute).
    """
    
    if type(dtime) == pd._libs.tslibs.timestamps.Timestamp:
        dtime = dtime.strftime("%Y%m%d%H%M")
        
    date = dtime[0:8]
    time = dtime[8:12]
    
    nday = dtconv.yyyymmdd_to_ndays(date, with_year=True) # day num since Jan1, 2020.
    fname = "G16V04.0.ATOMIC.%s.%s.PX.02K.nc" % tuple([nday, time])
    return xr.load_dataset(path_goesdir + fname), fname



def clip_region(data, latcen, loncen, dlat, dlon, 
                latkey='latitude', lonkey='longitude'):
    """
    Return subset of passed geospatial dataset in specified lat / lon region.
    
    Inputs
    ------
    data: xr.Dataset
    
    latcen, loncen: scalars.
        Center latitude, longitude of subset region.
        
    dlat, dlon: scalars.
        Half widths of desired latitude, longitude region.
        
    latkey, lonkey:
        Keys in 'data' for latitude and longitude.
    """
    data_subset = data.where(
        (data[lonkey]>(loncen-dlon)) & (data[lonkey]<(loncen+dlon)), 
        drop=True
        )
    data_subset = data_subset.where(
        (data_subset[latkey]>(latcen-dlat)) & (data_subset[latkey]<(latcen+dlat)), 
        drop=True
        )
    
    return data_subset



if __name__=="__main__":

    # P-3 cloud modules table (info on lat, lon, time of each module):    
    cldmodtab = pd.read_csv(path_cldmodtab)
    
    
    # Things that may want to be modified:
    dlat = 1.25 # half-width of latitude for desired GOES region.
    dlon = 1.25 # half-width of longitude for desired region.
    varsubset = [ # subset of GOES variables to save.
        'reflectance_vis', 'temperature_ir', 
        'cloud_top_height', 'cloud_top_temperature'
        ]
    

    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)   
    
    for i, row in cldmodtab.iterrows():
        
        ncld = str(row['num_cld_iop']).zfill(2)
        print("Getting goes image data for P-3 cloud module %s" % ncld)
        
        # Data nearest P-3 cloud module time:
        t1 = row['start_datetime']
        t2 = row['end_datetime']
        goesdata, fname_goes = goesfile_neartime(path_goesdir, t1, t2, varsubset)
        
        # Subset of data near cloud module sampling lat / lon:
        latcen = row['lat_mean']
        loncen = row['lon_mean']
        goesdata_sub = clip_region(goesdata, latcen, loncen, dlat, dlon)
        
        # IR threhold mask:
        lowT = 275 # IR temperature in K
        highT = 292
        goesdata_sub['ir_mask'] = threshmask.threshmask(
            goesdata_sub['temperature_ir'], low_thresh=lowT, high_thresh=highT)
        goesdata_sub['ir_mask'].attrs['long_name'] = (
            "Mask = 1 if IR temperature is between the low and high thresholds."
            )
        goesdata_sub['ir_mask'].attrs['low_threshold'] = lowT
        goesdata_sub['ir_mask'].attrs['high_threshold'] = highT
        
        # Save file:
        fname_save = fname_goes[:-3] + "_ncld" + ncld + ".NC"
        goesdata_sub.to_netcdf(path_savedir + fname_save)