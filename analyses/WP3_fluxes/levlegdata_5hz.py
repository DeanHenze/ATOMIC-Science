# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:37:07 2022

@author: Dean

Isolate and process data for a P-3 horizontal level leg, prepping it 
for eddy covariance and spectral analyses at 5 Hz sampling rate. Processed 
data is saved as a new file, one file per level leg. The following variables 
are collected into the file:
    latitude, longitude, altitude
    horizontal wind componenets
    vertical wind
    temperature
    water mass mixing ratio
    HDO mixing ratio
    
    
Current status: need to change fxns so they expect pd.DataFrames rather than 
xr.Datasets.
"""



# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code:
import xcorr
import iso_fxns as isofxn



#def process_5hz(t1, t2, fname_wind, fname_pic_mr, fname_pic_iso):
def process_5hz(wind_50hz, mr_5hz, iso_5hz, t1, t2):
    """ 
    Returns processed wind, water, and isotope ratio data at 5Hz, collected 
    into a single xarray dataset.   

    Inputs
    ------
    wind_50hz: pandas.Dataframe.
        Wind data at 50 Hz sampling frequency, with time as the index.
    mr_5hz, iso_5hz: pandas.Dataframe's.
        Water and HDO mixing ratio at 5 Hz sampling frequency, with time 
        as the index.
    t1, t2: floats.
        Start and end times for segment to isolate.
        
    Returns
    -------
    
    """
    
    ## Process wind data
    ##_________________________________________________________________________            
    windkeys_new = ["lon", "lat", "alt", "u", "v", "w", "T"]
    wind_pro = varsubset(
        wind_50hz, t1, t2, 
        ["lon", "lat", "alt", "eastward_wind", 
         "northward_wind", "vertical_wind", "total_air_temperature"
         ], 
        keys_new=windkeys_new
        )
    wind_pro = avgtofreq(wind_pro, 5, reportna=False)
    wind_pro = perturbations(wind_pro, windkeys_new, [k+"'" for k in windkeys_new])
    ##_________________________________________________________________________            
    ## Process wind data
        
    
    ## Process water mixing ratio and isotope data
    ##_________________________________________________________________________            
    # Remove missing time values then merge datasets:
    mr_5hz = mr_5hz.where(~mr_5hz['time'].isnull(), drop=True)
    iso_5hz = iso_5hz.where(~iso_5hz['time'].isnull(), drop=True)
    mriso = mr_5hz.merge(iso_5hz, join='exact')
    
    # Compute HDO mixing ratio 
    mriso['qD'] = isofxn.qD_dD_convert('dD2qD', mriso['mmr'], mriso['dD'])

    # Process the time dimension:
    #mriso = mriso.where(~mriso['time'].isnull(), drop=True) # remove missing vals.
    tnew = [convert_time(dt64) for dt64 in mriso.time.values]
    mriso_pro = mriso.assign_coords(time=tnew)  
    mriso_pro = mriso_pro.interp(time=wind_pro['time'], method='nearest')
    
    mriso_pro = varsubset(
        mriso_pro, t1, t2, 
        ["mmr","qD"], 
        keys_new=["q","qD"]
        )
    mriso_pro = perturbations(
        mriso_pro, ["q","qD"], 
        ["q'","qD'"]
        )    
    
    # Time shift:
        # Determine time shift using max cross-correlation:
    c, xcor, n = xcorr.correlation(mriso_pro["q'"], wind_pro["w'"], 0, 35)
    tshift = c[np.argmax(xcor)]/5
        # Verification plot        
    plt.figure(figsize=(6,5))
    plt.plot(c/5, xcor, label='before shift')
    plt.xlabel("Time lag q' wrt to w' (seconds)", fontsize=14)
    plt.ylabel("Cross-correlation", fontsize=14)
        # Apply shift and realign with wind time:
    mriso_pro = mriso_pro.assign_coords(time=mriso_pro['time']+tshift)  
    mriso_pro = mriso_pro.interp(
        time=wind_pro['time'], 
        method='nearest'
        )
    c, xcor, n = xcorr.correlation(mriso_pro["q'"], wind_pro["w'"], 0, 35)
    plt.plot(c/5, xcor, label='after shift')
    plt.legend(fontsize=12)
    ##_________________________________________________________________________            
    ## Process water mixing ratio and isotope data

    
    data_merged = wind_pro.merge(mriso_pro, join='exact')
    return data_merged



def varsubset(data_xrds, t1, t2, keys_old, reportna=False, keys_new=None):
    """
    Get a subset of the variables in a dataset over a specific time interval 
    with missing values interpolated. Option to rename variables.
    
    Inputs
    ------
    data_xrds: xarray.Dataset.
        Time series data with 'time' as a dimension.
    t1, t2: floats.
        Start, end times to get data for; t2>t1.
    keys_old: list of str's.
        Current keys in data_xrds for variables to keep.
    keys_new: Optional list of str's.
        Rename variables with these labels. If specified, should be same 
        length as keys_old.
        
    Returns
    -------
    data_tintv: xarray.Dataset.
    """
    # Subset of vars cut to values within the time interval:
    data_tintv = data_xrds.sel(time=slice(t1, t2))
    data_tintv = data_tintv[keys_old]
    
    # Interpolate any missing values:
    if reportna:
        print("Number of missing values\n"
              "--------------------------")
        print(data_tintv.isnull().sum())
    data_tintv = data_tintv.interpolate_na(dim='time', method='linear')
    
    if keys_new is not None:
        data_tintv = data_tintv.rename(dict(zip(keys_old, keys_new)))
    
    return data_tintv
    
    
    
def avgtofreq(data_xrds, frq, reportna=False):
    """
    Average a dataset to a specified frequency. Interpolate any NAN values.
    
    Inputs
    ------
    data_xrds: xarray.Dataset.
        Dataset with 'time' as a coordinate.
    frq: float.
        Frequency (Hz) to average data to. 
        e.g. a timestamp every 1/frq seconds.
    reportna: bool.
        If True, print number of NAN values after averaging but before 
        interpolation.
        
    Returns
    -------
    data_avg: xarray.Dataset.
    """
    # Dataset time values rounded to input frequency:
    dt_frq = 1/frq
    t_frq = np.round(data_xrds['time']/dt_frq)*dt_frq
    # Average duplicate time stamps to get desired result:
    data_avg = data_xrds.groupby(t_frq).mean()
        
    return data_avg



def perturbations(data_xrds, varkeys, purtkeys):
    """
    Computes purturbations from the mean for time series data. The mean is 
    taken over the entire time series. Purturbations are added as new 
    variables to the Dataset.

    Inputs
    ------
    data_xrds: xarray.Dataset.
        Timeseries data.
    varkeys, purtkeys: each a list of str's, same size.
        varkeys are the variable keys in data_xrds to get purturbations for. 
        purtkeys are the names to assign each purturbation variable.  

    Returns
    -------
    data_new: xarray.Dataset.
    """
    data_new = data_xrds.copy()
    for vk, pk in zip(varkeys, purtkeys):
        data_new[pk] = data_new[vk] - data_new[vk].mean()
    return data_new
    


def convert_time(dt64):
    """
    Convert a datetime value to seconds since midnight Jan 01, 2020.
    
    Inputs:
    -------
    dt64: numpy.datetime64 object.
        
    Returns:
    --------
    dt_secs2020jan01: float.
    """
    jan01_2020 = np.datetime64("2020-01-01")
    return (dt64 - jan01_2020).astype(float)/10**9
    







# Load data (currently for the 1/31 flight specifically):
fname_wind = "P3_FST_50Hz_20200131_135707.nc" # wind and temp.
fname_pic_mr = ("EUREC4A_ATOMIC_P3_Isotope-Analyzer_"
                "Water-Vapor-5Hz_20200131_v1.1.nc") # Picarro h2o.
fname_pic_iso = ("EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-Isotope-"
                 "Ratios-5Hz_20200131_v1.0.nc") # Picarro isotope ratios.
# Time intervals for two of the horizontal level legs:
tleg1 = (2667600, 2667913)
tleg2 = (2671900, 2672633)
tleg3 = (2670758, 2671573)
tleg4 = (2664350, 2665134)
tleg5 = (2663208, 2664007)
tleg6 = (2662429, 2662987)




tleg_list = (tleg1, tleg2, tleg3, tleg4, tleg5, tleg6)
levleg_str = ['levleg'+str(i) for i in range(1,7)]

for tleg, lls in zip(tleg_list, levleg_str):
    # 50 Hz files:
    data_50hz = process_50hz(tleg[0], tleg[1], fname_wind)
    data_50hz.to_netcdf('P-3_20200131_%s_50Hz.nc' % lls)
    # 5 hz files:
    data_5hz = process_5hz(
        tleg[0], tleg[1], 
        fname_wind, fname_pic_mr, fname_pic_iso
        )    
    data_5hz.to_netcdf('P-3_20200131_%s_5Hz.nc' % lls)








