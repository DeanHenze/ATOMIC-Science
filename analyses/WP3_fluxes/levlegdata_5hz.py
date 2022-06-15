# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:37:07 2022

@author: Dean


Isolate and process data for a P-3 horizontal level leg, prepping it for eddy 
covariance and spectral analyses at 5 Hz sampling rate. The following variables 
are collected into the file:
    latitude, longitude, altitude
    horizontal wind componenets
    vertical wind
    temperature
    water mass mixing ratio
    HDO mixing ratio
    
    
Current status: 
    The roll can be ~0 for the first few seconds of the time series and then 
    jump to > 5 after. Ideally, an algorithm will remove the first few 
    seconds along with the high roll interval.
"""



# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code:
import xcorr
import iso



def process_5hz(wind_50hz, mr_5hz, iso_5hz, roll_1hz, 
                t1, t2, timesync_method='xcorr', timesync_results=False):
    """ 
    Returns processed wind, water, and isotope ratio data at 5Hz, collected 
    into a single xarray dataset. Includes total values as well as 
    perturbations.

    Inputs
    ------
    wind_50hz: xarray.Dataset.
        Wind / temperature at 50 Hz sampling frequency, with time as the 
        single dimension.
    
    mr_5hz, iso_5hz: xarray.Dataset's.
        Water and HDO mixing ratio at 5 Hz sampling frequency, with time 
        as the single dimension.
    
    roll_1hz: xarray.Dataset.
        Includes aircraft roll data. Has 1Hz time as a single dimension.
    
    t1, t2: floats.
        Start and end times for segment to isolate.
        
    timesync_method: str.
        Either 'xcorr' or 'linear', see timesync function below.
        
    timesync_results: bool.
        If True, returns maximum cross correlation value and associate time 
        shift from syncing the water / iso data to wind.
        
    Returns
    -------
    data_merged: xarray.Dataset.
        Combined and processed wind, water, and isotope data.
        
    xcormax, tshift (optional): scalars.
        Returns these if timesync_results=True.
    """
    
    ## Process wind data
    ##_________________________________________________________________________            
    print("Processing wind data.")

    windkeys_new = ["u", "v", "w", "T"]
    wind_pro = varsubset(
        wind_50hz, t1, t2, 
        ["lon", "lat", "alt", "eastward_wind", 
         "northward_wind", "vertical_wind", "total_air_temperature"
         ], 
        keys_new=["lon", "lat", "alt"] + windkeys_new
        )
    wind_pro_df = wind_pro.to_dataframe() # Convert to pd for faster comps.
    wind_pro_df['time'] = wind_pro_df.index.values
    wind_pro_df = avgtofreq(wind_pro_df, 5, reportna=False)
    wind_pro_df = perturbations(wind_pro_df, windkeys_new, [k+"'" for k in windkeys_new])
    wind_pro = wind_pro_df.to_xarray() # Back to xarray for rest of fxn.
    ##_________________________________________________________________________            
    ## Process wind data
        
    
    ## Process water mixing ratio and isotope data
    ##_________________________________________________________________________            
    print("Processing water and isotope data.")

    # Remove missing time values then merge datasets:
    mr_5hz = mr_5hz.where(~mr_5hz['time'].isnull(), drop=True)
    iso_5hz = iso_5hz.where(~iso_5hz['time'].isnull(), drop=True)
    mriso = mr_5hz.merge(iso_5hz, join='exact')
    #mriso = mr_leg.merge(iso_leg, join='exact')    

    # Compute HDO mixing ratio 
    mriso['qD'] = iso.qD_dD_convert('dD2qD', mriso['mmr'], mriso['dD'])

    # Process the time dimension:
    tnew = [convert_time(dt64) for dt64 in mriso.time.values]
    mriso_pro = mriso.assign_coords(time=tnew)  
    interp_opts = {'bounds_error':None, 'fill_value':np.nan}
    mriso_pro = mriso_pro.interp(time=wind_pro['time'], method='nearest')
    
    # Variable subset in the cloud leg interval:
    mriso_pro = varsubset(
        mriso_pro, t1, t2, 
        ["mmr","qD"], 
        keys_new=["q","qD"]
        )
    mriso_pro = perturbations(
        mriso_pro, ["q","qD"], 
        ["q'","qD'"]
        )    
    
    # Time shift to wind data:     
    #mriso_pro, xcormax, tshift = timesync(
    tsyncresults = timesync(
        mriso_pro, wind_pro, "q'", "w'", 
        fs=5, method=timesync_method, leadmax=0, lagmax=7
        )
    if timesync_results: 
        mriso_pro, xcormax, tshift = tsyncresults  
    else:
        mriso_pro = tsyncresults
    ##_________________________________________________________________________            
    ## Process water mixing ratio and isotope data
    
    
    ## Aircraft roll
    ##_________________________________________________________________________            
    print("Processing roll data.")

    # If a dataset was passed, isolate the roll variable:
    if type(roll_1hz) == xr.core.dataset.Dataset:
        roll_1hz = roll_1hz['roll']  
    
    tnew_roll = [convert_time(dt64) for dt64 in roll_1hz.time.values]
    roll_1hz = roll_1hz.assign_coords(time=tnew_roll) 
    roll_leg = roll_1hz.sel(time=slice(t1, t2))
    roll_leg = roll_leg.interp(time=wind_pro['time'], method='nearest', kwargs=interp_opts)
    ##_________________________________________________________________________            
    ## Aircraft roll QC
    
    
    ## Merge all datasets and return:
    data_merged = xr.merge([wind_pro, mriso_pro, roll_leg], join='exact')
        # Sometimes there are leading or lagging NANs in mriso from time alignment:
    data_merged = data_merged.dropna(dim='time', how='any')
    
    if timesync_results:
        return data_merged, xcormax, tshift
    else:
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
    #data_tintv = data_xrds # NEW LINE IN TESTING !!!!!!
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
    
    
    
def avgtofreq(data, frq, reportna=False):
    """
    Average a dataset to a specified frequency. Interpolate any NAN values.
    
    Inputs
    ------
    data: pandas.DataFrame or xarray.Dataset.
        DataFrame / Dataset with 'time' as a coordinate or variable key.
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
    t_frq = np.round(data['time']/dt_frq)*dt_frq
    # Average duplicate time stamps to get desired result:
    data_avg = data.groupby(t_frq).mean()
        
    return data_avg



def perturbations(data, varkeys, purtkeys):
    """
    Computes purturbations from the mean for time series data. The mean is 
    taken over the entire time series. Purturbations are added as new 
    variables to the Dataset.

    Inputs
    ------
    data: pandas.DataFrame or xarray.Dataset.
        Timeseries data.
    varkeys, purtkeys: each a list of str's, same size.
        varkeys are the variable keys in data_xrds to get purturbations for. 
        purtkeys are the names to assign each purturbation variable.  

    Returns
    -------
    data_new: xarray.Dataset.
    """
    data_new = data.copy()
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



def timesync(ds1, ds2, k1, k2, fs=5, 
             method='xcorr', leadmax=0, lagmax=7, altkey='alt'):
    """
    Time shift variables in ds1 to ds2 using either max cross-correlation or 
    predetermined linear function of altitude. ds1 time stamps are interpolated 
    to ds2 timestamps after shifting.
    
    Inputs
    ------
    ds1, ds2: xarray.Dataset's.
        Each should have time as a dimension and the same sampling freqency. 
        If time_sync method is 'linear', ds2 should have altitude as a 
        variable.
        
    k1, k2: strs.
        Keys for the variables in ds1 and ds2 to cross-correlate.
        
    fs: scalar.
        Sampling frequency in Hz for ds1 and ds2.
        
    method: 'xcorr' or 'linear'
        If 'xcorr' use maximum cross-correlation method. If 'linear', use a 
        predetermined linear function of altitude.
        
    leadmax, lagmax: scalars.
        Maximum lead and lag (seconds) time shift to consider if using the 
        'xcorr' method.
        
    altkey: 'str'
        Key for altitude in ds1, applicable if method='linear'.

    Returns
    -------
    ds1: xarray.Dataset.
        ds1 time shifted.
        
    xcormax, tshift: scalars
        Returned if method='xcorr'. Value of maximum cross-correlation and 
        corresponding time shift.
    """
    # Compute time shift:
    if method=='xcorr':
        # Determine time shift using max cross-correlation:
        c, xcor, n = xcorr.correlation(ds1[k1], ds2[k2], leadmax*fs, lagmax*fs)
        imaxcor = np.argmax(xcor)
        tshift = c[imaxcor]/fs
        xcormax = xcor[imaxcor]
    
    if method=='linear':
        # Determine using preset linear function of altitude:
        m = 0.0002419
        b = -3.838
        tshift = m*ds2[altkey].mean().item() + b
        
    # Apply shift and interpolate to ds2 time values:
    ds1 = ds1.assign_coords(time=ds1['time']+tshift)  
    ds1 = ds1.interp(
        time=ds2['time'], 
        method='nearest'
        )
    
    if method=='xcorr':
        return ds1, xcormax, tshift 
    else:
        return ds1
















