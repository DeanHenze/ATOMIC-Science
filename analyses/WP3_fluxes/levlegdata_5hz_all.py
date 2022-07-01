# -*- coding: utf-8 -*-
"""
Created on Tue May  3 10:41:10 2022

@author: Dean

Processes all P-3 level leg fast data at 5 Hz. Creates a sepearate file for 
each leg. 

The following warning may be raised: 
    
"RuntimeWarning: invalid value encountered in multiply
    flat_num_dates.astype(np.float64) * _NS_PER_TIME_DELTA[delta]"
    
It occurs because there are NANs in the time variable for some of the Picarro 
files.
"""


import os

import pandas as pd
import xarray as xr

import levlegdata_5hz



# Load level leg time interval data and onvert start / end times to timestamps:
time_levlegs = pd.read_csv("p3_cloudmodule_legs.csv").copy()
time_levlegs['tstart_tstamp'] = time_levlegs['tstart_leg'].apply(pd.Timestamp)
time_levlegs['tend_tstamp'] = time_levlegs['tend_leg'].apply(pd.Timestamp)


# Paths to data sets:
path_fastwind = "../../data/WP3/fast_wind_temp/" # High freq wind / temp
path_fastwater = "../../data/WP3/fast_water_iso/" # high freq water / iso
path_flightlev1hz = "../../data/WP3/flight_level_1Hz/" # files with roll and temperature data


# Save results to these directories:
    # Processed 5Hz files:
dir_5hzfiles = "./levlegdata_5hz/"
if not os.path.isdir(dir_5hzfiles): os.makedirs(dir_5hzfiles)
    
    # Time sync results:
dir_tsync = "./other_data/"
if not os.path.isdir(dir_tsync): os.makedirs(dir_tsync)


# Store data for time synce results .csv file in these lists:
xcormax_list = []
tshift_list = []
altleg_list = [] # mean altitude of the level leg.
ncld = []
nleg = []


for i, row in time_levlegs.iterrows():
        
    date = str(row['flight_date'])
    
    print("Flight on, %s, cloud module %i, leg number %i"
          % tuple([date, row['num_cld_iop'], row['num_leg']])
          )
    
    # Load wind and temperature data:
    fnames_fastwind = os.listdir(path_fastwind)
    f = [f for f in fnames_fastwind     # Find filename for current date.
         if f.startswith('P3_FST_50Hz_%s' % date)
         ]
    if len(f)==0: continue # A couple of dates don't have files.
    fname_wind = path_fastwind + f[0] # There should be only one.
    wind = xr.load_dataset(fname_wind)
    
    
    # Load water mixing ratio (mr) and isotope data from the Picarro:
    fnamehead_mr = "EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-5Hz_"
    fname_mr = (path_fastwater + fnamehead_mr + date + "_v1.1.nc")
    mr = xr.load_dataset(fname_mr)
    
    fnamehead_iso = "EUREC4A_ATOMIC_P3_Isotope-Analyzer_Water-Vapor-Isotope-Ratios-5Hz_"
    fname_iso = (path_fastwater + fnamehead_iso + date + "_v1.0.nc")
    iso = xr.load_dataset(fname_iso)
    
    
    # Load 1 Hz flight level data, get roll and temperature:
    fnamehead_flev = "EUREC4A_ATOMIC_P3_Flight-level_"
    fname_flev = (path_flightlev1hz + fnamehead_flev + date + "_v1.0.nc")
    flightlev = xr.load_dataset(fname_flev)    
    roll = flightlev['roll']
    T = flightlev['Ta'].to_dataset()
    
    
    # Get processed data and water-wind time sync results:
    t1 = levlegdata_5hz.convert_time(row['tstart_tstamp'].to_datetime64())
    t2 = levlegdata_5hz.convert_time(row['tend_tstamp'].to_datetime64())
        # If using the original cross correlation time shift:
    #data_proc, xcormax, tshift = levlegdata_5hz.process_5hz(
    #    wind, mr, iso, roll, 
    #    t1, t2, timesync_results=True
    #    )
        # With new linear function of altitude:
    data_proc = levlegdata_5hz.process_5hz(
        wind, mr, iso, roll, T, 
        t1, t2, timesync_method='linear', timesync_results=False
        )
    
    
    # Save processed data as new file:
    ncld_str = str(int(row['num_cld_iop'])).zfill(2)
    data_proc.to_netcdf(dir_5hzfiles + "WP3_5hz_%s_cld%s_levleg%i.nc"
                        % tuple([date, ncld_str, row['num_leg']])
                        )
    
    
    # Append time sync results:
    #xcormax_list.append(xcormax)
    #tshift_list.append(tshift)
    altleg_list.append(data_proc['alt'].mean().values.item())
    ncld.append(row['num_cld_iop'])
    nleg.append(row['num_leg'])
    
    
    del wind, mr, iso, roll, data_proc
    
"""  
# Save time sync results:
tsync_results = pd.DataFrame({
    'tshift': tshift_list, 'xcormax': xcormax_list,
    'alt_leg': altleg_list, 'ncld': ncld, 'nleg': nleg 
    })
tsync_results.to_csv(dir_tsync + "timesync_wind-water_results.csv", index=False)
"""