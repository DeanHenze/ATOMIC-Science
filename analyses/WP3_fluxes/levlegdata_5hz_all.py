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



# Load level leg time interval data:
time_levlegs = pd.read_csv("p3_cloudmodule_legs.csv").copy()

#
time_levlegs['tstart_tstamp'] = time_levlegs['tstart_leg'].apply(pd.Timestamp)
time_levlegs['tend_tstamp'] = time_levlegs['tend_leg'].apply(pd.Timestamp)


# Paths to data sets:
path_fastwind = "../../data/WP3/fast_wind_temp/" # High freq wind / temp
path_fastwater = "../../data/WP3/fast_water_iso/" # high freq water / iso
path_flightlev1hz = "../../data/WP3/flight_level_1Hz/" # files with roll data


#row = time_levlegs.iloc[0]
#fname = "EUREC4A_ATOMIC_P3_Flight-level_20200124_v1.0.nc"
#roll = xr.load_dataset(path_flightlev1hz + fname)
#roll = roll[['lat','lon','alt','roll']]



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
    
    
    # Load 1 Hz flight level data and get roll var:
    fnamehead_flev = "EUREC4A_ATOMIC_P3_Flight-level_"
    fname_flev = (path_flightlev1hz + fnamehead_flev + date + "_v1.0.nc")
    flightlev = xr.load_dataset(fname_flev)    
    roll = flightlev['roll']
    
    
    # Get processed data and save as new file:
    t1 = levlegdata_5hz.convert_time(row['tstart_tstamp'].to_datetime64())
    t2 = levlegdata_5hz.convert_time(row['tend_tstamp'].to_datetime64())
    data_proc = levlegdata_5hz.process_5hz(wind, mr, iso, roll, t1, t2)
    data_proc.to_netcdf("./levlegdata_5hz/WP3_5hz_%s_cld%i_levleg%i.nc"
                        % tuple([date, row['num_cld_iop'], row['num_leg']])
                        )





"""
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
    print(lls)
    # 5 hz files:
    data_5hz = levlegdata_5hz.process_5hz(
        wind, mr, iso,
        tleg[0], tleg[1]
        )    
    data_5hz.to_netcdf('./levlegdata_5hz/P-3_20200131_%s_5Hz.nc' % lls)
"""