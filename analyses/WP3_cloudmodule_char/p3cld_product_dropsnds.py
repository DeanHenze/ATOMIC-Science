# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:53:38 2022

@author: Dean

Dropsonde data product for the P-3 cloud modules. Saves one file of 
dropsonde data per cloud module.
"""



# Built in:
import sys

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# My modules:
if r'../../' not in sys.path: sys.path.insert(0,r'../../')
from henze_python_modules import atomic_data_loader as adl



def dropsnds_cld_single(date, ncld, t1_cld, t2_cld):
    """
    date: int
    ncld: int, float
    t1_cld, t2_cld: Timestamps
    tab_hlegs: pandas.DataFrame
    """
    dropsnds = adl.dropsondes(version='Henze', level=3) # Load data.
    dropsnds_cld = dropsnds.where( # P-3 sondes in cld time interval.
        (dropsnds['launch_time']>t1_cld) 
        & (dropsnds['launch_time']<t2_cld)
        & (dropsnds['platform']=='P3'),
        drop=True
        )
    return dropsnds_cld



###############################################################################
# Remainder of script creates a separate data file for each cloud module.
###############################################################################

# Load cloud module info tables:
tab_cld = adl.p3_cld_modules(datetime_type=True) # Cloud modules

for i, row in tab_cld.iterrows():
    print(row['flight_date'])
    
    date = row['flight_date']
    ncld = row['num_cld']
    
    dropsnds_cld = dropsnds_cld_single(
        date, ncld, 
        row['start_datetime'], row['end_datetime'], 
        )
    
    
    # Save:
    ncld_int = str(int(ncld)).zfill(2)
    fname = r"../output/p3cld_dropsondes_%i_%s.nc" % tuple([date, ncld_int])
    dropsnds_cld.to_netcdf(fname)




