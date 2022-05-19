# -*- coding: utf-8 -*-
"""
Created on Thu May 19 12:21:25 2022

@author: Dean

Goal is to find the level leg closest to the estimated LCL.
"""


import os

import numpy as np
import pandas as pd
import xarray as xr


# Level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]


# Key altitude info:
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")



def lcl_legfinder(ncld, fnames_levlegs, keyalts_table):
    """
    """
    keyalts_ncld = keyalts_table.loc[keyalts_table['ncld']==ncld]
    fnames_levlegs_ncld = [f for f in fnames_levlegs if "_cld%i_"%ncld in f]

    z_lcl = keyalts_ncld['z_lcl'].values.item()
    z_legs = []
    for f in fnames_levlegs_ncld:
    
        data = xr.load_dataset(dir_5hzdata + f)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        z_legs.append(data_df['alt'].mean())
        
    z_diff = np.array(z_legs) - z_lcl
    #return z_diff[(z_diff < 0) & (z_diff > -200)]
    return z_diff[(z_diff < 0) | (abs(z_diff) < 600)]


z_diff = []
for n in np.arange(1,17,1): 
    z_diff.append(lcl_legfinder(n, fnames_levlegs, keyalts_table))
