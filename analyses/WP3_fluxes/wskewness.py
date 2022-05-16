# -*- coding: utf-8 -*-
"""
Created on Mon May 16 09:35:38 2022

@author: Dean

Compute vertical velocity skewness for each level leg.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal

# Local code
import rollqc
import thermo



# Save skewness results to this directory:
dir_wskew = "./wskew_levlegs/"
if not os.path.isdir(dir_wskew):
    os.makedirs(dir_wskew)
        

# All level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]
   
   
# Fxn to compute skewness:
def skewness(x): # x=array-like.
    xmean = np.nanmean(x)
    std = np.nanstd(x)
    return np.nanmean((x-xmean)**3)/std**3


varkeys = ["w'"]
wskew = []
ncld = []
nleg = []
altleg = []
for f in fnames_levlegs:
    
    data = xr.load_dataset(dir_5hzdata + f)
    data_df = data.to_dataframe() # Easier / faster to work with pandas.

    # Aircraft roll QC
    # Remove data where the roll was greater than 5 degrees:
    roll_crit = 5
    highroll = abs(data_df['roll']) > roll_crit
    data_dfqc = data_df.loc[~highroll]

    # Compute w skewness and add to list:
    wskew.append(skewness(data_dfqc["w'"]))
    
    # Add other data to lists:
    altleg.append(np.nanmean(data_dfqc['alt']))
    i_levleg = f.index("_levleg")
    ncld.append( int(f[f.index("_cld")+4: i_levleg]) )
    nleg.append( int(f[i_levleg+7]) )
    
    
# Collect info pandas df and save:
dfkeys = ['ncld', 'nleg', 'altleg', 'wskew']
wskew_df = pd.DataFrame(dict(zip(dfkeys, [ncld, nleg, altleg, wskew])))
wskew_df.to_csv(dir_wskew + "WP3_wskewness_levlegs.csv", index=False)