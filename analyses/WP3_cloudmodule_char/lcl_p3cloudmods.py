# -*- coding: utf-8 -*-
"""
Created on Tue May 17 13:59:41 2022

@author: Dean

Get LCL for each P-3 cloud module.

LCL calcs will use near-surface data from the dropsondes released during each 
cloud module.
"""


# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr

# Local code
import thermo


gamma = 9.8/1000 # dry adiabatic lapse rate (K/m)


# Cloud module dropsondes directory and filenames:
path_drpsnd = "./cldmod_datafiles/"
fnames = [f for f in os.listdir(path_drpsnd) 
          if f.startswith("p3cld_dropsondes") and f.endswith(".nc")
          ]


Tlcl_cldmods = [] # LCL temperatures appended here.
zlcl_cldmods = [] # LCL heights appended here.
ncld = np.arange(1,17,1) # 16 Cloud modules over the IOP
for n in ncld:
    
    # Load dropsondes for this cloud mod:
    f = [f for f in fnames if "_ncld%s"%str(n).zfill(2) in f][0]
    drpsnds = xr.load_dataset(path_drpsnd + f)
    
    # Get near surface mean thermo quantities:
    drpsnds_ns = drpsnds.sel(alt=slice(20,100)) # near surface portions.
    sfcmean = drpsnds_ns.mean()[['ta', 'p', 'q']]
    
    # Compute LCL temperature and height:
    Tlcl = (thermo.Tlcl(sfcmean['ta'], sfcmean['p'], sfcmean['q']).values)
    Tlcl_cldmods.append(Tlcl.item())
    zlcl = (sfcmean['ta'] - Tlcl)/gamma
    zlcl_cldmods.append(zlcl.values.item())


# Append LCL height as column in the key altitudes table:
path_keyaltstab = "./cldmod_keyaltitudes.csv"
keyalts_tab = pd.read_csv(path_keyaltstab)
if ((keyalts_tab['ncld'] == ncld).sum()) == len(ncld): # check that cols match
    keyalts_tab['alt_LCL'] = zlcl_cldmods
    keyalts_tab.to_csv(path_keyaltstab, index=False)
else:
    keyalts_tab.sort_values('ncld', axis=0, ascending=True, inplace=True)
    keyalts_tab['alt_LCL'] = zlcl_cldmods
    keyalts_tab.to_csv(path_keyaltstab, index=False)
















