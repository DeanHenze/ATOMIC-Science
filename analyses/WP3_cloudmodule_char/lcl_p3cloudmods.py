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



path_drpsnd = "./cldmod_datafiles/" # Path to cloud mod. dropsondes folder.
fnames = [f for f in os.listdir(path_drpsnd) 
          if f.startswith("p3cld_dropsondes") and f.endswith(".nc")
          ]

keyalts_tab = pd.read_csv("./cldmod_keyaltitudes.csv")

Tlcl = [] # LCL temperatures appended here.
ncld_list = np.arange(1,17,1) # 16 Cloud modules over the IOP
for n in ncld_list:
    
    f = [f for f in fnames if "_ncld%s"%str(n).zfill(2) in f][0]
    drpsnds = xr.load_dataset(path_drpsnd + f)
    
    drpsnds_ns = drpsnds.sel(alt=slice(30,100)) # near surface portions.
    sfcmean = drpsnds_ns.mean()[['ta', 'p', 'q']]
    
    Tlcl.append(thermo.Tlcl(sfcmean['ta'], sfcmean['p'], sfcmean['q']).values)





















