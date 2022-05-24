# -*- coding: utf-8 -*-
"""
Created on Tue May 24 08:12:19 2022

@author: Dean

Estimate cloud tops for each of the P-3 cloud modules. 
"""



import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm

# Local code
import thermo



def cldtop_single(p3data):
    """
    Returns estimation of cloud top height for passed xarray.dataset with 
    cloud top key "alt_CT". Cloud top height is estimated as mean of data 
    greater than the 3rd quartile.
    """
    
    # P-3 data for highest level leg:    
    nlegtop = p3data['levlegflag'].max() # Highest P-3 horizontal leg
    p3datatop = p3data.where(p3data['levlegflag']==nlegtop, drop=True)
    #if len(p3datatop['time']<(5*8)): # Use next lowest leg if top leg is too short.
    #    p3datatop = p3data.where(p3data['levlegflag']==(nlegtop-1), drop=True)

    # Mean cloud top:
    q3 = np.quantile(p3datatop["alt_CT"].dropna(dim='time'), 0.75)
    ct_overq3 = p3datatop["alt_CT"] >= q3
    
    return p3datatop["alt_CT"].where(ct_overq3).mean().values.item()



ctops = []
tab_cld = pd.read_csv("./cldmod_keyaltitudes.csv")
dir_p3clddata = "./cldmod_datafiles/"
fnames_insitueremote = [f for f in os.listdir(dir_p3clddata)
                        if "_insitu+remote_" in f
                        ]
for i, row in tab_cld.iterrows():
    
    ncld = str(int(row['ncld'])).zfill(2)

    # Find P-3 insitue + remote datafile name (should only be one) and load:
    f = [f for f in fnames_insitueremote if "_ncld%s" % ncld in f]
    if len(f)==0: continue
    #fname_insituremote = "p3cld_insitu+remote_%i_ncld%s.nc" % tuple([date, ncld])
    p3data = xr.load_dataset(dir_p3clddata + f[0])
    
    ctops.append(cldtop_single(p3data))




