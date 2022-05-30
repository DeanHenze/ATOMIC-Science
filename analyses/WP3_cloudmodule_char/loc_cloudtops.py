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
    Returns estimation of cloud top height (meters) for passed xarray.dataset 
    with cloud top key "alt_CT". Cloud top height is estimated as mean of data 
    greater than the 3rd quartile.
    """
    
    # P-3 data for highest level leg:    
    nlegtop = p3data['levlegflag'].max() # Highest P-3 horizontal leg
    p3datatop = p3data.where(p3data['levlegflag']==nlegtop, drop=True)
    #if len(p3datatop['time']<(5*8)): # Use next lowest leg if top leg is too short.
    #    p3datatop = p3data.where(p3data['levlegflag']==(nlegtop-1), drop=True)

    # Mean cloud top:
    q_50p = np.quantile(p3datatop["alt_CT"].dropna(dim='time'), 0.50)
    q_95p = np.quantile(p3datatop["alt_CT"].dropna(dim='time'), 0.95)
    ct_over50p = p3datatop["alt_CT"] >= q_50p
    ct_under95p = p3datatop["alt_CT"] <= q_95p
    
    ctmean = p3datatop["alt_CT"].where(ct_over50p & ct_under95p).mean().item()
    return q_50p*1000, q_95p*1000, ctmean*1000 # convert all values to km.



if __name__=="__main__":
    """
    Cloud top heights computed and appended to key altitudes table.
    
    Three quantites computed: cloud top 50th percentile, 95th percentile, 
    and mean of data in this percentile range.
    """
    
    # Cloud module key altitudes table:
    path_cldtab = "./cldmod_keyaltitudes.csv"
    tab_cld = pd.read_csv(path_cldtab)
    
    # P-3 filenames for insitu+remote:
    dir_p3clddata = "./cldmod_datafiles/"
    fnames_insitueremote = [f for f in os.listdir(dir_p3clddata)
                            if "_insitu+remote_" in f
                            ]
    
    # Get cloud top height percentiles and means:
    ct_50p = []
    ct_95p = []
    ctmean = [] # Mean of the 50th - 95th percentile of cloud tops.
    for i, row in tab_cld.iterrows():
        
        ncld = str(int(row['ncld'])).zfill(2)
        print("Working on cloud module %s" % ncld)
    
        # Find P-3 insitue + remote datafile name (should only be one) and load:
        f = [f for f in fnames_insitueremote if "_ncld%s" % ncld in f]
        if len(f)==0: continue
        p3data = xr.load_dataset(dir_p3clddata + f[0])
        
        ctresults = cldtop_single(p3data)
        ct_50p.append(ctresults[0])
        ct_95p.append(ctresults[1])
        ctmean.append(ctresults[2])
        
    # Append as new columns and save:
    tab_cld['z_ct50p'] = np.round(ct_50p)
    tab_cld['z_ct95p'] = np.round(ct_95p)
    tab_cld['z_ctmean_50p95p'] = np.round(ctmean)
    tab_cld.to_csv(path_cldtab)




