# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:37:26 2022

@author: Dean

Mean altitude for each P-3 cloud-module level-leg is computed and the results 
are saved as a new .csv file.
"""



# Built in
import os

# Third party
import pandas as pd
import xarray as xr



if __name__=="__main__":
    
    # Level leg data filenames:
    dir_5hzdata = "./levlegdata_5hz/"
    fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]
    fnames_levlegs.sort()
    
    
    # Get cloud module number, level leg number, and mean altitude in lists:
    altmean_list = []
    ncld_list = []
    nlevleg_list = []
    for f in fnames_levlegs:
        
        # Get cloud module and level leg numbers from filename:
        icldmod = f.index("_cld")
        ncld = int(f[icldmod+4: icldmod+6])
        ilevleg = f.index("_levleg")
        nlevleg = int(f[ilevleg+7: ilevleg+8])
        ncld_list.append(ncld)
        nlevleg_list.append(nlevleg)
        
        # Load P-3 insitu data and get mean altitude:
        data = xr.load_dataset(dir_5hzdata + f)
        altmean_list.append(data['alt'].mean().round().item())
                
     
    # Save lists as columns in a .csv file:
    results = pd.DataFrame({
        'ncld':ncld_list, 'nlevleg':nlevleg_list, 
        'altmean':altmean_list
        })
    results.to_csv("p3cld_levleg_altitudes.csv", index=False)