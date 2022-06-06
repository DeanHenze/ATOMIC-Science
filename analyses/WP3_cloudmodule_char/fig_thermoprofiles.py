# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:03 2022

@author: Dean
"""


# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import rangescaler
import profileplotter



## Data directories and filenames.
##_____________________________________________________________________________    
path_keyaltstable = "../WP3_cloudmodule_char/cldmod_keyaltitudes.csv"

path_clddrpsnds_dir = "./cldmod_datafiles/"
path_cldinsituremote_dir = "./cldmod_datafiles/"
fnames_clddrpsnds = [f for f in os.listdir(path_clddrpsnds_dir) 
                     if "p3cld_dropsondes_" in f]
fnames_cldinsituremote = [f for f in os.listdir(path_cldinsituremote_dir) 
                     if "p3cld_insitu+remote_" in f]
##_____________________________________________________________________________
## Data directories and filenames.



# Cloud module number groupings:
ncld_g1 = [7, 9, 11, 10, 12, 6]
ncld_g2 = [15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4, 8]

ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()



def drpsnd_meanprfs(ncld_list, fpaths, keyalts_table, scale_altkeys):
    """
    """

    # Averaged sounding for each cloud module, with scaled altitude:
    drpsndmeans_list = []
    for ncld, f in zip(ncld_list, fpaths):
        print(f[-14:-3])
        drpsnds = xr.load_dataset(f)
        drpsndmean = drpsnds.mean(dim='sounding')
    
        # Scale altitude and append to dataset:
        keyalts = keyalts_table.loc[keyalts_table['ncld']==ncld]
        alt_scalepoints = [0] + [keyalts[k].item() for k in scale_altkeys]
        altscaled = rangescaler.piecewise_linscale(
            drpsndmean['alt'].values, 
            alt_scalepoints, np.arange(len(alt_scalepoints))
            )
        drpsndmean = drpsndmean.assign_coords(
            alt_scaled = xr.DataArray(
                data=altscaled,
                dims=["alt"],
                coords=dict(alt=drpsndmean['alt']),
                )
            )
        #drpsndmean = drpsndmean.swap_dims({'alt':'alt_scaled'})

        #drpsndmean = drpsndmean.assign(
        #    alt_scaled = xr.DataArray(
        #        data=altscaled,
        #        dims=["alt"],
        #        coords=dict(alt=drpsndmean['alt']),
        #        )
        #    )
    
        
        drpsndmeans_list.append(drpsndmean)
        
        
    # Collect into an xr dataset:    
    ncld_sorted = ncld_list.copy()
    ncld_sorted.sort()
    ncld_idx = pd.Index(ncld_sorted, name='ncld')
    drpsndmeans_xr = xr.concat(drpsndmeans_list, ncld_idx, data_vars='all', 
                      coords='all', compat='identical', join='outer', 
                      combine_attrs='override')
            
    return drpsndmeans_xr
    


# Get dropsondes grouped by cloud characteristics. Each group is a dataset 
# and all datasets are collected in a list:
keyalts_table = pd.read_csv(path_keyaltstable)
scale_altkeys = ['z_lcl', 'z_tib']    
drpsnds_grouped = [] # Will be list of xr.Datasets:
for ncld_list in [ncld_g2, ncld_g3, ncld_g1]:
    fnames_cldgroup = [f for f in fnames_clddrpsnds if int(f[-5:-3]) in ncld_list]
    fpaths = [path_clddrpsnds_dir+f for f in fnames_cldgroup]
    drpsnds_grouped.append(
        drpsnd_meanprfs(
            ncld_list, fpaths, 
            keyalts_table, scale_altkeys
            )
        )
    


## Plot:
fig = plt.figure(figsize=(6.5, 2.5))
w = 0.215
lmarge = 0.1
axset = [
    fig.add_axes([x, 0.2, w, 0.78]) 
    for x in [lmarge, lmarge + w + 0.01, lmarge + 2*w + 2*0.01, lmarge + 3*w + 3*0.01]
    ]

pltcolors=['blue','red','grey']
for drpsnds, c in zip(drpsnds_grouped, pltcolors):

    for ax, varkey in zip(axset, ['theta','q','rh']):

        profileplotter.plotprf_singlevar(
            drpsnds[varkey].to_pandas().T, 
            ax, 
            alt_binwidth=100, pcolor=c
            #alt_binwidth=0.25, pcolor=c
            )
        























