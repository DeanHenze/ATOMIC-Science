# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:03 2022

@author: Dean
"""


# Built in
import sys
import os
from os import listdir
from os.path import join

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

# Local code
import rangescaler



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
        drpsndmean = drpsndmean.assign(
            alt_scaled = xr.DataArray(
                data=altscaled,
                dims=["alt"],
                coords=dict(alt=drpsndmean['alt']),
                )
            )
    
        
        drpsndmeans_list.append(drpsndmean)
        
        
    # Collect into an xr dataset:    
    ncld_sorted = ncld_list.copy()
    ncld_sorted.sort()
    ncld_idx = pd.Index(ncld_sorted, name='ncld')
    drpsndmeans_xr = xr.concat(drpsndmeans_list, ncld_idx, data_vars='all', 
                      coords='all', compat='identical', join='outer', 
                      combine_attrs='override')
            
    return drpsndmeans_xr
      

for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
    keyalts_table = pd.read_csv(path_keyaltstable)
    scale_altkeys = ['z_lcl', 'z_tib']
    #scale_altkeys = ['z_lcl', 'z_ctmean_50p95p']
    #ncld_list = ncld_g2
    fnames_cldgroup = [f for f in fnames_clddrpsnds if int(f[-5:-3]) in ncld_list]
    fpaths = [path_clddrpsnds_dir+f for f in fnames_cldgroup]
    drpsnds = drpsnd_meanprfs(ncld_list, fpaths, keyalts_table, scale_altkeys)
    
    
    plt.figure()
    for n in ncld_list:
        test = drpsnds.sel(ncld=n)
        plt.plot(test.theta, test.alt_scaled, label=n)
    plt.legend()
   
""" 
plt.figure()
for n in ncld_list:
    test = drpsnds.sel(ncld=n)
    plt.plot(test.theta, test.alt, label=n)
plt.legend()
"""


def prf_cldgroup(ncldgroup):
    """
    Get mean and std for all dropsondes in a subset of the cloud modules.
    """
    
    fnames_cldgroup = [f for f in fnames_clddrpsnds if int(f[-5:-3]) in ncldgroup]
    drpsnds = drpsnd_meanprfs([path_clddrpsnds_dir+f for f in fnames_cldgroup])
    mean_drpsnd = drpsnds.mean(dim=['cloud_module', 'sounding'])
    std_drpsnd = drpsnds.std(dim=['cloud_module', 'sounding'])
    return mean_drpsnd, std_drpsnd
    
"""    
drpsndresults = []
for ncldgroup in [ncld_g1, ncld_g2, ncld_g3]:
        drpsndresults.append(prf_cldgroup(ncldgroup))

plt.figure()
ax = plt.axes() 
colors = ['grey', 'blue', 'red']   
for results, c in zip(drpsndresults, colors):
    
    mean = results[0]
    std = results[1]
    ax.plot(mean['theta'], mean['alt'], color=c)
    ax.fill_betweenx(
        mean['alt'], mean['theta']-std['theta'], x2=mean['theta']+std['theta'], 
        alpha=0.25, color=c, edgecolor='none'
        )
"""


def plot_prf_envelop(data_xr, axset, varset, color):
    
    means = data_xr.mean(dim='cloud_module')
    maxvals = data_xr.max(dim='cloud_module')
    minvals = data_xr.min(dim='cloud_module')

    for ax, v in zip(axset, varset):
        ax.plot(means[v], means['alt'], color=color)
        ax.fill_betweenx(
            means['alt'], minvals[v], x2=maxvals[v], 
            alpha=0.25, color=color, edgecolor='none'
            )
        
"""
    # Plot:
fig = plt.figure(figsize=(6.5, 2.5))
w = 0.215
lmarge = 0.1
axset = [
    fig.add_axes([x, 0.2, w, 0.78]) 
    for x in [lmarge, lmarge + w + 0.01, lmarge + 2*w + 2*0.01, lmarge + 3*w + 3*0.01]
    ]


plot_prf_envelop(p3other, axset, ['theta','mr','rh','dD'], 'blue')
plot_prf_envelop(p3nightconv, axset, ['theta','mr','rh','dD'], 'grey')
plot_prf_envelop(p3dayconv, axset, ['theta','mr','rh','dD'], 'red')


ax1 = axset[0]
ax1.set_ylim(0, 4000)
ax1.set_xlim(297, 320)
ax1.set_xlabel(r'$\theta$ (K)', fontsize=12)
ax1.set_ylabel(r'altitude (m)', fontsize=12)

ax2 = axset[1]
ax2.set_ylim(0, 4000)
ax2.set_xlim(0, 16)
ax2.set_xlabel(r'q (g/kg)', fontsize=12)

ax3 = axset[2]
ax3.set_ylim(0, 4000)
ax3.set_xlabel(r'RH (%)', fontsize=12)

ax4 = axset[3]
ax4.set_ylim(0, 4000)
ax4.set_xlabel(r'$\delta D$ (permil)', fontsize=12)


for ax in axset[1:]:
    ax.set_yticklabels(['' for t in axset[0].get_xticklabels()])
"""

#fig.savefig(r'../figures/day-night_prfs.png')























