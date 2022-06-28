# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:02:35 2022

@author: Dean

Figure comparing mean dD profiles to flux_dD profiles for each P-3 cloud group.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Local code
import profileplotter




## Data directories and filenames.
##_____________________________________________________________________________    
path_keyaltstable = "../WP3_cloudmodule_char/cldmod_keyaltitudes.csv"

path_cldinsituremote_dir = "../WP3_cloudmodule_char/cldmod_datafiles/"
fnames_cldinsituremote = [f for f in os.listdir(path_cldinsituremote_dir) 
                     if "p3cld_insitu+remote_" in f]
path_fluxdir = "../WP3_fluxes/fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(path_fluxdir) 
               if f.endswith(".nc")]
##_____________________________________________________________________________
## Data directories and filenames.



# Cloud module number groupings:
ncld_g1 = [1, 5, 4, 8]
ncld_g2 = [7, 9, 11, 10, 12, 6]
ncld_g3 = [15, 3, 2, 13, 16, 14]


ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()



## For scaling altitude by e.g. LCL, trade-inversion, if needed:
keyalts_table = pd.read_csv(path_keyaltstable)
scale_altkeys = ['z_lcl', 'z_tib']   



def insitu_meanprfs(ncld_list, fpaths, keyalts_table, scale_altkeys):
    """
    Returns mean profiles derived from insitu for a subset of the cloud P-3
    modules. Returned as an xarray Dataset.
    """

    # Averaged profile for each cloud module and append to a list:
    insitumeans_list = []
    for ncld, f in zip(ncld_list, fpaths):
        print(f[-14:-3])
        insitu = xr.load_dataset(f)
        
        # Get mean profile and append to list:
        insitumeans_list.append(
            insitu.groupby(100*np.round(insitu['alt']/100)).mean())
        
    # Collect into an xr dataset:    
    ncld_sorted = ncld_list.copy()
    ncld_sorted.sort()
    ncld_idx = pd.Index(ncld_sorted, name='ncld')
    insitumeans_xr = xr.concat(
        insitumeans_list, ncld_idx, data_vars='all', 
        coords='all', compat='identical', join='outer', 
        combine_attrs='override'
        )
                             
    return insitumeans_xr




def flux_scatter(ncld_list, path_fluxdir, fnames_flux, ax, **pltkwargs):
    """
    """
    for n in ncld_list:
        f = [f for f in fnames_flux if "_cld%s" % str(n).zfill(2) in f]
        f = f[0]
        data = xr.load_dataset(path_fluxdir + f)
        ax.scatter(data['dD_flux'], data['alt'], **pltkwargs)



## Get insitu mean profiles as a list of xr.Datasets, one for each cloud group: 
insitu_grouped = [] # Will be list of xr.Datasets:
for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
    fnames_cldgroup = [f for f in fnames_cldinsituremote if int(f[-5:-3]) in ncld_list]
    fpaths = [path_cldinsituremote_dir+f for f in fnames_cldgroup]
    insitu_grouped.append(
        insitu_meanprfs(
            ncld_list, fpaths, 
            keyalts_table, scale_altkeys
            )
        )
    
    

fig, axset = plt.subplots(1, 3, figsize=(6.5, 4))
#plt.figure()
#ax = plt.axes()
pltcolors=['red','grey','blue']
for insitu, c, ax in zip(insitu_grouped, pltcolors, axset):

        profileplotter.plotprf_singlevar(
            insitu['dD'].to_pandas().T, 
            ax, 
            altbinwidth=100, pcolor=c
            #alt_binwidth=0.25, pcolor=c
            )
        
        

for ax, nclds, c in zip(axset, [ncld_g1, ncld_g2, ncld_g3], pltcolors):
    flux_scatter(nclds, path_fluxdir, fnames_flux, ax, color=c, s=10)
    
    
for ax in axset:
    ax.set_ylim(-100, 3300)
    
yticks = np.arange(0, 3300, 500)
for ax in axset[1:]:
    ax.set_yticks(yticks)
    ax.set_yticklabels(['' for t in yticks])
    

axlabels = ['cg1', 'cg2', 'cg3']
for ax, lab in zip(axset, axlabels):
    ax.text(
        0.5, 1.01, lab, fontsize=12, 
        ha='center', va='bottom', transform=ax.transAxes
        )
    
    
axset[1].set_xlabel(r'$\delta D$'+u'(\u2030)', fontsize=12)
axset[0].set_ylabel('altitude (m)', fontsize=12)



legend_elements = [
    Line2D([0], [0], color='grey', lw=4, label='bulk'),
    Line2D([0], [0], marker='o', color='grey', label='flux',
           markerfacecolor='grey', linestyle='none', markersize=8)
    ]

# Create the figure
axset[1].legend(handles=legend_elements, loc='lower left', fontsize=9)


fig.savefig("./fig_dDbulkvsflux.png")












