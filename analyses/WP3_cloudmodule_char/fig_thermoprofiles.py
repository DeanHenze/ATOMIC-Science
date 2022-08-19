# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:03 2022

@author: Dean

Figure of mean profiles for potential temperature, humidity, relative humidity 
and water isotope ratio dD for P-3 cloud modules. Cloud modules are 
grouped by cloud properties.

Use in-situ measurements for dD profiles and dropsonde data for the other 
three variables.
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
import rangescaler
import profileplotter
import thermo



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
ncld_g1 = [1, 5, 4]
ncld_g2 = [8, 7, 9, 11, 10, 6]
ncld_g3 = [12, 15, 3, 2, 13, 16, 14]

ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()



def drpsnd_meanprfs(ncld_list, fpaths, keyalts_table, scale_altkeys):
    """
    Returns mean profiles derived from dropsondes for a subset of the P-3 cloud 
    modules. Returned as an xarray Dataset.
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
        
        drpsndmeans_list.append(drpsndmean)
        
        
    # Collect into an xr dataset:    
    ncld_sorted = ncld_list.copy()
    ncld_sorted.sort()
    ncld_idx = pd.Index(ncld_sorted, name='ncld')
    drpsndmeans_xr = xr.concat(
        drpsndmeans_list, ncld_idx, data_vars='all', 
        coords='all', compat='identical', join='outer', 
        combine_attrs='override'
        )
            
    return drpsndmeans_xr



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
        
        # Append potential temperature:
        insitu['theta'] = thermo.theta(insitu['Ta'], insitu['press'], qv=0)
        
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
     


## For scaling altitude by e.g. LCL, trade-inversion, if needed:
keyalts_table = pd.read_csv(path_keyaltstable)
scale_altkeys = ['z_lcl', 'z_tib']   



## Get insitu mean profiles as a list of xr.Datasets, one for each cloud group: 
insitu_grouped = [] # Will be list of xr.Datasets:
for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
    fnames_cldgroup = [f for f in fnames_cldinsituremote if int(f[-5:-3]) in ncld_list]
    fpaths = [path_clddrpsnds_dir+f for f in fnames_cldgroup]
    insitu_grouped.append(
        insitu_meanprfs(
            ncld_list, fpaths, 
            keyalts_table, scale_altkeys
            )
        )
    
    

## Get dropsonde mean profiles as a list of xr.Datasets, one for each cloud group:  
drpsnds_grouped = [] # Will be list of xr.Datasets:
for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
    fnames_cldgroup = [f for f in fnames_clddrpsnds if int(f[-5:-3]) in ncld_list]
    fpaths = [path_clddrpsnds_dir+f for f in fnames_cldgroup]
    drpsnds_grouped.append(
        drpsnd_meanprfs(
            ncld_list, fpaths, 
            keyalts_table, scale_altkeys
            )
        )
    
    
## Save mean profiles:
def get_meanprfs(data_grouped, varkey, avg_over):
    """
    data_grouped: list (length=3) of xarray datasets.
    varkey: key of variable to get profiles for
    avg_over: key of dimension to take mean over in the datasets.
    """
    # Mean profile for each cloud group as a dataframe with altitude as index:
    ng_cloudgroups = [1, 2, 3]
    meanprfs_list = [] # Mean profiles go here.
    for data, ng in zip(data_grouped, ng_cloudgroups):
        meanprf_ng = data[varkey].mean(dim=avg_over).to_dataframe()
        meanprf_ng.rename(columns={varkey:'cg%i' % ng}, inplace=True)
        meanprfs_list.append(meanprf_ng)
    
    # Return profiles merged into a single dataframe:
    meanprfs = pd.merge(meanprfs_list[0], meanprfs_list[1], how='outer', 
                        left_index=True, right_index=True)
    meanprfs = pd.merge(meanprfs, meanprfs_list[2], how='outer', 
                        left_index=True, right_index=True)
    return meanprfs

path_savedir = "./mean_profiles/"
if not os.path.isdir(path_savedir): os.mkdir(path_savedir)
    # In-situe dD:  
dD_meanprf = get_meanprfs(insitu_grouped, 'dD', 'ncld')
dD_meanprf.to_csv(path_savedir + "meanprf_dD_WP3.csv")
    # Dropsonde T, theta, q, RH:
for varkey in ['ta', 'theta', 'q', 'rh']:
    meanprf = get_meanprfs(drpsnds_grouped, varkey, 'ncld')
    meanprf.to_csv(path_savedir + "meanprf_%s_WP3.csv" % varkey)
    


## Plot:
fig = plt.figure(figsize=(6.5, 3))
w = 0.21
lmarge = 0.12
axset = [
    fig.add_axes([x, 0.2, w, 0.78]) 
    for x in [lmarge, lmarge + w + 0.01, lmarge + 2*w + 2*0.01, lmarge + 3*w + 3*0.01]
    ]

pltcolors=['red','grey', 'blue']
for drpsnds, c in zip(drpsnds_grouped, pltcolors):

    for ax, varkey in zip(axset, ['theta','q','rh']):

        profileplotter.plotprf_singlevar(
            drpsnds[varkey].to_pandas().T, 
            ax, 
            alt_binwidth=100, pcolor=c
            #alt_binwidth=0.25, pcolor=c
            )
        
        
for insitu, c in zip(insitu_grouped, pltcolors):

        profileplotter.plotprf_singlevar(
            insitu['dD'].to_pandas().T, 
            axset[3], 
            alt_binwidth=100, pcolor=c
            #alt_binwidth=0.25, pcolor=c
            )


## Axes limits, labels, legend, and save:
axset[0].set_xlim(294, 315)
axset[0].set_xticks(np.arange(295, 315, 5))
axset[0].set_xticklabels(axset[0].get_xticks().astype(str), fontsize=9)
axset[0].set_xlabel(r'$\theta$ (K)', fontsize=12)


axset[1].set_xlim(-0.0006, 0.02)
axset[1].set_xticks(np.arange(0, 0.02, 0.005))
axset[1].set_xticklabels(
    [(1000*t).astype(int).astype(str) for t in axset[1].get_xticks()],
    fontsize=9
    ) # factor of 1000 to convert to g/kg.
axset[1].set_xlabel(r'$q$ (g/kg)', fontsize=12)


axset[2].set_xlim(-0.02, 1.05)
axset[2].set_xticks(np.arange(0, 1.01, 0.25))
axset[2].set_xticklabels(
    [(100*t).astype(int).astype(str) for t in axset[2].get_xticks()],
    fontsize=9
    ) # factor of 100 to convert to %.
axset[2].set_xlabel(r'$RH$ (%)', fontsize=12)


axset[3].set_xlim(-285, -40)
axset[3].set_xticks(np.arange(-250, -49, 50))
axset[3].set_xticklabels(axset[3].get_xticks().astype(str), fontsize=9)
axset[3].set_xlabel(r'$\delta D$ '+u'(\u2030)', fontsize=12)


for ax in axset: 
    ax.set_ylim(-100, 3200)
    ax.set_yticks(np.arange(0, 3001, 500))   
axset[0].set_yticklabels(axset[0].get_yticks().astype(str), fontsize=9)        
for ax in axset[1:]: 
    ax.set_yticklabels(['' for e in ax.get_yticks()], fontsize=9)
axset[0].set_ylabel('altitude (m)', fontsize=12)


legend_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='blue', lw=4)]
ax.legend(legend_lines, ['LCLT', 'LCHT', 'HCHT'])


fig.savefig("./fig_thermoprofiles.png")