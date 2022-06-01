# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 10:04:03 2022

@author: Dean
"""


# Built in:
import sys
import os
from os import listdir
from os.path import join

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

# My modules:
if r'../../' not in sys.path: sys.path.insert(0,r'../../')
from henze_python_modules import atomic_data_loader as adl
from henze_python_modules import iso_fxns
from henze_python_modules import atms_physics as ap




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
ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4]


"""
    # Filenames for all P-3 cloud module in-situ data:
datadir = r"../output/"
fnames_insitremote = [
    join(datadir, f) for f in listdir(datadir) 
    if f.startswith("p3cld_insitu+remote") and f.endswith(".nc")
    ]
dates_dayconv = ['20200124']
dates_nightconv = ['20200209','20200210','20200211']
fnames_dayconv = [f for f in fnames_insitremote if f[-14:-6] in dates_dayconv]
fnames_nightconv = [f for f in fnames_insitremote if f[-14:-6] in dates_nightconv]
fnames_other = [f for f in fnames_insitremote 
                if f[-14:-6] not in dates_dayconv + dates_nightconv]
"""

def collect_p3data(datadir, fnames):
    p3all = None
    for f in fnames:
        print(f[-14:-3])
        """
        p3 = xr.load_dataset(f)
        if p3all is not None:
            p3all = xr.concat([p3all, p3], 'time', data_vars='all', 
                              coords='all', compat='identical', join='outer', 
                              combine_attrs='override')
        else:
            p3all = p3
        """
        p3 = xr.load_dataset(datadir+f)
            # Add potential temperature:
        p3['theta'] = ap.ptemp(p3['Ta'], p3['press'])
        p3_binned = p3.groupby(100*np.round(p3['alt']/100)).mean()
        if p3all is not None:
            p3all = xr.concat([p3all, p3_binned], 'cloud_module', data_vars='all', 
                              coords='all', compat='identical', join='outer', 
                              combine_attrs='override')
        else:
            p3all = p3_binned           

            
    return p3all
        

fnames_cldgroup1 = [f for f in fnames_cldinsituremote if int(f[-5:-3]) in ncld_g1]
data_g1 = collect_p3data(path_cldinsituremote_dir, fnames_cldgroup1)
fnames_cldgroup2 = [f for f in fnames_cldinsituremote if int(f[-5:-3]) in ncld_g2]
data_g2 = collect_p3data(path_cldinsituremote_dir, fnames_cldgroup2)
fnames_cldgroup3 = [f for f in fnames_cldinsituremote if int(f[-5:-3]) in ncld_g3]
data_g3 = collect_p3data(path_cldinsituremote_dir, fnames_cldgroup3)

#p3other = collect_p3data(fnames_other)
#p3dayconv = collect_p3data(fnames_dayconv)
#p3nightconv = collect_p3data(fnames_nightconv)



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
        

    # Plot:
fig = plt.figure(figsize=(6.5, 2.5))
w = 0.215
lmarge = 0.1
axset = [
    fig.add_axes([x, 0.2, w, 0.78]) 
    for x in [lmarge, lmarge + w + 0.01, lmarge + 2*w + 2*0.01, lmarge + 3*w + 3*0.01]
    ]


plot_prf_envelop(data_g1, axset, ['theta','mr','rh','dD'], 'grey')
plot_prf_envelop(data_g2, axset, ['theta','mr','rh','dD'], 'blue')
plot_prf_envelop(data_g3, axset, ['theta','mr','rh','dD'], 'red')
#plot_prf_envelop(p3other, axset, ['theta','mr','rh','dD'], 'blue')
#plot_prf_envelop(p3nightconv, axset, ['theta','mr','rh','dD'], 'grey')
#plot_prf_envelop(p3dayconv, axset, ['theta','mr','rh','dD'], 'red')


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


#fig.savefig(r'../figures/day-night_prfs.png')





















