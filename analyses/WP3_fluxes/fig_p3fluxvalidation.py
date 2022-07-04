# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:58:44 2022

@author: Dean
"""


import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


## I/O paths, filenames
##_____________________________________________________________________________
path_atrfluxdir = "../../data/ATR/"
fnames_atrflux = [f for f in os.listdir(path_atrfluxdir) if f.endswith(".nc")]
path_p3fluxdir = "./fluxes_levlegs/"
fnames_p3flux = [f for f in os.listdir(path_p3fluxdir) if f.endswith(".nc")]
##_____________________________________________________________________________
## I/O paths, filenames


## Collect ATR data into a pandas df
##_____________________________________________________________________________
alt_atr = np.array([])
wq_atr = np.array([])
wt_atr = np.array([])
ww_atr = np.array([])
uu_atr = np.array([])
vv_atr = np.array([])
for f in fnames_atrflux:
    data = xr.load_dataset(path_atrfluxdir + f)
    alt_atr = np.append(alt_atr, data['alt'].values)
    wq_atr = np.append(wq_atr, data['COVAR_WMR'].values)
    wt_atr = np.append(wt_atr, data['COVAR_WT'].values)
    ww_atr = np.append(ww_atr, data['VAR_W'].values)
    uu_atr = np.append(uu_atr, data['VAR_U'].values)
    vv_atr = np.append(vv_atr, data['VAR_V'].values)
dfkeys = ('alt', 'wq', 'wt', 'ww', 'vv', 'uu')    
atr_df = pd.DataFrame(
    dict(
        zip(
            dfkeys, 
            (alt_atr, wq_atr, wt_atr, ww_atr, vv_atr, uu_atr)
            )
        )
    )
atr_df['tke'] = 0.5*(atr_df['ww'] + atr_df['vv'] + atr_df['uu'])
##_____________________________________________________________________________
## Collect ATR data into a pandas df


## Collect P-3 data into a pandas df
##_____________________________________________________________________________
alt_p3 = np.array([])
wq_p3 = np.array([])
wt_p3 = np.array([])
ww_p3 = np.array([])
uu_p3 = np.array([])
vv_p3 = np.array([])
for f in fnames_p3flux:
    data = xr.load_dataset(path_p3fluxdir + f)
    alt_p3 = np.append(alt_p3, data['alt'].values)
    wq_p3 = np.append(wq_p3, data["q'w'_bar"].values)
    wt_p3 = np.append(wt_p3, data["T'w'_bar"].values)
    ww_p3 = np.append(ww_p3, data["w'w'_bar"].values)
    uu_p3 = np.append(uu_p3, data["u'u'_bar"].values)
    vv_p3 = np.append(vv_p3, data["v'v'_bar"].values)
p3_df = pd.DataFrame(
    dict(
        zip(
            dfkeys,
            (alt_p3, wq_p3, wt_p3, ww_p3, vv_p3, uu_p3)
            )
        )
    )
p3_df['tke'] = 0.5*(p3_df['ww'] + p3_df['vv'] + p3_df['uu'])
##_____________________________________________________________________________
## Collect P-3 data into a pandas df
    

## Scatter plots
##_____________________________________________________________________________
pltvarkeys = ('ww', 'tke', 'wq', 'wt')    
fig = plt.figure(figsize=(6.5, 3))
ax1 = fig.add_axes([0.1, 0.175, 0.21, 0.75])
ax2 = fig.add_axes([0.325, 0.175, 0.21, 0.75])
ax3 = fig.add_axes([0.55, 0.175, 0.21, 0.75])
ax4 = fig.add_axes([0.775, 0.175, 0.21, 0.75])
axset = [ax1, ax2, ax3, ax4]
for k, ax in zip(pltvarkeys, axset):
    ax.scatter(
        atr_df[k], atr_df['alt'], 
        color='lightgrey', s=2, label='ATR'
        )
    ax.scatter(
        p3_df[k], p3_df['alt'], 
        color='black', s=2, label='P-3'
        )
##_____________________________________________________________________________
## Scatter plots


## Axes limits, labels, legend, ...
##_____________________________________________________________________________
axset[0].set_xlim(-0.025, 0.84)
axset[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
axset[0].set_xticklabels(axset[0].get_xticks().astype(str), fontsize=9)
axset[0].set_xlabel(r"$\bar{w'w'}$ (m$^{2}$ s$^{-2}$)", fontsize=12)

axset[1].set_xlim(-0.025, 1.3)
axset[1].set_xticks([0, 0.25, 0.5, 0.75, 1, 1.25])
axset[1].set_xticklabels(['0', '', '0.5', '', '1', ''], fontsize=9)
axset[1].set_xlabel(r"TKE (m$^{2}$ s$^{-2}$)", fontsize=12)

axset[2].set_xlim(-0.05, 0.255)
axset[2].set_xticks([0, 0.05, 0.1, 0.15, 0.2])
axset[2].set_xticklabels(['0', '', '0.1', '', '0.2'], fontsize=9)
axset[2].set_xlabel(r"$\bar{w'q'}$ (g kg$^{-1}$ m s$^{-1}$)", fontsize=12)

axset[3].set_xlim(-0.03, 0.03)
axset[3].set_xticks([-0.025, -0.0125, 0, 0.0125, 0.025])
axset[3].set_xticklabels(['-0.025', '', '0', '', '0.025'], fontsize=9)
axset[3].set_xlabel(r"$\bar{w'\theta'}$ (K m s$^{-1}$)", fontsize=12)

for ax in axset: 
    ax.set_yticks(np.arange(0, 3001, 500))
    ax.set_ylim(-100, 3200)
axset[0].set_yticklabels(axset[0].get_yticks(), fontsize=9)
for ax in axset[1:]: 
    ax.set_yticklabels(['' for t in ax.get_yticks()])
axset[0].set_ylabel('altitude (m)', fontsize=12)

axset[2].legend(loc='upper right', fontsize=9, markerscale=2)
##_____________________________________________________________________________
## Axes limits, labels, legend, ...


fig.savefig("./fig_p3fluxvalidation.png")


"""
altbin = 200
atr_binavg = atr_df.groupby(
    np.round(atr_df['alt']/altbin)*altbin, axis=0, as_index=True)
atr_binavg = atr_binavg.mean()
p3_binavg = p3_df.groupby(
    np.round(p3_df['alt']/altbin)*altbin, axis=0, as_index=True)
p3_binavg = p3_binavg.mean()

atr_binavg = atr_binavg.rolling(
        window=4, min_periods=1, win_type='hamming', 
        center=True, closed='neither'
        )
atr_binavg = atr_binavg.mean()
p3_binavg = p3_binavg.rolling(
        window=4, min_periods=1, win_type='hamming', 
        center=True, closed='neither'
        )
p3_binavg = p3_binavg.mean()


plt.plot(atr_binavg['wq'], atr_binavg['alt'], color='black', linewidth=4)
plt.plot(p3_binavg['wq'], p3_binavg['alt'], color='lightgrey', linewidth=4)
"""
