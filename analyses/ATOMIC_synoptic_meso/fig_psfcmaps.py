# -*- coding: utf-8 -*-
"""
Created on Mon May 30 07:46:23 2022

@author: Dean

Surface pressure anomaly maps for Jan 22-24.
P-3 cloud sampling locations overlain.
"""



# Built in
import os

# Third party
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs



## Set basemap
##_____________________________________________________________________________
def set_basemapaxes(fig, axpos):
    ax = fig.add_axes(axpos, projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_xlim()
    ax.set_extent([-60, -48, 9, 24.5], crs=None)
    #gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    #gl.xlines=False
    #gl.ylines=False
    return ax
##_____________________________________________________________________________
## Set basemap



## Surface pressure anomaly contour maps
##_____________________________________________________________________________
def dpsfc_contour(ax, fname_dpsfc):
    path_era5psfc = "./era5_Psfc/"
    dpsfc = xr.load_dataset(path_era5psfc + fname_dpsfc)
    cntr_dp = ax.contour(
        dpsfc['longitude'], dpsfc['latitude'],
        dpsfc['sp'], 
        colors='black', vmin=-8, vmax=0
        )
    ax.clabel(cntr_dp, fmt='%0.1f', fontsize=9)
##_____________________________________________________________________________



## Plot P-3 cloud module locations
##_____________________________________________________________________________
def p3_cldlocs(ax):
    path_p3insitu = "../WP3_cloudmodule_char/cldmod_datafiles/"
    fnames_p3insitu = [f for f in os.listdir(path_p3insitu) 
                       if '_insitu+remote_' in f]

    ncld_list = [str(i).zfill(2) for i in np.arange(1, 17, 1)]
    for ncld in ncld_list:
        
        # Load p3 cloud module data:
        fname_p3 = [f for f in fnames_p3insitu if "ncld%s" % ncld in f][0]
        p3 = xr.load_dataset(path_p3insitu + fname_p3)
    
        # Plot mean location:
        p3mean = p3.mean(dim='time')
        ax.scatter(p3mean['lon'], p3mean['lat'], marker='^', color='black', s=15)
##_____________________________________________________________________________    
## Plot P-3 cloud module locations
    


## Axes ticks, labels
##_____________________________________________________________________________
def ticks_labels(ax):
    xticks = [-60 , -57, -54 , -51, -48]
    ax.set_xticks(xticks)
    ax.set_xticklabels([t + u'\N{DEGREE SIGN}' for t in ax.get_xticks().astype(str)])
    
    yticks = [10, 12, 14, 16, 18, 20, 22, 24]
    ax.set_yticks(yticks)
    ax.set_yticklabels([t + u'\N{DEGREE SIGN}' for t in ax.get_yticks().astype(str)])
##_____________________________________________________________________________
## Axes ticks, labels



## Make and save figures
##_____________________________________________________________________________
fig1 = plt.figure(figsize=(2, 2))
ax1pos = [0.15, 0.15, 0.8, 0.8]
ax1 = set_basemapaxes(fig1, ax1pos)
fname_dpsfc_jan22 = "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200122.nc"
dpsfc_contour(ax1, fname_dpsfc_jan22)
p3_cldlocs(ax1)
ax1.text(0.05, 0.05, 'Jan. 22', fontsize=11, transform=ax1.transAxes, color='b')
ticks_labels(ax1)
fig1.savefig("./fig_dpsfc_jan22.png")


fig2 = plt.figure(figsize=(2, 2))
ax2pos = [0.15, 0.15, 0.8, 0.8]
ax2 = set_basemapaxes(fig2, ax2pos)
fname_dpsfc_jan23 = "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200123.nc"
dpsfc_contour(ax2, fname_dpsfc_jan23)
p3_cldlocs(ax2)
ax2.text(0.05, 0.05, 'Jan. 23', fontsize=11, transform=ax2.transAxes, color='b')
ticks_labels(ax2)
fig2.savefig("./fig_dpsfc_jan23.png")


fig3 = plt.figure(figsize=(2, 2))
ax3pos = [0.15, 0.15, 0.8, 0.8]
ax3 = set_basemapaxes(fig3, ax3pos)
fname_dpsfc_jan24= "ECMWF_ERA5_psfcmap_ATOMIC-region_mean_20200124.nc"
dpsfc_contour(ax3, fname_dpsfc_jan24)
p3_cldlocs(ax3)
ax3.text(0.05, 0.05, 'Jan. 24', fontsize=11, transform=ax3.transAxes, color='b')
ticks_labels(ax3)
fig3.savefig("./fig_dpsfc_jan24.png")
##_____________________________________________________________________________











