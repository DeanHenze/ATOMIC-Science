# -*- coding: utf-8 -*-
"""
Created on Sun May 29 12:41:53 2022

@author: Dean

Figure of mososcale and synoptic characteristics for each P-3 cloud module.
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
import thermo



path_era5dir = "./cldmod_datafiles/"
path_p3insitu = "./cldmod_datafiles/"
fnames_p3insitu = [f for f in os.listdir(path_p3insitu) 
                   if '_insitu+remote_' in f]


varkeys = ['LTS', 'wspd_sfc', 'wdir_sfc', 
           'wdiv_950hPa', 'wdiv_850hPa', 'wdiv_750hPa'
           ]
mesomean = dict(zip(varkeys, [[] for e in varkeys]))
mesostd = dict(zip(varkeys, [[] for e in varkeys]))
nearp3 = dict(zip(varkeys, [[] for e in varkeys]))


windspeed = lambda u, v: (u**2 + v**2)**0.5
winddirection = lambda u, v: np.arctan(u/v)*(180/np.pi)
winddiff = lambda u1, u2, v1, v2: ((u2-u1)**2 + (v2-v1)**2)**0.5


ncld_g1 = [7, 9, 11, 10, 12, 6]
ncld_g2 = [8, 15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4]
#ncld_list = np.arange(1, 17, 1)
ncld_list = list(np.arange(1, 17, 1))

#ncld_list = [str(i).zfill(2) for i in np.arange(1, 17, 1)]
for ncld in ncld_list:
    
    # Load data:
    ncld_str = str(ncld).zfill(2)
    fname_era5 = "p3cld_ECMWF_ERA5_plevs_hourly_ncld%s.nc" % ncld_str
    era5 = xr.load_dataset(path_era5dir + fname_era5)
    fname_p3 = [f for f in fnames_p3insitu if "ncld%s" % ncld_str in f][0]
    p3 = xr.load_dataset(path_p3insitu + fname_p3)
    
    # Compute additional vars for ERA5:
    era5['theta'] = thermo.theta(era5['t'], era5['level']*100)
    era5_ns = era5.sel(level=1000) # near surface
    era5['wdiv_950hPa'] = era5['d'].sel(level=950)
    era5['wdiv_850hPa'] = era5['d'].sel(level=850)
    era5['wdiv_750hPa'] = era5['d'].sel(level=750)
    era5['LTS'] = era5['theta'].sel(level=700) - era5_ns['theta']
    era5['wspd_sfc'] = windspeed(era5_ns['u'], era5_ns['v'])
    era5['wdir_sfc'] = winddirection(era5_ns['u'], era5_ns['v'])
    
    # ERA5 grid point nearest the mean lat/lon of P-3 cloud module:
    p3mean = p3.mean(dim='time')
    era5_nearp3 = era5.sel(
        longitude=p3mean['lon'], latitude=p3mean['lat'],
        method='nearest'
        )
    era5_nearp3 = era5_nearp3.mean(dim='time')
    
    era5mean = era5.mean(dim=['longitude','latitude','time'])
    era5std = era5.std(dim=['longitude','latitude','time'])
    
    for k in varkeys:
        mesomean[k].append(era5mean[k].item())
        mesostd[k].append(era5std[k].item())
        nearp3[k].append(era5_nearp3[k].item())
    
    
for k in nearp3.keys():
    nearp3[k] = np.array(nearp3[k])

#for dc in mesoscalemean

def mesoscale_shading(x, mesomean_dict, mesostd_dict, varkey, ax):
    """
    """
    ymean = np.array(mesomean_dict[varkey])
    ystd = np.array(mesostd_dict[varkey])
    y1 = ymean-ystd
    y2 = ymean+ystd
    ax.fill_between(x, y1, y2, color='grey', alpha=0.25)


fig_stab = plt.figure(figsize=(6.5, 4))
ax_lts = fig_stab.add_axes([0.15, 0.15, 0.7, 0.8])
mesoscale_shading(ncld_list, mesomean, mesostd, 'LTS', ax_lts)
colors=['grey', 'blue', 'red']
for nlist, c in zip([ncld_g1, ncld_g2, ncld_g3], colors):
    i_group = [ncld_list.index(i) for i in nlist]
    ax_lts.scatter(
        nlist, nearp3['LTS'][i_group], marker='o', 
        color=c, label='near p3'
        )
    
ax_sfcwspd = ax_lts.twinx()
mesoscale_shading(ncld_list, mesomean, mesostd, 'wspd_sfc', ax_sfcwspd)
for nlist, c in zip([ncld_g1, ncld_g2, ncld_g3], colors):
    i_group = [ncld_list.index(i) for i in nlist]
    ax_sfcwspd.scatter(
        nlist, nearp3['wspd_sfc'][i_group], marker='^', 
        color=c, label='near p3'
        )
ax_lts.set_xticks(ncld_list)
ax_lts.set_xlabel("cloud module", fontsize=12)
ax_lts.set_ylabel("LTS (K)", fontsize=12)
ax_sfcwspd.set_ylabel(r"$U_{sfc}$ (m/s)", fontsize=12)



fig_div = plt.figure(figsize=(6.5, 4))
ax_div = fig_div.add_axes([0.15, 0.15, 0.7, 0.8])
mesoscale_shading(ncld_list, mesomean, mesostd, 'wdiv_950hPa', ax_div)

for nlist, c in zip([ncld_g1, ncld_g2, ncld_g3], colors):
    i_group = [ncld_list.index(i) for i in nlist]
    ax_div.scatter(
        nlist, nearp3['wdiv_950hPa'][i_group], marker='o', 
        color=c, label='near p3'
        )
ax_div.set_xticks(ncld_list)
ax_div.hlines(0, 1, 16, colors='black', linestyles='--')
ax_div.set_xlabel("cloud module", fontsize=12)
ax_div.set_ylabel(r"$D_{950}$ ($s^{-1}$)", fontsize=12)


fig_wdir = plt.figure(figsize=(6.5, 4))
ax_wdir = fig_wdir.add_axes([0.15, 0.15, 0.7, 0.8])
mesoscale_shading(ncld_list, mesomean, mesostd, 'wdir_sfc', ax_wdir)
ax_wdir.scatter(
    ncld_list, nearp3['wdir_sfc'], marker='o', 
    color='black', label='near p3'
    )
for nlist, c in zip([ncld_g1, ncld_g2, ncld_g3], colors):
    i_group = [ncld_list.index(i) for i in nlist]
    ax_wdir.scatter(
        nlist, nearp3['wdir_sfc'][i_group], marker='o', 
        color=c, label='near p3'
        )
ax_wdir.set_xticks(ncld_list)
ax_wdir.hlines(0, 1, 16, colors='black', linestyles='--')
ax_wdir.set_xlabel("cloud module", fontsize=12)
ax_wdir.set_ylabel(r"$wdir_{sfc}$ (deg)", fontsize=12)



"""
ax_lts.errorbar(
    ncld_list, mesoscalemean['LTS'], yerr=mesoscalestd['LTS'], 
    marker='o', color='black', label='LTS meso'
    )
ax_lts.scatter(
    ncld_list, nearp3['LTS'], marker='*', 
    color='black', label='LTS near p3'
    )


fig_wind = plt.figure(figsize=(6.5, 4))
ax_sfcwspd = fig_wind.add_axes([0.15, 0.15, 0.7, 0.8])
ax_sfcwdir = ax_sfcwspd.twinx()
ax_sfcwspd.errorbar(
    ncld_list, mesoscalemean['wspd_sfc'], yerr=mesoscalestd['wspd_sfc'], 
    marker='o', color='black'
    )
ax_sfcwspd.scatter(
    ncld_list, nearp3['wspd_sfc'], marker='*', 
    color='black'
    )
"""

"""
fig_lsforce = plt.figure(figsize=(6.5, 4))
ax_lts = fig_lsforce.add_axes([0.15, 0.15, 0.7, 0.8])
ax_sfcwspd = ax_lts.twinx()

ax_lts.scatter(
    ncld_list, lts, marker='o', 
    color='black', label='LTS'
    )
ax_lts.errorbar(
    ncld_list, lts, yerr=lts_std, 
    marker='o', color='black', label='LTS'
    )
#ax_sfcwspd.scatter(
#    ncld_list, sfcwspd, 
#    marker='^', color='grey', label=r'$U_{sfc}$'
#    )

#ax_lts.set_ylim(8.5, 13)
ax_sfcwspd.set_ylim(5, 13)
ax_lts.set_xlabel("cloud module #", fontsize=12)
ax_lts.set_ylabel("LTS (K)", fontsize=12)
ax_sfcwspd.set_ylabel(r"$U_{sfc}$ (m/s)", fontsize=12)
ax_lts.legend(loc='upper left')
ax_sfcwspd.legend(loc='upper right')

fig_lsforce.savefig("./fig_lsforce_cldmods.png")
"""

#fig = plt.figure(figsize=(6.5, 4))
#ax_qlapse = fig.add_axes([0.15, 0.15, 0.8, 0.8])
#ax_qlapse.scatter(ncld_list, qlapse, marker='o', color='black')


#fig = plt.figure(figsize=(6.5, 4))
#ax_tsfc = fig.add_axes([0.15, 0.15, 0.8, 0.8])
#ax_tstd = ax_tsfc.twinx()
#ax_tsfc.scatter(ncld_list, sfctemp_mean, marker='o', color='black')
#ax_tstd.scatter(ncld_list, sfctemp_std, marker='^', color='grey')


#fig = plt.figure(figsize=(6.5, 4))
#ax_wnddif = fig.add_axes([0.15, 0.15, 0.8, 0.8])
#ax_wnddif.scatter(ncld_list, wnddif, marker='o', color='black')

















