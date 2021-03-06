# -*- coding: utf-8 -*-
"""
Created on Sun May 29 12:41:53 2022

@author: Dean

Figure of mososcale characteristics for each P-3 cloud module.
Expressed as pair plots of mesoscale circulation variables.

Current status:
--------------
- Pending final code cleanup (including removing code for time series plots, 
  which likely won't end up in the manuscript.)
- question of whether or not to add pair plot of 600 hPa vertical velocity 
  vs. something.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# Local code
import thermo



path_era5dir = "./era5_cldmods/"
path_p3insitu = "../WP3_cloudmodule_char/cldmod_datafiles/"
fnames_p3insitu = [f for f in os.listdir(path_p3insitu) 
                   if '_insitu+remote_' in f]


varkeys = ['LTS', 'wspd_sfc', 'wdir_sfc', 
           'wdiv_950hPa', 'wdiv_850hPa', 'wdiv_750hPa'
           ]


# Wind related functions:
windspeed = lambda u, v: (u**2 + v**2)**0.5
winddirection = lambda u, v: np.arctan(u/v)*(180/np.pi)
windshear = lambda u1, u2, v1, v2, dz: (1/dz)*((u2-u1)**2 + (v2-v1)**2)**0.5


#ncld_g1 = [1, 5, 4]
#ncld_g2 = [8, 7, 9, 11, 10, 6]
#ncld_g3= [12, 15, 3, 2, 13, 16, 14]
ncld_g1 = [1, 5, 4, 8]
ncld_g2 = [7, 9, 11, 10, 6]
ncld_g3= [12, 15, 13, 16, 14]
ncld_g3_lowp = [3, 2]


def lsforce_tseries(ncld, pltcolor, axset):
    """
    """
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
    era5['LTS'] = era5['theta'].sel(level=650) - era5_ns['theta']
    era5['qlapse'] = era5['q'].sel(level=700) - era5_ns['q']
    era5['wspd_sfc'] = windspeed(era5_ns['u'], era5_ns['v'])
    era5['wdir_sfc'] = winddirection(era5_ns['u'], era5_ns['v'])
        # Windshear between surface and ~3km:
    era5['windshear_sfc700hpa'] = windshear(
        era5['u'].sel(level=1000), era5['u'].sel(level=700), 
        era5['v'].sel(level=1000), era5['v'].sel(level=700), 
        (1000-700)*10
        )
    
    
    # ERA5 grid point nearest the mean lat/lon of P-3 cloud module:
    #p3mean = p3.mean(dim='time')
    #era5_nearp3 = era5.sel(
    #    longitude=p3mean['lon'], latitude=p3mean['lat'],
    #    method='nearest'
    #    )
    #era5_nearp3 = era5_nearp3.mean(dim='time')
    

    era5mean_ts = era5.mean(dim=['longitude','latitude'])
    dt = np.flip(era5mean_ts['time'][-1] - era5mean_ts['time'])/(60*10**9)
    dt.values = dt.values.astype(float)
    era5mean_ts = era5mean_ts.where(dt<=120, drop=True)
    dt = dt[dt<=120]
    
        
    axset[0].plot(dt, era5mean_ts['wdiv_950hPa'], color=pltcolor)
    axset[1].plot(dt, era5mean_ts['wspd_sfc'], color=pltcolor)
    axset[2].plot(dt, era5mean_ts['windshear_sfc700hpa'], color=pltcolor)
    axset[3].plot(dt, era5mean_ts['LTS'], color=pltcolor)
    
    
    
def lsforce_pairplot(ncld, pltcolor, axset, pltmarker='o'):
    """
    """
    # Load data:
    ncld_str = str(ncld).zfill(2)
    fname_era5 = "p3cld_ECMWF_ERA5_plevs_hourly_ncld%s.nc" % ncld_str
    era5 = xr.load_dataset(path_era5dir + fname_era5)
    
    # Compute additional vars for ERA5:
    era5['theta'] = thermo.theta(era5['t'], era5['level']*100)
    era5_ns = era5.sel(level=1000) # near surface
    era5['wdiv_950hPa'] = era5['d'].sel(level=950)
    era5['w_600hPa'] = era5['w'].sel(level=600)
    era5['LTS'] = era5['theta'].sel(level=650) - era5_ns['theta']
    era5['qlapse'] = era5['q'].sel(level=700) - era5_ns['q']
    era5['wspd_sfc'] = windspeed(era5_ns['u'], era5_ns['v'])
    era5['wdir_sfc'] = winddirection(era5_ns['u'], era5_ns['v'])
        # Windshear between surface and ~3km:
    era5['windshear_sfc700hpa'] = windshear(
        era5['u'].sel(level=1000), era5['u'].sel(level=700), 
        era5['v'].sel(level=1000), era5['v'].sel(level=700), 
        (1000-700)*10
        )
    
    
    era5mean = era5.mean(dim=['longitude','latitude','time'])
    #era5std = era5.std(dim=['longitude','latitude','time'])  
    axset[0].plot(
        era5mean['wspd_sfc'], era5mean['wdiv_950hPa'], 
        marker='o', markersize=4, color=pltcolor
        )
    axset[1].plot(
        era5mean['windshear_sfc700hpa'], era5mean['LTS'], 
        marker=pltmarker, markersize=4, color=pltcolor
        )
    axset[2].plot(
        era5mean['LTS'], era5mean['wspd_sfc'], 
        marker=pltmarker, markersize=4, color=pltcolor
        )  
    axset[3].plot(
        era5mean['wspd_sfc'], era5mean['w_600hPa'], 
        marker=pltmarker, markersize=4, color=pltcolor
        )
    
    
    return (
        [era5mean['wspd_sfc'].item(), era5mean['wdiv_950hPa'].item(), 
        era5mean['LTS'].item(), era5mean['w_600hPa'].item()],
        
        (r'$|U|_{sfc}$ (m/s)', r'$D_{950}$ (s$^{-1}$ * 10$^{-5}$)', 
         'LTS (K)', r'$\omega_{600}$ (Pa s$^{-1}$)')
        )

    
    
    
"""
## Time series plot
##_____________________________________________________________________________
# Plot:
fig, axset = plt.subplots(4, 1, figsize=(6.5, 7))
ncldgroups = [ncld_g1, ncld_g2, ncld_g3]
pltcolors = ['red', 'grey', 'blue']
for ngroup, c in zip(ncldgroups, pltcolors):
    for n in ngroup:
        lsforce_tseries(n, c, axset)


# Axes labels, limits, legend:
for ax in axset: 
    ax.invert_xaxis()
    ax.set_xticklabels((ax.get_xticks()*-1).astype(int).astype(str))

axset[-1].set_xlabel('hours before P-3 sampling', fontsize=12)

axset[0].set_ylabel(r'$div_{950}$ (s$^{-1}$ * 10$^{-5}$)', fontsize=12)
axset[0].set_yticks(10**(-5)*np.array([-1, 0, 1]))
axset[0].set_yticklabels((axset[0].get_yticks()*10**5).astype(int).astype(str))
axset[1].set_ylabel(r'$|U|_{sfc}$ (m/s)', fontsize=12)
axset[2].set_ylabel(r'$(dU/dz)_{h}$ (s$^{-1}$)', fontsize=12)
axset[3].set_ylabel(r'LTS (K)', fontsize=12)

legend_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='blue', lw=4)]
axset[0].legend(
    legend_lines, ['cg1', 'cg2', 'cg3'], 
    loc='lower left', bbox_to_anchor= (0.5, 1.01), ncol=3,
    borderaxespad=0, frameon=False
    )


fig.savefig("./fig_mesoscale_timeseries.png")
##_____________________________________________________________________________
## Timeseries plot
"""


## Pair plot
##_____________________________________________________________________________
# Plot:
fig, axset = plt.subplots(1, 4, figsize=(6.5, 2))
#fig = plt.figure(figsize=(6.5, 2.5))
#ax1 = fig.add_axes([0.1, 0.2, 0.225, 0.7])
#ax2 = fig.add_axes([0.425, 0.2, 0.225, 0.7])
#ax3 = fig.add_axes([0.75, 0.2, 0.225, 0.7])
#axset = (ax1, ax2, ax3)

results = []
n_cldgrouplabels = ['1', '1', '1', '1', '2', '2', '2', '2', '2', 
              '3', '3', '3', '3', '3', '3*', '3*']
ncldgroups = [ncld_g1, ncld_g2, ncld_g3, ncld_g3_lowp]
vkeys = None
pltcolors = ['red', 'grey', 'blue', 'blue']
pltmarkers = ['o', 'o', 'o', 'x']
for ngroup, c, m in zip(ncldgroups, pltcolors, pltmarkers):
    for n in ngroup:
        results_n, vkeys = lsforce_pairplot(n, c, axset, pltmarker=m)
        results.append(results_n)
        

# Axes labels, limits, legend:
axset[0].set_xlabel(r'$|U|_{sfc}$ (m/s)', fontsize=12)
axset[0].set_xticks([6, 8, 10, 12])
axset[0].set_xticklabels(['6', '8', '10', '12'], fontsize=9)
axset[0].set_ylabel(r'$div_{950}$ (s$^{-1}$ * 10$^{-5}$)', fontsize=12)
axset[0].set_yticks(10**(-5)*np.array([-1, 0, 1]))
axset[0].set_yticklabels((axset[0].get_yticks()*10**5).astype(int).astype(str))

axset[1].set_xlabel(r'$(dU/dz)_{h}$ (s$^{-1}$)', fontsize=12)
axset[1].set_ylabel(r'LTS (K)', fontsize=12)
axset[1].set_xlim(0.0005, 0.0032)
axset[1].set_xticks([0.001, 0.002, 0.003])
axset[1].set_xticklabels(axset[1].get_xticks().astype(str), fontsize=9)
axset[1].set_yticks(np.arange(15, 21, 1))
axset[1].set_yticklabels(axset[1].get_yticks().astype(str), fontsize=9)
"""
axset[2].set_xlabel(r'LTS (K)', fontsize=12)
axset[2].set_ylabel(r'$\omega_{600}$ (Pa * s$^{-1}$)', fontsize=12)
axset[2].set_xticks(np.arange(15, 20, 2))
axset[2].set_xticklabels(axset[2].get_xticks().astype(str), fontsize=9)
axset[2].set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
axset[2].set_yticklabels(['-0.1','','0','','0.1'], fontsize=9)
"""
axset[2].set_xlabel(r'LTS (K)', fontsize=12)
axset[2].set_ylabel(r'$|U|_{sfc}$ (m/s)', fontsize=12)
axset[2].set_xticks(np.arange(15, 20, 2))
axset[2].set_xticklabels(axset[2].get_xticks().astype(str), fontsize=9)
axset[2].set_yticks([6, 8, 10, 12])
axset[2].set_yticklabels(['6', '8', '10', '12'], fontsize=9)

legend_lines = [Line2D([0], [0], color='red', lw=4),
                Line2D([0], [0], color='grey', lw=4),
                Line2D([0], [0], color='blue', lw=4)]
axset[0].legend(
    legend_lines, ['cg1', 'cg2', 'cg3'], 
    loc='lower left', bbox_to_anchor= (1.0, 1.01), ncol=3,
    borderaxespad=0, frameon=False
    )


fig.savefig("./fig_mesoscale_pairplot.png")
##_____________________________________________________________________________
## Pair plot



## Try out a seaborn pairplot:
results_df = pd.DataFrame(results, columns=vkeys)  
results_df['CG#'] = n_cldgrouplabels

palette = ['tab:red', 'tab:grey', 'tab:blue', 'tab:blue']
#sns.set(rc={"figure.figsize":(2, 2)}) #width=8, height=4
#fig = plt.figure(figsize=(4, 4))
g = sns.pairplot(
    data=results_df, hue='CG#', 
    markers=['o', 'o', 'o', 'D'], palette=palette, 
    diag_kind='None', corner=True, 
    #height=6, aspect=1.5, #figsize=(4, 4)
    )
def hide_current_axis(*args, **kwds):
    ax = plt.gca()
    ax.set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
g.map_diag(hide_current_axis)

# Divergence ticks, labels:
ticks_D = np.array([-1.5e-05, -1.0e-05, -5.0e-06,  
                    0.0e+00,  5.0e-06, 1.0e-05,  1.5e-05])
g.axes[1,0].set_yticks(ticks_D)
g.axes[1,0].set_yticklabels((ticks_D/10**-5).astype(str), fontsize=10)
g.axes[3,1].set_xticks(ticks_D)
g.axes[3,1].set_xticklabels((ticks_D/10**-5).astype(str), fontsize=10)

# LTS ticks, labels:
ticks_lts = np.array([15, 17, 19])
g.axes[2,0].set_yticks(ticks_lts)
g.axes[2,0].set_yticklabels((ticks_lts).astype(str), fontsize=10)
g.axes[3,2].set_xticks(ticks_lts)
g.axes[3,2].set_xticklabels((ticks_lts).astype(str), fontsize=10)

# Subsidence rate ticks, labels:
ticks_omega = np.array([-0.1, -0.05, 0, 0.05, 0.1])
g.axes[3,0].set_yticks(ticks_omega)
g.axes[3,0].set_yticklabels((ticks_omega).astype(str), fontsize=10)

# Surface windspeed ticks, labels:
ticks_lts = np.array([6, 8, 10, 12])
g.axes[3,0].set_xticks(ticks_lts)
g.axes[3,0].set_xticklabels((ticks_lts).astype(str), fontsize=10)




## Correlation heatmap:
fig_corr = plt.figure(figsize=(3, 3))
ax = fig_corr.add_axes([0.2, 0.2, 0.75, 0.75])
sns.heatmap(
    results_df.corr(), annot=True, cbar=False, ax=ax, 
    xticklabels=(r'$|U|_{sfc}$', r'$D_{950}$', 'LTS', r'$\omega_{600}$'),
    yticklabels=(r'$|U|_{sfc}$', r'$D_{950}$', 'LTS', r'$\omega_{600}$')    
    )
fig_corr.savefig("./fig_lsforce_correlations.png")

        

