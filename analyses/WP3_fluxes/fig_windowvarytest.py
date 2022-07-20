# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:05:24 2022

@author: Dean

Test to see if window size changes flux computations.
"""



import os

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt



## flux file path and filenames:
path_fluxdir = "./fluxes_levlegs_windowvary/"
fnames_flux = os.listdir(path_fluxdir)


results_normed_arr = np.array([])
results_normed_list1 = []
results_normed_list2 = []
ncld_list = [str(i).zfill(2) for i in range(1, 17)]


#fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 3))
for ncld in ncld_list:

    fnames_cld = [f for f in fnames_flux if "_cld%s" % ncld in f] 
    fnames_cld.sort()


    results = []
    for f in fnames_cld:
        window = f[f.index('_window') + 1: f.index('.nc')] 
        data = xr.load_dataset(path_fluxdir + f)
        if data["q'w'_bar"].mean() < 4e-06: continue
        results.append(data["q'w'_bar"].values)
    results = np.array(results)
    
    
    results_normed1  = (results[:6, :] - results[0, :])/results[0, :]
    x = np.arange(0, results_normed1.shape[0], 1)
    for j in range(results_normed1.shape[1]):
        #ax1.scatter(x, results_normed1[:, j], color='darkgrey')    
        results_normed_list1.append(results_normed1[:, j])


    results_normed2  = (results[5:, :] - results[5, :])/results[5, :]
    x = np.arange(0, results_normed2.shape[0], 1)
    for j in range(results_normed2.shape[1]):
        #ax2.scatter(x, results_normed2[:, j], color='darkgrey')
        results_normed_list2.append(results_normed2[:, j])

    
    results_normed_list = np.append(
        results_normed_arr, 
        results_normed2.flatten()
        )


windows = [int(f[i1:i2]) for f, i1, i2 
           in zip(
               fnames_flux, 
               [f.index('_window') + 7 for f in fnames_flux], 
               [f.index('s.nc') for f in fnames_flux]
               )
          ]
windows_unique = np.unique(windows)
windows_unique.sort()

"""
# PDF of deviations with fit:
from scipy.stats import norm
plt.figure()
plt.hist(results_normed_arr, bins=50, density=True, color='darkgrey')

    # Fit a normal distribution and plot:
mu, sig = norm.fit(results_normed_list)
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, sig)
plt.plot(x, p, 'k', linewidth=2)

ax1.set_xticks(range(0, 6))
ax1.set_xticklabels([str(i) for i in windows_unique[0:6]], fontsize=10)
"""

fig = plt.figure(figsize=(6.5, 2.5))
ax1 = fig.add_axes([0.1, 0.2, 0.375, 0.75])
ax2 = fig.add_axes([0.6, 0.2, 0.375, 0.75])

ax1.boxplot(np.array(results_normed_list1), whis=(5, 95), showfliers=False)
ax1.hlines(0.3, 0.5, 6.5, colors='black', linestyles='dashed')

ax2.boxplot(np.array(results_normed_list2), whis=(5, 95), showfliers=False)

ax1.set_xticks(range(1, 7))
ax1.set_xticklabels([str(i) for i in windows_unique[0:6]], fontsize=10)
ax1.set_xlabel("W (s)", fontsize=11)
ax1.set_ylim(-1.5, 2.1)
ax1.set_ylabel(r"$f_{30s}(W, \bar{w'q'})$", fontsize=11)

ax2.set_xticks(range(1, 14))
ax2.set_xticklabels(
    [str(i) if i%2==0 else '' for i in windows_unique[6:]], 
    fontsize=10
    )
ax2.set_xlabel("W (s)", fontsize=11)
ax2.set_ylim(-1.5, 2.1)
ax2.set_ylabel(r"$f_{105s}(W, \bar{w'q'})$", fontsize=11)


fig.savefig("./fig_spectrawindowvary.png")





