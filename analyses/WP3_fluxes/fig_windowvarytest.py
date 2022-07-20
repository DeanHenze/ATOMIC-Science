# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 11:05:24 2022

@author: Dean

Test to see if window size changes flux computations.
"""

import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import profileplotter
import rangescaler
import prfrestructure
import thermo_deSzoeke
import oversampler



## flux file path and filenames:
path_fluxdir = "./fluxes_levlegs_windowvary/"
fnames_flux = os.listdir(path_fluxdir)


#ncld = '03'
ncld_list = [str(i).zfill(2) for i in range(1, 17)]
for ncld in ncld_list:
    fnames_cld = [f for f in fnames_flux if "_cld%s" % ncld in f] 
    fnames_cld.sort()

    results = []
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 8))
    for f in fnames_cld:
    
        window = f[f.index('_window') + 1: f.index('.nc')] 
        
        data = xr.load_dataset(path_fluxdir + f)
        results.append(data["w'w'_bar"].values)
        #ax1.scatter(data["w'w'_bar"], data["alt"], label=window)
        #ax2.scatter(data["q'w'_bar"], data["alt"], label=window)
    
    #plt.legend()
    
    results = np.array(results)
    results_normed  = (results - results[0, :])/results[0, :]
    
    x = np.arange(0, results.shape[0], 1)
    for j in range(results_normed.shape[1]):
        plt.plot(x, results_normed[:, j])
    #x = np.arange(1, results.shape[0], 1)
    #dt = 15
    #dydt = (results[1:, :] - results[:-1, :])/dt
    #for j in range(dydt.shape[1]):
    #    plt.plot(x*dt, dydt[:, j])











