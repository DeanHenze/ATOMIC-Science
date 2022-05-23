# -*- coding: utf-8 -*-
"""
Created on Thu May 19 10:29:51 2022

@author: Dean

Goal is to have some summary plots of vertical velocity PDFs for sub-cloud, 
cloud-base, in-cloud, and cloud top heights across all cloud modules.
"""



import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from scipy.stats import norm

# Local code
import thermo
import rangescaler



ncld = 4


# Level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]
fnames_levlegs_ncld = [f for f in fnames_levlegs if "_cld%i"%ncld in f]


# Key altitude info:
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
keyalts_ncld = keyalts_table.loc[keyalts_table['ncld']==ncld]

n_figrows = int(np.ceil(len(fnames_levlegs_ncld)/3))
fig, axes = plt.subplots(n_figrows, 3, figsize=(12,8))

for f, ax in zip(fnames_levlegs_ncld, axes.flatten()):
    
    data = xr.load_dataset(dir_5hzdata + f)
    data_df = data.to_dataframe() # Easier / faster to work with pandas.

    # Aircraft roll QC
    # Remove data where the roll was greater than 5 degrees:
    roll_crit = 5
    highroll = abs(data_df['roll']) > roll_crit
    data_dfqc = data_df.loc[~highroll]
    
    #x = data_dfqc["w'"]*data_dfqc["q'"]
    x = data_dfqc["w'"]
    
    # Level leg histogram:
    ax.hist(x, bins=25, density=True)
    
    # Gaussian for reference:
    def plot_gauss(x, ax):
        std = x.std()
        xmin = x.min()
        xmax = x.max()
        domain_bound = np.max([abs(xmin), xmax])
        xgauss = np.linspace(-domain_bound, domain_bound, 100) # symmetric domain
        xpdfgauss = norm.pdf(xgauss, loc=0, scale=std)
        ax.plot(xgauss, xpdfgauss)
        
    plot_gauss(x, ax)
        
    # 
    #ax.set_xlim(-2.5, 2.5)
    ax.text(
        0.9, 0.9, "z = %i m" % round(data_dfqc["alt"].mean()), 
        ha='right', va='top', transform=ax.transAxes
        )
    
    
    
    
    
    
    
    
    
    
    
    
    