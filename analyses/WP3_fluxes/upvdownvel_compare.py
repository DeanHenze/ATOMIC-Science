# -*- coding: utf-8 -*-
"""
Created on Sat May 14 16:49:36 2022

@author: Dean

Compare updraft and downdraft q, dD distributions.

Vertical wind PDF is not symmetric, the upward tail is larger than downward 
tail. Maybe this is capturing updrafts. 


Next steps: 
    - try taking the top and bottom w' quartiles for comparison.
    - Could also do a bandpass filter for the frequencies with the most 
        cospectral power.
    - Could try seeing if the larger tails for upward wind are only for the 
    day time clouds. Maybe night time is reversed.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# Local code
import iso



# All level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]


def plots(ncld):
    
    fnames_levlegs_cld = [f for f in fnames_levlegs if "_cld%i"%ncld in f]
   
    for f in fnames_levlegs_cld:
        
        data = xr.load_dataset(dir_5hzdata + f)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        
        # Add dD column:
        data_df["dD"] = iso.qD_dD_convert('qD2dD', data_df["q"], data_df["qD"])
        
        # Remove data where the roll was greater than 5 degrees:
        roll_crit = 5
        highroll = abs(data_df['roll']) > roll_crit
        data_dfqc = data_df.loc[~highroll]
        
        # Split into rows of upward vs downward velocities:
        dataup = data_dfqc.loc[data_dfqc["w'"]>0]
        datadown = data_dfqc.loc[data_dfqc["w'"]<0]
        
        # Plots:
        #plt.figure()
        #plt.scatter(data_dfqc["w'"], data_dfqc["q"])
        
        #plt.figure()
        #plt.scatter(data_dfqc["w'"], data_dfqc["q'"], s=1)
        
        #plt.figure()
        #plt.hist(dataup['q'], bins=20, histtype='step')
        #plt.hist(datadown['q'], bins=20, histtype='step')
        
        #plt.figure()
        #plt.hist(dataup["w'"], bins=20, histtype='step')
        #plt.hist(datadown["w'"], bins=20, histtype='step')
    
        print("percentage upward = %0.2f" % (100*len(dataup)/len(data_dfqc)))
        print("percentage downward = %0.2f" % (100*len(datadown)/len(data_dfqc)))
        
        
        plt.figure()
        ax1 = plt.subplot(1,2,1)
        qdD_pdf(dataup, 'blue', ax1)
        qdD_pdf(datadown, 'red', ax1)
        
        
        # Same as above execpt for data above / below the upper / lower quartile:
        q1, q3 = np.nanquantile(data_dfqc["w'"], [0.25, 0.75])
        dataup = data_dfqc.loc[data_dfqc["w'"]>q3]
        datadown = data_dfqc.loc[data_dfqc["w'"]<q1]
        ax2 = plt.subplot(1,2,2)
        qdD_pdf(dataup, 'blue', ax2)
        qdD_pdf(datadown, 'red', ax2)
        
        
        
def qdD_pdf(data, color, ax, scatter=False):
        
    
        if scatter: ax.scatter(data["q"], data["dD"], s=1, c=color)
        kernel = gaussian_kde([data["q"].values, data["dD"].values])
        qq, dDdD = np.meshgrid(
            np.linspace(np.min(data["q"]), np.max(data["q"])),
            np.linspace(np.min(data["dD"]), np.max(data["dD"]))
            )
        prob = np.reshape(kernel([qq.ravel(), dDdD.ravel()]).T, qq.shape)
        ax.contour(qq, dDdD, prob, colors=color)
        



plots(4)

 
    
    
