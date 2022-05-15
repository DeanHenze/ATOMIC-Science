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
        ax = plt.axes()
        qdD_pdf(dataup, 'blue', ax)
        qdD_pdf(datadown, 'red', ax)
        
        
        # Same as above execpt for data above / below the upper / lower quartile:
        q1, q3 = np.nanquantile(data_df["w'"], [0.25, 0.75])
        dataup = data_dfqc.loc[data_dfqc["w'"]>q3]
        datadown = data_dfqc.loc[data_dfqc["w'"]<q1]
        
        
        
def qdD_pdf(data, color, ax):
        
        ax.scatter(data["q"], data["dD"], s=1, c=color)
        kernel = gaussian_kde([data["q"].values, data["dD"].values])
        qq, dDdD = np.meshgrid(
            np.linspace(np.min(data["q"]), np.max(data["q"])),
            np.linspace(np.min(data["dD"]), np.max(data["dD"]))
            )
        prob = np.reshape(kernel([qq.ravel(), dDdD.ravel()]).T, qq.shape)
        ax.contour(qq, dDdD, prob, colors=color)
        



plots(4)

"""
for i in (8, 9, 10, 11):
    fname = fnames_levlegs[i]
    print(fname)
    data = xr.load_dataset(dir_5hzdata + fname)
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

    plt.scatter(dataup["q"], dataup["dD"], s=1, c='blue')
    kernelup = gaussian_kde([dataup["q"].values, dataup["dD"].values])
    qq, dDdD = np.meshgrid(
        np.linspace(np.min(data_dfqc["q"]), np.max(data_dfqc["q"])),
        np.linspace(np.min(data_dfqc["dD"]), np.max(data_dfqc["dD"]))
        )
    probup = np.reshape(kernelup([qq.ravel(), dDdD.ravel()]).T, qq.shape)
    plt.contour(qq, dDdD, probup, colors='blue')
    
    plt.scatter(datadown["q"], datadown["dD"], s=1, c='red')
    kerneldown = gaussian_kde([datadown["q"].values, datadown["dD"].values])
    probdown = np.reshape(kerneldown([qq.ravel(), dDdD.ravel()]).T, qq.shape)
    plt.contour(qq, dDdD, probdown, colors='red')    
"""    
    
    
