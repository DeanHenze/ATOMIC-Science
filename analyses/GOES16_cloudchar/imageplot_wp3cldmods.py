# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:55:39 2021

@author: Dean


Produces one representative GOES-16 images for each P-3 cloud module. 
The GOES images are taken at the nearest 20 minute timestamp to the mean 
time of the P-3 cloud modules. Plot either reflectance for daytime or IR 
temperature for night (in progress). Overlays P-3 flight track during cloud 
module.
"""


# Built in:
import sys
import os
#from os import listdir
#from os.path import join

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.legend import Legend
import matplotlib.lines as mlines
import seaborn as sns

# My modules:
if r'../../' not in sys.path: sys.path.insert(0,r'../../')
from henze_python_modules import atomic_data_loader as adl
from henze_python_modules import iso_fxns



def plt_goes16(goesdata, p3data_cld, ax, radtype='vis'):
    """ pcolor plot of GOES-16 visible reflectance or IR. Overlay P-3 flight 
    track. radtype can be 'vis' or 'ir'.
    """
    latmean = p3data_cld['lat'].mean()
    lonmean = p3data_cld['lon'].mean()
    goesdata = goesdata.where((goesdata['longitude']>(lonmean-2.5)) 
                              & (goesdata['longitude']<(lonmean+2.5)), 
                              drop=True)
    goesdata = goesdata.where((goesdata['latitude']>(latmean-2.5)) 
                              & (goesdata['latitude']<(latmean+2.5)), 
                              drop=True)
    
    if radtype=='vis': 
        varkey = 'reflectance_vis'
        vmin = 0; vmax = 0.8
        cmap = 'gray'
    if radtype=='ir': 
        varkey = 'temperature_ir'
        vmin = 272; vmax = 300
        cmap = 'gist_yarg'
   
    pc = ax.pcolor(goesdata['longitude'], goesdata['latitude'], goesdata[varkey], 
                   cmap = cmap, vmin=vmin, vmax=vmax)
    ax.plot(p3data_cld['lon'], p3data_cld['lat'], 'r-')
    
    return pc



def get_goes(t1_cld, t2_cld):
    """
    Get a single GOES .nc data file for the mean time of t1_cld and t2_cld
    (pd.Timestamp objects), to the nearest 20 minutes.

    Returns
    -------
    xarray.dataset

    """    
    # Make sure times are pd.Timestamps:
    t1_cld = pd.Timestamp(t1_cld)
    t2_cld = pd.Timestamp(t2_cld)
    
    # Mean datetime of cloud module rounded to nearest 20 min:
    dtime_cld = (t1_cld + (t2_cld - t1_cld)/2).round('min')
    dtime_cld += pd.Timedelta(10, 'min') # Add 10 min so next line floors instead of rounds.
    dtime_cld -= pd.Timedelta(dtime_cld.minute % 20, 'min') # nearest 20 min.
    
    # Load corresponding GOES-16 file and isolate a few variables:
    goes16 = adl.goes16_p3cld(dtime_cld)
    goes16 = goes16[['reflectance_vis', 'temperature_ir', 'cloud_top_height']]
    
    return goes16, dtime_cld



###############################################################################
# Rest of script loops through cloud modules datetimes and calls above fxns 
# make make figures.
###############################################################################


## P-3 data filenames each cloud module:
datadir = r"../output/"
fnames = [os.path.join(datadir, f) for f in os.listdir(datadir) 
          if f.startswith("p3cld_") and f.endswith(".nc")
          ]

for f in fnames:
    
    # Load P-3 data during cloud module:
    p3data = xr.load_dataset(f) 
    
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])

    # Get GOES data and plot image with P-3 flight track:
    t1_cld = p3data.time.values[0]
    t2_cld = p3data.time.values[-1]
    goes, dtime_cld_20min = get_goes(t1_cld, t2_cld)
    radtype = 'ir'
    pc = plt_goes16(goes, p3data, ax, radtype=radtype)
    
    #fig.colorbar(pc, ax=ax)
    
    # Figure aesthetics:
    ax.tick_params(axis='both', labelsize=8.5)
    
    date = f[-14:-6]; n_cld = f[-5:-3] # Date and cloud module number
    ax.set_title(str(dtime_cld_20min), fontsize=12)
    
    
    figsdir = r"../figures/goes16/"
    if not os.path.isdir(figsdir):
        os.makedirs(figsdir)
    fname_fig = "goes16%s_%s_cld%s" % tuple([radtype, date, n_cld])
    fig.savefig(figsdir + fname_fig)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    
