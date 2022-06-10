# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 13:53:05 2022

@author: Dean
"""



# Built in
import os

# Third party
import xarray as xr
import matplotlib.pyplot as plt

    
    
## I/O paths
##_____________________________________________________________________________
path_cldmodtab = "../WP3_cloudmodule_char/p3_cloudmodules.csv"
path_goescldmoddir = "./goes16_WP3cldmods/"
path_savedir = "./images/"
##_____________________________________________________________________________
## I/O paths



fnames_goescldmod = os.listdir(path_goescldmoddir)    
 
for f in fnames_goescldmod:
    
    goesdata = xr.load_dataset(os.path.join(path_goescldmoddir, f)) 
    ncld_str = f[-9:-3]
    
    n_cloudy = goesdata['ir_mask'].sum().item()
    n_total = goesdata['ir_mask'].size
    cc = n_cloudy/n_total
    print("CC for cloud module %s = %0.2f" % tuple([ncld_str, cc]))

    
    
    
    
    
    
    