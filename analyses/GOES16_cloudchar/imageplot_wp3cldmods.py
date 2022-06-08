# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:55:39 2021

@author: Dean


Plot and save one representative GOES-16 image for each P-3 cloud module. GOES 
image timestamp is near the P-3 mean sampling time and lat, lon window is 
centered on mean P-3 sampling location. 
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
    

 
def createimages(path_goescldmoddir, path_savedir, varkey='temperature_ir'):  
 
    fnames_goescldmod = os.listdir(path_goescldmoddir)    
 
    for f in fnames_goescldmod:
        
        goesdata = xr.load_dataset(os.path.join(path_goescldmoddir, f)) 
        
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_axes([0.2, 0.2, 0.7, 0.7])
        
        if varkey=='reflectance_vis': 
            vmin = 0; vmax = 0.8
            cmap = 'gray'
        if varkey=='temperature_ir': 
            vmin = 272; vmax = 300
            cmap = 'gist_yarg'
 
        pc = ax.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata[varkey], 
            cmap = cmap, vmin=vmin, vmax=vmax
            )
        
        ncld_str = f[-9:-3]
        fnamesave = "G16V04.0.ATOMIC.PX.02K_%s_%s.png" % tuple([ncld_str, varkey])
        fig.savefig(os.path.join(path_savedir, fnamesave))
            


if __name__=="__main__":
    
    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)   
    createimages(path_goescldmoddir, path_savedir, varkey='temperature_ir')

    
    
    
     
    
