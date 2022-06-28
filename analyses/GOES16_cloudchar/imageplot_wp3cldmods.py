# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 12:55:39 2021

@author: Dean


Plot and save GOES-16 images for each P-3 cloud module. 
Image data saved:
    - Visual reflectance
    - IR temperature
    - IR threshold mask (binary)

Timestamp of GOES image data is near the P-3 mean sampling time and lat, lon 
window is centered on mean P-3 sampling location. 
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
 

def set_imageticklabels(ax):
    """
    """
    ax.set_xticklabels(ax.get_xticks().astype(str), fontsize=6)
    ax.set_yticklabels(ax.get_yticks().astype(str), fontsize=6)    

 
def createimages(path_goescldmoddir, path_savedir):
 
    fnames_goescldmod = os.listdir(path_goescldmoddir)    
 
    for f in fnames_goescldmod:
        
        goesdata = xr.load_dataset(os.path.join(path_goescldmoddir, f)) 
        
        ncld_str = f[-9:-3]

                
        #fig_ir = plt.figure(figsize=(1.25, 1.25))
        fig_ir = plt.figure(figsize=(2, 2))
        ax_ir = fig_ir.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_ir.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['temperature_ir'], 
            cmap = 'gist_yarg', vmin=272, vmax=300
            )
        set_imageticklabels(ax_ir)
        fnameir_save = "G16V04.0.ATOMIC.PX.02K_IRtemp_%s.png" % ncld_str
        fig_ir.savefig(os.path.join(path_savedir, fnameir_save), bbox_inches='tight')        
        plt.close(fig=fig_ir)
        
        
        fig_vis = plt.figure(figsize=(1.25, 1.25))
        ax_vis = fig_vis.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_vis.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['reflectance_vis'], 
            cmap = 'gray', vmin=0, vmax=0.8
            )
        set_imageticklabels(ax_vis)
        fnamevis_save = "G16V04.0.ATOMIC.PX.02K_vis-reflectance_%s.png" % ncld_str
        fig_vis.savefig(os.path.join(path_savedir, fnamevis_save), bbox_inches='tight')
        plt.close(fig=fig_vis)
        
        
        fig_irmask = plt.figure(figsize=(1.25, 1.25))
        ax_irmask = fig_irmask.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_irmask.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['ir_mask'], cmap='gist_yarg'
            )
        set_imageticklabels(ax_irmask)
        fnameirmask_save = "G16V04.0.ATOMIC.PX.02K_IRmask_%s.png" % ncld_str
        fig_irmask.savefig(os.path.join(path_savedir, fnameirmask_save), bbox_inches='tight')   
        plt.close(fig=fig_irmask)
    


if __name__=="__main__":
    
    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)   
    createimages(path_goescldmoddir, path_savedir)

    
    
    
     
    
