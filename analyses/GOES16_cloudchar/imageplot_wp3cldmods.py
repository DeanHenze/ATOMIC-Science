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
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

    
    
## I/O paths
##_____________________________________________________________________________
path_cldmodtab = "../WP3_cloudmodule_char/p3_cloudmodules.csv"
path_goescldmoddir = "./goes16_WP3cldmods/"
path_savedir = "./images/"
path_levlegsdir = "../WP3_fluxes/levlegdata_5hz/"
##_____________________________________________________________________________
## I/O paths   
 

def set_imageticklabels(ax):
    """
    """
    ax.set_xticklabels(ax.get_xticks().astype(str), fontsize=6)
    ax.set_yticklabels(ax.get_yticks().astype(str), fontsize=6)  

 
def createimages(path_goescldmoddir, path_levlegsdir, path_savedir):
 
    ncld = [str(n).zfill(2) for n in np.arange(1, 17, 1)] # Cloud module numbers
    # goes and level leg filenames:
    fnames_goescldmod = os.listdir(path_goescldmoddir)
    fnames_levlegs = [f for f in os.listdir(path_levlegsdir) if f.endswith(".nc")]
 
    for n in ncld:
        
        # GOES filename for this cloud module - should only be one:
        f_goes = [f for f in fnames_goescldmod if "_ncld%s" % n in f][0]
        # P-3 level leg filenames - multiple
        f_levlegs = [f for f in fnames_levlegs if "_cld%s" % n in f]
        
        
        # Load GOES data and plot image:
        goesdata = xr.load_dataset(os.path.join(path_goescldmoddir, f_goes))
                
        latmin = goesdata['latitude'].min()
        latmax = goesdata['latitude'].max()
        lonmin = goesdata['longitude'].min()
        lonmax = goesdata['longitude'].max()
        
        fig_ir = plt.figure(figsize=(2, 2))
        ax_ir = fig_ir.add_axes([0.2, 0.2, 0.7, 0.7])
        axim = ax_ir.imshow(
            goesdata['temperature_ir'], 
            cmap = 'gist_yarg', vmin=272, vmax=300, 
            extent =[lonmin, lonmax, latmin, latmax],
            interpolation ='nearest', origin ='upper'
            )
        #ax_ir.pcolor(
        #    goesdata['longitude'], goesdata['latitude'], 
        #    goesdata['temperature_ir'], 
        #    cmap = 'gist_yarg', vmin=272, vmax=300, 
        #    extent =[lonmin, lonmax, latmin, latmax],
        #             interpolation ='nearest', origin ='lower'
        #    )
        #ax_ir.set_xlim(lonmin, lonmax)
        #ax_ir.set_ylim(latmin, latmax)
        #set_imageticklabels(ax_ir)
        ax_ir.axis('off')
        ax_ir.set_title('cld %s' % n, fontsize=10)
        


        # Overlay P-3 level leg:
        for f in f_levlegs[1:-1]:
            p3data = xr.load_dataset(os.path.join(path_levlegsdir, f))
            ax_ir.plot(p3data['lon'], p3data['lat'], color='red', linewidth=0.5)


        # For the first cloud module, draw a 50 km line for reference:
        if n=='01':
            lons_refline = (lonmax-1.-0.85, lonmax-1.) # 0.85 deg. ~50km at 15N lat.
            lats_refline = (latmin+0.25, latmin+0.25)
            ax_ir.plot(lons_refline, lats_refline, color='orange', linewidth=3)
            ax_ir.text(
                np.mean(lons_refline), lats_refline[0], '50 km', 
                ha='center', va='bottom', fontsize=8.5, color='orange'
                )
        

        # Save image:
        fnameir_save = "G16V04.0.ATOMIC.PX.02K_IRtemp_%s.png" % n
        fig_ir.savefig(os.path.join(path_savedir, fnameir_save), bbox_inches='tight')        
        plt.close(fig=fig_ir)


        """
        fig_vis = plt.figure(figsize=(1.25, 1.25))
        ax_vis = fig_vis.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_vis.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['reflectance_vis'], 
            cmap = 'gray', vmin=0, vmax=0.8
            )
        set_imageticklabels(ax_vis)
        fnamevis_save = "G16V04.0.ATOMIC.PX.02K_vis-reflectance_%s.png" % n
        fig_vis.savefig(os.path.join(path_savedir, fnamevis_save), bbox_inches='tight')
        plt.close(fig=fig_vis)
        
        
        fig_irmask = plt.figure(figsize=(1.25, 1.25))
        ax_irmask = fig_irmask.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_irmask.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['ir_mask'], cmap='gist_yarg'
            )
        set_imageticklabels(ax_irmask)
        fnameirmask_save = "G16V04.0.ATOMIC.PX.02K_IRmask_%s.png" % n
        fig_irmask.savefig(os.path.join(path_savedir, fnameirmask_save), bbox_inches='tight')   
        plt.close(fig=fig_irmask)
        
        
        fig_ctt = plt.figure(figsize=(1.25, 1.25))
        ax_ctt = fig_ctt.add_axes([0.2, 0.2, 0.7, 0.7])
        ax_ctt.pcolor(
            goesdata['longitude'], goesdata['latitude'], 
            goesdata['cloud_top_temperature'], 
            cmap = 'gray', vmin=270, vmax=292
            )
        set_imageticklabels(ax_ctt)
        fnamectt_save = "G16V04.0.ATOMIC.PX.02K_cloudtop-temperature_%s.png" % n
        fig_ctt.savefig(os.path.join(path_savedir, fnamectt_save), bbox_inches='tight')
        plt.close(fig=fig_ctt)
        """


if __name__=="__main__":
    
    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)   
    createimages(path_goescldmoddir, path_levlegsdir, path_savedir)

    
    
    
     
    
