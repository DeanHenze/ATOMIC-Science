# -*- coding: utf-8 -*-
"""
Created on Sun May  8 16:08:31 2022

@author: Dean

Extension of the 'cospectra.py' script for varying window sizes. Cospectra 
for each leg saved with different window sizes as different files.
"""


# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal

# Local code
import rollqc
import thermo
import cospectra



def cospectrafiles(dt_window):
    """
    Compute cospectra with a specified window size (seconds) for all level leg 
    files and save results to new folder. One cospectra file per level leg.
    """
    # All level leg data filenames:
    dir_5hzdata = "./levlegdata_5hz/"
    fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]
    
    
    # Save cospectra to this directory:
    dir_cospectra = "./cospectra_levlegs_windowvary/"
    if not os.path.isdir(dir_cospectra):
        os.makedirs(dir_cospectra)
    
    
    # Keys of variables to get cospectra for:
    varkeys = ["u","v","w","T","q","qD"]
    varkeys = [k+"'" for k in varkeys]
    #varkeypairs = [("u'","u'"), ("v'","v'")] + [(vk, "w'") for vk in varkeys] 
    varkeypairs = [(vk, vk) for vk in varkeys] + [(vk, "w'") for vk in varkeys] 
        
    
    for fname in fnames_levlegs:
    
        # Cospectra:
        data = xr.load_dataset(dir_5hzdata + fname)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        print("Computing cospectra for %s" % fname)
        cospec_df, N_tot, N_imputed = \
            cospectra.cospectra(data_df, varkeypairs, dt_window=dt_window)
        
        
        # Additional info to include in the data file:
            # Level leg number for this cloud module:
        i_levleg = fname.index('_levleg')+7
        n_levleg = int(fname[i_levleg]) # Works b/c there are < 10 lev legs.
            # Mean altitude, temperature, pressure (estimated):
        alt_mean = data['alt'].mean().round(decimals=0).values
        T_mean = data['T'].mean().round(decimals=0).values 
        P_mean = thermo.P_est(alt_mean/1000)
        q_mean = data['q'].mean().round(decimals=3).values
        reftime = data['time'].mean().values
        
        
        # Convert back to xarray and add some vars/coords/attributes 
        # before saving:
        cospec_df.set_index('freq', drop=True, inplace=True, append=False)
        cospec_xr = cospec_df.to_xarray()
            
            # Coords:        
        cospec_xr = cospec_xr.assign_coords({'n_levleg':n_levleg})
        
            # Vars:
        addvar_keys = ['alt_mean', 'T_mean', 'P_mean', 'q_mean', 'reftime']
        addvars = [alt_mean, T_mean, P_mean, q_mean, reftime]
        for k, v in zip(addvar_keys, addvars):
            da = xr.DataArray(
                data=[v],
                dims=["n_levleg"],
                coords=dict(n_levleg=[n_levleg]),
                )
            cospec_xr = cospec_xr.assign({k:da})
        
            # Attributes:
        cospec_xr.attrs = dict(
            title="Cospectra with vertical wind computed with 5Hz data from "
                "P-3 level legs.", 
            Ntot_ts=N_tot,
            Nimputed_ts=N_imputed
            )

        
        fnamesave =  ("WP3_%s_cospectra_window%ss.nc" 
                      % tuple([fname[8:-3], str(dt_window).zfill(3)]))
        cospec_xr.to_netcdf(dir_cospectra + fnamesave)




if __name__=="__main__":
    
    #dt_window_list = [30, 1*60, 2*60, 3*60, 4*60, 5*60]
    dt_windows = np.arange(30, 5*60, 15)
    for dt_window in dt_windows: cospectrafiles(dt_window)


