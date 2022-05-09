# -*- coding: utf-8 -*-
"""
Created on Sun May  8 16:08:31 2022

@author: Dean

Compute and save cospectra and fluxes for each P-3 level leg.
"""


# Built in
import os

# Third party
import pandas as pd
import xarray as xr
from scipy import signal

# Local code
import rollqc
import thermo



def wcospectra(data_df, varkeys):
    """
    Get cospectra of one or more variables with vertical wind perturbation. 
    
    QC data for high aircraft roll first.
    
    Inputs
    ------
    data_df: pandas.DataFrame
        Time series data for a single level leg. Has key "w'" for the vertical 
        wind perturbation.
        
    varkeys: list of str's
        Variable keys in data_df to get cospectra for.

    Returns
    -------
    cospectra: pandas.DataFrame
    """

    # Aircraft roll QC    
    data_df = rollqc.meanimpute(data_df, varkeys, roll_crit=5) # Impute

    # Cospectra of variables with vertical wind:
    cospectra = pd.DataFrame({}) # Collect results here.
    for vk in varkeys:
        
        fs = 5     # Sampling frequency        
        dt_window=3*60 # time interval of window for spectral decomp, seconds

        f, Pcross = signal.csd(    # Cross-spectrum
            data_df[vk], data_df["w'"],  
            fs=fs, nperseg=int(dt_window*fs), 
            window='hamming', detrend='linear'
            )
        Pco = Pcross.real  # Co-spectrum
        f, Pco = f[1:], Pco[1:] # Remove 0 frequency.

        cospectra[vk+"w'"] = Pco
        cospectra['freq'] = f
        
    # Correction for imputed values:
    N_tot = len(data_df)
    N_imputed = data_df['imputed'].sum()
    for vk in varkeys:
        cospectra[vk+"w'"] = cospectra[vk+"w'"]*((N_tot-N_imputed)/N_tot)

    return cospectra, N_tot, N_imputed



if __name__=="__main__":
    """
    Compute cospectra for all level leg files and save results to new 
    folder. One cospectra file per level leg.
    """

    # All level leg data filenames:
    dir_5hzdata = "./levlegdata_5hz/"
    fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]
    
    
    # Save cospectra to this directory:
    dir_cospectra = "./wcospectra_levlegs/"
    if not os.path.isdir(dir_cospectra):
        os.makedirs(dir_cospectra)
    
    
    # Keys of variables to get cospectra / fluxes for:
    varkeys = ["u","v","w","T","q","qD"]
    varkeys = [k+"'" for k in varkeys]        
        
    
    for fname in fnames_levlegs:
    
        # Cospectra:
        data = xr.load_dataset(dir_5hzdata + fname)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        print("Computing cospectra for %s" % fname)
        wcospec_df, N_tot, N_imputed = wcospectra(data_df, varkeys)
        
        
        # Additional info to include in the data file:
            # Level leg number for this cloud module:
        i_levleg = fname.index('_levleg')+7
        n_levleg = int(fname[i_levleg]) # Works b/c there are < 10 lev legs.
            # Mean altitude, temperature, pressure (estimated):
        alt_mean = data['alt'].mean().round(decimals=0).values
        T_mean = data['T'].mean().round(decimals=0).values 
        P_mean = thermo.P_est(alt_mean/1000)
        reftime = data['time'].mean().values
        
        
        # Convert back to xarray and add some vars/coords/attributes 
        # before saving:
        wcospec_df.set_index('freq', drop=True, inplace=True, append=False)
        wcospec_xr = wcospec_df.to_xarray()
            
            # Coords:        
        wcospec_xr = wcospec_xr.assign_coords({'n_levleg':n_levleg})
        
            # Vars:
        addvar_keys = ['alt_mean', 'T_mean', 'P_mean', 'reftime']
        addvars = [alt_mean, T_mean, P_mean, reftime]
        for k, v in zip(addvar_keys, addvars):
            da = xr.DataArray(
                data=[v],
                dims=["n_levleg"],
                coords=dict(n_levleg=[n_levleg]),
                )
            wcospec_xr = wcospec_xr.assign({k:da})
        
            # Attributes:
        wcospec_xr.attrs = dict(
            title="Cospectra with vertical wind computed with 5Hz data from "
                "P-3 level legs.", 
            Ntot_ts=N_tot,
            Nimputed_ts=N_imputed
            )

        
        fnamesave = dir_cospectra + ("WP3_%s_wcospectra.nc" % fname[8:-3])
        wcospec_xr.to_netcdf(fnamesave)




