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



def cospectra(data_df, varkeypairs):
    """
    Get cospectra of one or more set of variables. 
    QC data for high aircraft roll first.
    
    Inputs
    ------
    data_df: pandas.DataFrame
        Time series data for a single level leg. Has key "w'" for the vertical 
        wind perturbation.
        
    varkeypairs: list of 2-tuples of str's
        Variable key pairs in data_df to get cospectra for.

    Returns
    -------
    cospectra: pandas.DataFrame. Column are the cospectra for each pair in 
        varkeypairs.
    """

    # Aircraft roll QC    
    data_df = rollqc.meanimpute(data_df, varkeys, roll_crit=5) # Impute
    #data_df = rollqc.gaussianimpute(data_df, varkeys, roll_crit=5) # Impute


    # Cospectra:
    cospectra = pd.DataFrame({}) # Collect results here.
    for vpair in varkeypairs:
        
        vk1 = vpair[0]
        vk2 = vpair[1]
        
        fs = 5     # Sampling frequency        
        dt_window=3*60 # time interval of window for spectral decomp, seconds

        f, Pcross = signal.csd(    # Cross-spectrum
            data_df[vk1], data_df[vk2],  
            fs=fs, nperseg=int(dt_window*fs), 
            window='hamming', detrend='linear'
            )
        Pco = Pcross.real  # Co-spectrum
        f, Pco = f[1:], Pco[1:] # Remove 0 frequency.

        cospectra[vk1+vk2] = Pco
        cospectra['freq'] = f
        
        
    # Correction for imputed values:
    N_tot = len(data_df)
    N_imputed = data_df['imputed'].sum()
    for vpair in varkeypairs:
        vk1 = vpair[0]
        vk2 = vpair[1]
        cospectra[vk1+vk2] = cospectra[vk1+vk2]*(N_tot/(N_tot-N_imputed))

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
    dir_cospectra = "./cospectra_levlegs/"
    if not os.path.isdir(dir_cospectra):
        os.makedirs(dir_cospectra)
    
    
    # Keys of variables to get cospectra for:
    varkeys = ["u","v","w","T","q","qD"]
    varkeys = [k+"'" for k in varkeys]
    varkeypairs = [("u'","u'"), ("v'","v'")] + [(vk, "w'") for vk in varkeys] 
        
    
    for fname in fnames_levlegs:
    
        # Cospectra:
        data = xr.load_dataset(dir_5hzdata + fname)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        print("Computing cospectra for %s" % fname)
        cospec_df, N_tot, N_imputed = cospectra(data_df, varkeypairs)
        
        
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

        
        fnamesave = dir_cospectra + ("WP3_%s_cospectra.nc" % fname[8:-3])
        cospec_xr.to_netcdf(fnamesave)




