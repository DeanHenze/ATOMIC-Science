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

        cospectra[vk] = Pco
        cospectra['freq'] = f

    return cospectra


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

    data = xr.load_dataset(dir_5hzdata + fname)
    data_df = data.to_dataframe() # Easier / faster to work with pandas.
    print("Computing cospectra for %s" % fname)
    wcospec = wcospectra(data_df, varkeys)
    
    # Convert back to xarray and save:
    fnamesave = dir_cospectra + ("WP3_%s_wcospectra.nc" % fname[8:-3])
    wcospec.to_xarray().to_netcdf(fnamesave)




