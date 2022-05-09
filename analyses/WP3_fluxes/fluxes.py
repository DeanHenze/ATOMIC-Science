# -*- coding: utf-8 -*-
"""
Created on Sun May  8 17:28:45 2022

@author: Dean

Fluxes for P-3 level legs, computed from cospectra. 
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
from scipy import signal
from scipy.integrate import trapz



# Save fluxes to this directory:
dir_fluxes = "./fluxes_levlegs/"
if not os.path.isdir(dir_fluxes):
    os.makedirs(dir_fluxes)


# All cospectra filenames:
dir_wcospec = "./wcospectra_levlegs/"
fnames_wcospec = [f for f in os.listdir(dir_wcospec) if f.endswith(".nc")]


# All unique cloud module number cases:
i1_list = [f.index('_cld') + 4 for f in fnames_wcospec]
i2_list = [f.index('_levleg') for f in fnames_wcospec]
ncldmod = [f[i1:i2] for f, i1, i2 
           in zip(fnames_wcospec, i1_list, i2_list)
          ]
ncldmod_unique = np.unique(ncldmod)


for n in ncldmod_unique:
    
    fnames_cldmod = [f for f in fnames_wcospec if ("_cld%s" % n) in f]
    
    
    # Objects to store results in before combining into a dataset:
    covarkeys = ["u'w'", "v'w'", "w'w'", "T'w'", "q'w'", "qD'w'"]
    fluxes = dict(zip(covarkeys, [[] for e in covarkeys]))
    frac_imputed = [] # Fraction of imputed data points for each leg.
    altmean = [] # Mean altitude of each leg.
    n_levleg = [] # Level leg number for the cloud module.


    # Results for each leg:
    for fn in fnames_cldmod:
        wcospec = xr.load_dataset(dir_wcospec + fn)
        
        # Fluxes from integral of power spectra:
        for vk in covarkeys:
            flux_vk = trapz(wcospec[vk], x=wcospec['freq'])
            fluxes[vk].append(flux_vk)

        # Other info:
        frac_imputed.append(wcospec.attrs['Nimputed_ts']/wcospec.attrs['Ntot_ts'])
        altmean.append(wcospec.attrs['alt_mean'])
        i_levleg = fn.index('_levleg')+7
        n_levleg.append(int(fn[i_levleg])) # Works b/c there are < 10 lev legs.
        
        
    # Combine into xarray.Dataset:
    ds = xr.Dataset(
        data_vars=dict(
            zip(
                ["flux_"+vk for vk in covarkeys], 
                [(["alt"], fluxes[flux_vk]) for flux_vk in fluxes]
                )
            ),
        coords=dict(
            alt=altmean,
            n_levleg=(['alt'], n_levleg)
            #reference_time=reference_time,
            ),
        attrs=dict(
            description="Fluxes from P-3 level legs during ATOMIC, "
                "computed from 5Hz data."
            )
        )
    
    
    # Also include sensible and latent heat fluxes:
    

        

    # Correction for number of imputed data points:
        #...

