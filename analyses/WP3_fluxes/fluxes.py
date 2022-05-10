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

# Local code
import thermo



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
    
    print("Working on cloud module # %s" % n)
    
    fnames_cldmod = [f for f in fnames_wcospec if ("_cld%s" % n) in f]
    
    
    # Objects to store results in before combining into a dataset:
    covarkeys = ["u'w'", "v'w'", "w'w'", "T'w'", "q'w'", "qD'w'"]
    fluxkeys = ["%s_bar" % k for k in covarkeys]
    fluxes = dict(zip(fluxkeys, [[] for e in fluxkeys]))
    frac_imputed = [] # Fraction of imputed data points for each leg.
    addvar_keys = ['T_mean', 'P_mean', 'reftime'] # Additional variables
    addvars = dict(zip(addvar_keys, [[] for e in addvar_keys]))
    coord_keys = ['n_levleg', 'alt_mean'] # Coordinates for the xr.Dataset
    coords = dict(zip(coord_keys, [[] for e in coord_keys]))
    

    # Results for each leg:
    for fn in fnames_cldmod:
        wcospec = xr.load_dataset(dir_wcospec + fn)
        
        # Fluxes from integral of power spectra:
        for vk, fk in zip(covarkeys, fluxkeys):
            flux_vk = trapz(wcospec[vk], x=wcospec['freq'])            
            # Correct for imputed data points before appending:
            #N_tot = wcospec.attrs['Ntot_ts']
            #N_imputed = wcospec.attrs['Nimputed_ts']
            #flux_vk = flux_vk*(N_tot/(N_tot-N_imputed))
            
            fluxes[fk].append(flux_vk)

        # Additional vars:
        frac_imputed.append(wcospec.attrs['Nimputed_ts']/wcospec.attrs['Ntot_ts'])
        for k in addvar_keys:
            addvars[k].append(wcospec[k].values.item())
            
        # Info that will be used as coordinates:
        for k in coord_keys:
            coords[k].append(wcospec[k].values.item())        
        
        
    # Combine into xarray.Dataset:
    dsvars = {**fluxes, **addvars}
    flux_ds = xr.Dataset(
        data_vars=dict(
            zip(
                [k for k in dsvars], 
                [(["n_levleg"], dsvars[k]) for k in dsvars]
                )
            ),
        coords=dict(
            n_levleg=coords['n_levleg'],
            alt=(['n_levleg'], coords['alt_mean'])
            ),
        attrs=dict(
            description="Fluxes from P-3 level legs during ATOMIC, "
                "computed from 5Hz data."
            )
        )
    
    # Also add percentage of imputed data points:
    flux_ds = flux_ds.assign(
        frac_imputed = xr.DataArray(
            data=frac_imputed,
            dims=["n_levleg"],
            coords=dict(n_levleg=coords['n_levleg']),
            )
        )
           
    
    # Also include sensible and latent heat fluxes:
    flux_ds['flux_sh'] = thermo.flux_sh(
        flux_ds["T'w'_bar"], flux_ds["T_mean"], flux_ds["P_mean"]
        )
    flux_ds['flux_lh'] = thermo.flux_lh(
        flux_ds["q'w'_bar"]/1000, flux_ds["T_mean"], flux_ds["P_mean"]
        )

    
    # Save:
    date = fnames_cldmod[0][4:12]
    flux_ds.to_netcdf(dir_fluxes + "WP3_%s_cld%s_fluxes.nc" % tuple([date, n]))

