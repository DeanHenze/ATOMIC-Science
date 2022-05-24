# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:53:38 2022

@author: Dean

In-situ data product for the P-3 cloud modules.


Current status: 
Only the iso the flightlevel datasets are merged currently.
The iso the flightlevel datasets have compatible time 
dimensions. The remote dataset is not at 1 Hz frequency (> 1 Hz?). The 
microphysics dataset is only for a subset of the flights.

###############
Add as attribute 'platform' = which dataset I got each variable from
#################
"""



# Built in:
import os

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt



def create_files(path_flightlev_dir, path_iso_dir, path_remote_dir, 
                 path_cldmodtable, path_hlegstable, dirpath_save):
    """
    Create a separate data file for each cloud module.
    
    Inputs
    -------
    path_drpsnds: str/path.
        Path (rel or abs) to table of P-3 cloud module times (.csv file).
    
    path_cldmodtable: str/path.
        Path (rel or abs) to dropsondes .nc file.

    dirpath_save: str/path.
        Path (rel or abs) to directory in which to save output.
    """
        
    # Load cloud module table and dropsonde dataset:
    tab_cldmod = table_cldmodules(path_cldmodtable, datetime_type=True)
    tab_hlegs = pd.read_csv(path_hlegstable)
    
    
    for i, row in tab_cldmod.iterrows():
        print(row['flight_date'])
        
        date = row['flight_date']
        ncld = row['num_cld_iop']
        
        
        # Load P-3 data sets:
        fname_flightlev = "EUREC4A_ATOMIC_P3_Flight-level_%i_v1.0.nc" % date
        flightlev = xr.load_dataset(path_flightlev_dir + fname_flightlev)
        fname_iso = "EUREC4A_ATOMIC_P3_IsotopeAnalyzer_%i_v0.0.nc" % date
        iso = xr.load_dataset(path_iso_dir + fname_iso)
        fname_remote = "EUREC4A_ATOMIC_P3_Remote-sensing_%i_v1.1.nc" % date
        remote = xr.load_dataset(path_remote_dir + fname_remote)


        insituremote_cld = cld_single(
            flightlev, iso, remote, 
            date, ncld, row['start_datetime'], row['end_datetime'], 
            tab_hlegs
            )
            
            #drpsnds, date, 
            #row['start_datetime'], row['end_datetime'], 
            #)
        
        # Save:
        ncld_str = str(int(ncld)).zfill(2)
        fname = r"p3cld_insitu+remote_%i_%s.nc" % tuple([date, ncld_str])
        insituremote_cld.to_netcdf(dirpath_save + fname)
        
        
        verification_plot(insituremote_cld)
        #fname = "/p3cld_dropsondes_%i_ncld%s.nc" % tuple([date, ncld_str])
        #insituremote_cld.to_netcdf(dirpath_save + fname)
        
        

def cld_single(flightlev, iso, remote, date, ncld, t1_cld, t2_cld, tab_hlegs):
    """
    date: int
    ncld: scalar. Cloud module number for entire IOP.
    t1_cld, t2_cld: Timestamps
    tab_hlegs: pandas.DataFrame
    """
    
    ## Load P-3 datasets for the flight date:
    #mphys = adl.p3_microphys(date) # Microphysics
    #path_flightlev = ("../../data/WP3/flight_level_1Hz/"
    #                  "EUREC4A_ATOMIC_P3_Flight-level_%i_v1.0.nc" % date)
    #flightlev = xr.load_dataset(path_flightlev) # Flight level data (T, P, etc.)
    #path_iso = ("../../data/WP3/water_iso_1Hz/"
    #            "EUREC4A_ATOMIC_P3_IsotopeAnalyzer_%i_v0.0.nc" % date)
    #iso = xr.load_dataset(path_iso) # Isotope ratios
    #path_remote = ("../../data/WP3/remote_sensing/"
    #               "EUREC4A_ATOMIC_P3_Remote-sensing_%i_v1.1.nc" % date)
    #remote = xr.load_dataset(path_remote) # Remote sensing products.


    ## The remote products dataset in not 1 Hz so fix this:
    t_rem = remote['time'].values
    tref_1Hz = np.arange(t_rem[0], t_rem[-1], 
                         dtype='datetime64[s]') # Reference time array @ 1 Hz.
    remote_fixed = remote.interp(time=tref_1Hz, method='nearest')

    
    cld = xr.merge(
        [
            iso[['mr', 'dD', 'rh', 'lat', 'lon', 'alt']], 
            flightlev[['Ta', 'press', 'wvel', 'RH']], 
            remote_fixed[['alt_CT', 'SST_IR_est']]
            ], 
        join='left'
        )           
    cld = cld.sel(time = slice(t1_cld, t2_cld))
    
    
    ## Add level leg flag as a variable:
        # Empty data array:
    levlegflag = xr.DataArray(
        data = np.zeros(len(cld['time'])), 
        dims = ['time'], 
        coords = dict(time=cld['time']), 
        attrs = dict(long_name = "Integer flag for horizontal leg number "
                         "of the cloud module.", 
                     units = "0 if not in a leg, otherwise an integer.")
        )    
        # Row in the table for this cloud module:        
    hlegs_cld = tab_hlegs.loc[(tab_hlegs['flight_date']==date)
                              & (tab_hlegs['num_cld_iop']==ncld)]
        # Fill in flag with integers corresponding to the level leg intervals:
    for i, row in hlegs_cld.reset_index().iterrows():
        levlegflag.loc[row['tstart_leg']: row['tend_leg']] = i+1
    
    cld = cld.assign(levlegflag = levlegflag)
    
    
    ## Include horizontal cloud type as a dataset attribute:
    cld.attrs['cloud_type'] = ""

    
    return cld



def table_cldmodules(path_cldmodtable, datetime_type=True):
    """
    Returns .csv file of P-3 cloud module time interval data.
    
    datetime_type: bool.
        Default = True, in which case all datetime columns are converted to 
        datetime-type.
    """
    cldmods = pd.read_csv(path_cldmodtable)
    
    if datetime_type:
        dtimekeys = ['start_datetime', 'end_datetime']
        cldmods[dtimekeys] = cldmods[dtimekeys].apply(pd.to_datetime)
    
    return cldmods



def verification_plot(cld):
    """
    cld: xarray.Dataset
    """
    fig = plt.figure()
    plt.plot(cld.time, cld.alt, color='black')
    for i in np.unique(cld['levlegflag'])[1:]:
        temp = cld.where(cld['levlegflag'] == i, drop=True)
        plt.plot(temp.time, temp.alt)
        
    return fig
   


###############################################################################
# Remainder of script creates a separate data file for each cloud module.
###############################################################################
"""
# Load cloud module info tables:
tab_cld = pd.read_csv('p3_cloudmodules.csv')
tab_hlegs = pd.read_csv('p3_cloudmodule_legs.csv')

for i, row in tab_cld.iterrows():
    print(row['flight_date'])
    
    date = row['flight_date']
    ncld = row['num_cld_iop']
    
    data_cld = cld_single(
        date, ncld, 
        row['start_datetime'], row['end_datetime'], 
        tab_hlegs # Needed for the level leg flag.
        )
    
    # Add cloud classification as a dataset attribute:
    data_cld.attrs['cloud_type'] = row['cloud_type']
    
    # Save:
    ncld_str = str(int(ncld)).zfill(2)
    fname = r"./cldmod_datafiles/p3cld_insitu+remote_%i_%s.nc" % tuple([date, ncld_str])
    data_cld.to_netcdf(fname)
    
    #fig = verification_plot(data_cld)
    #fig.savefig(r"../scratch/levleg_altcheck_%i_%s_v2.png" % tuple([date, ncld_int]))
"""



if __name__=="__main__":
    """
    Specify paths / directories and call fxn's to get all cloud module 
    dropsonde files.
    """
    
    # P-3 data file directory paths:
    path_flightlev_dir = ("../../data/WP3/flight_level_1Hz/")
    path_iso_dir = ("../../data/WP3/water_iso_1Hz/")
    path_remote_dir = ("../../data/WP3/remote_sensing/")
    
    # P-3 cloud module table paths:
    path_cldmodtable = "./p3_cloudmodules.csv"
    path_hlegstable = "./p3_cloudmodule_legs.csv"
    
    dirpath_save = "./cldmod_datafiles/"
    if not os.path.isdir(dirpath_save):
        os.mkdir(dirpath_save)
    
    create_files(path_flightlev_dir, path_iso_dir, path_remote_dir, 
                 path_cldmodtable, path_hlegstable, dirpath_save)


