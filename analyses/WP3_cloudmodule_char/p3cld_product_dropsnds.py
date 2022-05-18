# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 10:53:38 2022

@author: Dean

Dropsonde data product for the P-3 cloud modules. Saves one file of 
dropsonde data per cloud module.
"""



# Built in
import os

# Third party
import pandas as pd
import xarray as xr



def create_dropsondefiles(path_drpsnds, path_cldmodtable, dirpath_save):
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
    drpsnds = xr.load_dataset(path_drpsnds)
    
    for i, row in tab_cldmod.iterrows():
        print(row['flight_date'])
        
        date = row['flight_date']
        ncld = row['num_cld_iop']
        
        dropsnds_cld = dropsnds_cld_single(
            drpsnds, date, ncld, 
            row['start_datetime'], row['end_datetime'], 
            )
        
        # Save:
        ncld_str = str(int(ncld)).zfill(2)
        fname = "/p3cld_dropsondes_%i_ncld%s.nc" % tuple([date, ncld_str])
        dropsnds_cld.to_netcdf(dirpath_save + fname)
    
  

def dropsnds_cld_single(dropsnds, date, ncld, t1_cld, t2_cld):
    """
    dropsnds: xarray.Dataset
    date: int
    ncld: int, float
    t1_cld, t2_cld: Timestamps
    tab_hlegs: pandas.DataFrame
    """
    dropsnds_cld = dropsnds.where( # P-3 sondes in cld time interval.
        (dropsnds['launch_time']>t1_cld) 
        & (dropsnds['launch_time']<t2_cld)
        & (dropsnds['platform']=='P3'),
        drop=True
        )
    return dropsnds_cld



def table_cldmodules(path_cldmodtable, datetime_type=True):
    """
    datetime_type: bool.
        Default = True, in which case all datetime columns are converted to 
        datetime-type.
    """
    cldmods = pd.read_csv(path_cldmodtable)
    
    if datetime_type:
        dtimekeys = ['start_datetime', 'end_datetime']
        cldmods[dtimekeys] = cldmods[dtimekeys].apply(pd.to_datetime)
    
    return cldmods



if __name__=="__main__":
    path_drpsnds = ("../../data/sondes/EUREC4A_JOANNE_Dropsonde-RD41_"
                    "Level_3_v0.7.0+2.g4a878b3.dirty.nc"
                    )
    path_cldmodtable = "./p3_cloudmodules.csv"
    dirpath_save = "./cldmod_datafiles/"
    if not os.path.isdir(dirpath_save):
        os.mkdir(dirpath_save)
    
    create_dropsondefiles(path_drpsnds, path_cldmodtable, dirpath_save)