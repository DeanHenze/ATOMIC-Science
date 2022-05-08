# -*- coding: utf-8 -*-
"""
Created on Sun May  8 14:18:40 2022

@author: Dean

QC on the level leg data files with regard to air craft roll. High roll likely 
means the high frequency wind and temperature data are poor. 

Solutions to try include:
    imputation
    completely remove bad data segments on either end of a leg
    
Current status:
    Includes only imputation with the mean.
"""


# Built in
import os

# Third party:
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt



def meanimpute(data_df, varkeys_impute, roll_crit=5):
    """
    Imput high aircraft roll time segments with the mean. Appends a columns 
    flaggin rows which have been imputed.
    
    Inputs
    ------
    data_df: pandas.DataFrame.
        Time series data. Includes roll data as a column with key 'roll'.
        
    varkeys_imput: list of str's.
        Keys in data_df to impute.
        
    roll_crit: float.
        Critical roll angle above which data will be imputed. Should be >= 0.
    """
    highroll = abs(data_df['roll']) > roll_crit
    mean = data_df.loc[~highroll].mean()
        
    # Impute:
    vals_impute = mean[varkeys_impute]
    for k in varkeys_impute:
        data_df.loc[highroll, k] = vals_impute[k] 
        
        # Flag for data points which have been imputed:
    data_df['imputed'] = highroll
    
    return data_df



def gaussianimpute(data_df, varkeys_impute):
    """
    In progress...
    """
    highroll = abs(data_df['roll']) > 5
    mean = data_df.loc[~highroll].mean()
    std = data_df.loc[~highroll].std()
    
    
    # Impute:
    #vals_impute = mean[varkeys_impute]
    #for k in varkeys_impute:
    #    data_df.loc[highroll, k] = vals_impute[k] 
        
        # Flag for data points which have been imputed:
    #data_df['imputed'] = highroll
    
    #return data_df
    


"""
fnames_levlegs = os.listdir("./levlegdata_5hz/")
fnames_levlegs = [f for f in fnames_levlegs if f.endswith(".nc")]

# Keys of variables to impute:
varkeys = ["u","v","w","T","q","qD"]
varkeys = varkeys + [k+"'" for k in varkeys]

for fname in fnames_levlegs:
    data = xr.load_dataset(fname)
    data_df = data.to_dataframe() # Easier / faster to work with pandas.
    data_df = meanimpute(data_df, varkeys, roll_crit=5) # Impute
    #data_df.to_xarray().to_netcdf(fname) # Back to xarray object to save.
"""














