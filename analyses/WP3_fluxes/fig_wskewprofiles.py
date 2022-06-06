# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:15:49 2022

@author: Dean

Mean profiles of vertical velocity skewness for P-3 cloud modules. Separate 
vertical profiles for each cloud grouping.
"""



# Third party:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Local code:
import profileplotter
import rangescaler



# Cloud groups:
ncld_g1 = [7, 9, 11, 10, 12, 6]
ncld_g2 = [15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4, 8]

ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()



## Load and merge data tables:
wskew = pd.read_csv("./wskew_levlegs/WP3_wskewness_levlegs.csv")
keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
data = wskew.merge(
    keyalts_table, 
    left_on="ncld", right_on="ncld", how="left"
    )


# Restructure data into a dataframe where each column is a separate profile 
# and the index is scaled altitude:
def restruct_wksewprofiles(data, ncld_list, keyalts_table, scale_altkeys=[]):
    """
    Inputs
    ------
    data: pandas.DataFrame.
    
    ncld_list: iterable of ints.
        Subset of cloud module numbers to include in returned results. 
    
    keyalts_table: pandas.DataFrame.
        If scaling altitude by a subset of columns in this table, those 
        columns cannot have NAN's.
        
    scale_altkeys: iterable of str's.
        Optional. Keys in keyalts_table to scale altitude by. If not empty,
        altitudes will be mapped as:
            0 <-- 0
            1 <-- 1st key in scale_alt
            2 <-- 2nd key in scale_alt
            3 ...
    """
    wskew_prfs = pd.DataFrame({})
    
    for ncld, data_cld in data.groupby(by='ncld'):
    
        if ncld in ncld_list:
            
            # Scale altitude:
            keyalts = keyalts_table.loc[keyalts_table['ncld']==ncld]
            alt_scalepoints = [0] + [keyalts[k].item() for k in scale_altkeys]
            data_cld['altleg_scaled'] = rangescaler.piecewise_linscale(
                data_cld['altleg'].values, 
                alt_scalepoints, np.arange(len(alt_scalepoints))
                )
            
            # Append:
            data_cld.set_index('altleg_scaled', inplace=True)
            data_cld = data_cld.rename(columns={"wskew": "wskew_cld%i"%ncld})
            wskew_prfs = pd.merge(
                    wskew_prfs, data_cld[["wskew_cld%i"%ncld]], 
                    #how='outer', left_on='altleg_scaled', right_on='altleg_scaled'
                    how='outer', left_index=True, right_index=True
                    )   
            
    return wskew_prfs
        
    

if __name__=="__main__":
        
    ## Figure where altitude is scaled by LCL and cloud top height.
    ##_________________________________________________________________________
    # Get wskew profiles for each cloud group as a list of Dataframes:
    wskewprfs_lclctscaling = []
    for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
        wskewprfs_lclctscaling.append(
            restruct_wksewprofiles(
                data, ncld_list, 
                keyalts_table, scale_altkeys=['z_lcl', 'z_ctmean_50p95p']
                )
            )
        

    pltcolors = ['grey', 'blue', 'red'] # Plot colors for each cloud group.
    
    
    # Plot:
    fig = plt.figure(figsize=(4,8))
    ax = fig.add_axes([0.3, 0.1, 0.65, 0.8])
    for prfs, c in zip(wskewprfs_lclctscaling, pltcolors):  
            profileplotter.plotprf_singlevar(prfs, ax, pcolor=c)    
       
        # Vertical line at wskew=0 for visual reference:
    ax.plot([0,0], [0,4], 'k-', linewidth=1)
    
    
    # Figure limits, labels, etc:
    ax.set_xlim(-1.5, 2.5)
    ax.set_ylim(0, 3)
    yticks = [0,1,2,3]
    yticklabels = ["0", r"$z_{LCL}$", r"$z_{CT}$", 
                   r"$z_{LCL}$+2$\Delta z_{(CT-LCL)}$"
                   ] 
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel(r"<w'w'w'>", fontsize=12)
    
    
    fig.savefig("./fig_wskewprofiles_LCLCTscaling.png")
    ##_________________________________________________________________________    
    ## Figure where altitude is scaled by LCL and cloud top height.
    


    ## Figure where altitude is scaled by LCL and trade inversion bottom.
    ##_________________________________________________________________________
    # Get wskew profiles for each cloud group as a list of Dataframes:
    wskewprfs_lcltibscaling = []
    for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
        wskewprfs_lcltibscaling.append(
            restruct_wksewprofiles(
                data, ncld_list, 
                keyalts_table, scale_altkeys=['z_lcl', 'z_tib']
                )
            )
        

    pltcolors = ['grey', 'blue', 'red'] # Plot colors for each cloud group.
    
    
    # Plot:
    fig = plt.figure(figsize=(4,8))
    ax = fig.add_axes([0.3, 0.1, 0.65, 0.8])
    for prfs, c in zip(wskewprfs_lcltibscaling, pltcolors):  
            profileplotter.plotprf_singlevar(prfs, ax, pcolor=c)    
       
        # Vertical line at wskew=0 for visual reference:
    ax.plot([0,0], [0,4], 'k-', linewidth=1)
    
    
    # Figure limits, labels, etc:
    ax.set_xlim(-1.5, 2.5)
    ax.set_ylim(0, 3.5)

    yticks = [0,1,2,3]
    yticklabels = ["0", r"$z_{LCL}$", r"$z_{IB}$", 
                   r"$z_{LCL} + 2\Delta z_{(IB-LCL)}$", 
                   ]  
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel(r"<w'w'w'>", fontsize=12)
    
    
    fig.savefig("./fig_wskewprofiles_LCLTIBscaling.png")
    ##_________________________________________________________________________    
    ## Figure where altitude is scaled by LCL and trade inversion bottom.