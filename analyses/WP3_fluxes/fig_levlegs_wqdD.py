# -*- coding: utf-8 -*-
"""
Created on Sat May 14 16:49:36 2022

@author: Dean

Save PDFs of q' vs w' for cloud groupings and altitude groupings.

Current status:
--------------
- For some of my profiles, the downdrafts are enriched in comparison to the 
  updrafts. Two possibilities:
      1. Precipitation removing heavy isotopes from updrafts then evaporating into 
         downdrafts.
      2. Condensate, removing heavy isotopes from the vapor which the aircraft 
         may not pick up. Then in the downdrafts, the vapor re-evaporates. 
"""



# Built in
import os
import itertools

# Third party
import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import gaussian_kde
import seaborn as sns
import matplotlib.pyplot as plt



## I/O paths / fnames
##_____________________________________________________________________________
# Level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]

# Level leg altitude table:
path_levlegalt = "./p3cld_levleg_altitudes.csv"

# Save directory:
path_savedir = "./levelleg_PDFs"
##_____________________________________________________________________________
## I/O paths / fnames



def create_qwPDFfiles(dir_5hzdata, fnames_levlegs, path_levlegalt, path_savedir):
    """
    Create and save joint w', q' PDFs for each cloud and altitude group.
    """
    
    # Load level leg altitude table:
    levlegalts = pd.read_csv(path_levlegalt)
    
    
    # Append column for cloud group flag to table:
    ncld_g1 = [1, 5, 4, 8]
    ncld_g2 = [7, 9, 11, 10, 12, 6]
    ncld_g3 = [15, 3, 2, 13, 16, 14]
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    
    flag_cldgroup = [1, 2, 3]
    in_cg1 = [x in ncld_g1 for x in levlegalts['ncld']]
    in_cg2 = [x in ncld_g2 for x in levlegalts['ncld']]
    in_cg3 = [x in ncld_g3 for x in levlegalts['ncld']]
    levlegalts['cldgroup'] = np.zeros(len(levlegalts.index))
    levlegalts.loc[in_cg1, 'cldgroup'] = flag_cldgroup[0]
    levlegalts.loc[in_cg2, 'cldgroup'] = flag_cldgroup[1]
    levlegalts.loc[in_cg3, 'cldgroup'] = flag_cldgroup[2]
    
    
    # Append column for altitude group flag to table:
    alt_thresh1 = 300
    alt_thresh2 = 800
    alt_thresh3 = 1600
    alt_thresh4 = 2500
    
    flag_altgroup = [1, 2, 3, 4, 5]
    altgroup = []
    for x in levlegalts['altmean']:
        if x < alt_thresh1:
            altgroup.append(flag_altgroup[0]); continue
        elif x < alt_thresh2:
            altgroup.append(flag_altgroup[1]); continue
        elif x < alt_thresh3:
            altgroup.append(flag_altgroup[2]); continue
        elif x < alt_thresh4:
            altgroup.append(flag_altgroup[3]); continue
        else:
            altgroup.append(flag_altgroup[4])
    levlegalts['altgroup'] = altgroup
    
    
    varkeys = ["w'","q'","roll","dD'"]

    
    # Create and save files:
    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)
    for pair in itertools.product(flag_altgroup, flag_cldgroup):
        
        print("Working on PDF for altitude group %i, cloud group %i" % pair)
        
        tableinfo_pair = levlegalts.loc[
            (levlegalts['altgroup']==pair[0])
            & (levlegalts['cldgroup']==pair[1])
            ]
        pdf = test_scatter1(
            dir_5hzdata, fnames_levlegs, 
            tableinfo_pair['ncld'], tableinfo_pair['nlevleg'], 
            varkeys
            )
        
        #fname_save = "jointPDF_wprimeqprime_altgrp%i_cldgrp%i.csv" % pair
        #pdf.to_csv(os.path.join(path_savedir, fname_save))



def test_scatter1(dir_5hzdata, fnames_levlegs, 
                    ncld_list, nlevleg_list, varkeys):
    """
    Generate PDF (using KDE method) for all P-3 level leg data for a subset of 
    cloud module and level leg numbers.
    
    Returns
    -------
    Joint-PDF as a pd.DataFrame with q' as index and w' as columns.
    If no data available for the set of level legs passed, returns an empty 
    DataFrame.
    """
    # All level leg data for cloud + altitude group:
    fnames_grp = fnamesubset(fnames_levlegs, ncld_list, nlevleg_list)    
    data_grp = multincdata_todf(dir_5hzdata, fnames_grp, varkeys)
    
    
    if len(data_grp.index) == 0: return pd.DataFrame({}) # No data for this group.
        
    
    # Remove data where the roll was greater than 5 degrees:
    roll_crit = 5
    highroll = abs(data_grp['roll']) > roll_crit
    data_grpqc = data_grp.loc[~highroll]
    
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(data_grpqc["w'"], data_grpqc["q'"], c=data_grpqc["dD'"])
    plt.colorbar()
        
    

def test_scatter2(dir_5hzdata, fnames_levlegs, nclds, 
                  qrange1=(0, 0.025), qrange2=(0.975, 1.)):
    """
    Generate PDF (using KDE method) for all P-3 level leg data for a subset of 
    cloud module and level leg numbers.
    
    Returns
    -------
    Joint-PDF as a pd.DataFrame with q' as index and w' as columns.
    If no data available for the set of level legs passed, returns an empty 
    DataFrame.
    """
    nlevleg_list = [1,2,3,4,5,6]
    #nclds = np.arange(1,17,1)


    q_up = []
    q_down = []
    q = []
    dD_up = []
    dD_down = []
    dD = []
    alt = []
    ncld_list = []
    
    #results = []
    
    
    for n in nclds:
    # All level leg data for cloud + altitude group:
        try:
            fnames_grp = fnamesubset(
                fnames_levlegs, (np.ones(6)*n).astype(int), nlevleg_list)    
        except IndexError:
            try:
                fnames_grp = fnamesubset(
                    fnames_levlegs, (np.ones(5)*n).astype(int), nlevleg_list[:-1])    
            except IndexError:
                fnames_grp = fnamesubset(
                    fnames_levlegs, (np.ones(4)*n).astype(int), nlevleg_list[:-2])    
            

        
        for f in fnames_grp:
            
            data = xr.load_dataset(dir_5hzdata + f)
            data = data.to_dataframe()
        
            
            # Remove data where the roll was greater than 5 degrees:
            roll_crit = 5
            highroll = abs(data['roll']) > roll_crit
            data_qc = data.loc[~highroll]
        
            if qrange1[1]=="w'=0":
                downdraft = (
                    (data_qc["w'"] > np.quantile(data_qc["w'"], qrange1[0])) & 
                    (data_qc["w'"] < 0)
                    )
                updraft = (
                    (data_qc["w'"] > 0) & 
                    (data_qc["w'"] < np.quantile(data_qc["w'"], qrange2[1]))
                    )                
        
            else:
                downdraft = (
                    (data_qc["w'"] > np.quantile(data_qc["w'"], qrange1[0])) & 
                    (data_qc["w'"] < np.quantile(data_qc["w'"], qrange1[1]))
                    )
                updraft = (
                    (data_qc["w'"] > np.quantile(data_qc["w'"], qrange2[0])) & 
                    (data_qc["w'"] < np.quantile(data_qc["w'"], qrange2[1]))
                    )

            data_down = data_qc.loc[downdraft]
            data_up = data_qc.loc[updraft]
            env = (
                (data_qc["w'"] > np.quantile(data_qc["w'"], 0.05)) & 
                (data_qc["w'"] < np.quantile(data_qc["w'"], 0.95))
                )
            data_env = data_qc[env]
            
            q_up.append(data_up["q"].mean())
            q_down.append(data_down["q"].mean())
            q.append(data_env["q"].mean())
    
            dD_up.append(data_up["dD"].mean())
            dD_down.append(data_down["dD"].mean())
            dD.append(data_env["dD"].mean())
            
            alt.append(data_qc["alt"].mean())
            ncld_list.append(np.ones(len(fnames_grp))*n)
            
    
    
    results = {
        'qup':np.array(q_up), 'qdown':np.array(q_down), 'q':np.array(q), 
        'dDup':np.array(dD_up), 'dDdown':np.array(dD_down), 'dD':np.array(dD),  
        'alt':np.array(alt)                       
        } 
    return results

        


def multincdata_todf(pathdir, fnames, varkeys):
    """
    Get all .nc file data for fnames located in directory pathdir. Return as 
    a single pandas df where columns correspond to varkeys (subset of keys in 
    the .nc files). Also append a column indicating which rows belong to the 
    same data file.
    
    Assumes that the .nc dataset is in a format convertable to a 
    pandas DataFrame with time as the dimension (e.g. will be converted to 
    the index).
    """
    data_alldf = pd.DataFrame({})
    ndatafile=1
    for f in fnames:
        data_nc = xr.load_dataset(pathdir + f)
        data_df = data_nc[varkeys].to_dataframe()
        data_df['ndatafile'] = ndatafile*np.ones(len(data_df.index))
        ndatafile += 1
        data_alldf = data_alldf.append(data_df)
    return data_alldf



def fnamesubset(fnames_levlegs, ncld_list, nlevleg_list):
    """
    Return subset of level leg filenames for each pair of (cloud module, level 
    leg number) in lists (ncld_list, nlevleg_list).
    """
    fnames = []
    for ncld, nlevleg in zip(ncld_list, nlevleg_list):
        
        ncld_str = str(ncld).zfill(2)
        fname = [f for f in fnames_levlegs 
                 if "_cld%s" % ncld_str in f 
                 and "_levleg%i" % nlevleg in f]
        fnames.append(fname[0]) # Should only find one filename per pair.    
    return fnames

    
   
if __name__=="__main__":
    #create_qwPDFfiles(dir_5hzdata, fnames_levlegs, path_levlegalt, path_savedir)
    
    
    # Cloud module number groupings:
    ncld_g1 = [1, 5, 4, 8]
    ncld_g2 = [7, 9, 11, 10, 12, 6]
    ncld_g3 = [15, 3, 2, 13, 16, 14]
    
    
    results_g1 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g1)
    results_g2 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g2)
    results_g3 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g3)
    
    
    #control_g1 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g1, 
    #                           qrange1=(0.05, 0.5), qrange2=(0.5, 0.95))
    #control_g2 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g2, 
    #                           qrange1=(0.05, 0.5), qrange2=(0.5, 0.95))
    #control_g3 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g3, 
    #                           qrange1=(0.05, 0.5), qrange2=(0.5, 0.95))
    control_g1 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g1, 
                               qrange1=(0.05, "w'=0"), qrange2=(0.5, 0.95))
    control_g2 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g2, 
                               qrange1=(0.05, "w'=0"), qrange2=(0.5, 0.95))
    control_g3 = test_scatter2(dir_5hzdata, fnames_levlegs, ncld_g3, 
                               qrange1=(0.05, "w'=0"), qrange2=(0.5, 0.95))
    
    
    ## Plot deviations of updrafts and downdrafts from mean
    ##_________________________________________________________________________
    fig = plt.figure(figsize=(6.5, 3.5))
    ax1 = fig.add_axes([0.11, 0.15, 0.275, 0.75])
    ax2 = fig.add_axes([0.41, 0.15, 0.275, 0.75])
    ax3 = fig.add_axes([0.71, 0.15, 0.275, 0.75])
    axset = [ax1, ax2, ax3]
    scattersize = 15
    axset[0].scatter(
        results_g1['dDup']-results_g1['dD'], results_g1['alt'], 
        color='blue', s=scattersize, label='updraft'
        )
    axset[0].scatter(
        results_g1['dDdown']-results_g1['dD'], results_g1['alt'], 
        color='magenta', s=scattersize, label='downdraft'
        ) 
    axset[0].scatter(
        control_g1['dDup']-control_g1['dD'], control_g1['alt'], 
        color='grey', s=scattersize, marker='x', label='upward \nturb.', 
        )
    axset[0].scatter(
        control_g1['dDdown']-control_g1['dD'], control_g1['alt'], 
        color='grey', s=scattersize, marker='.', label='downward \nturb.', 
        )
    
    
    axset[1].scatter(
        results_g2['dDup']-results_g2['dD'], results_g2['alt'], 
        color='blue', s=scattersize
        )
    axset[1].scatter(
        results_g2['dDdown']-results_g2['dD'], results_g2['alt'], 
        color='magenta', s=scattersize
        )    
    axset[1].scatter(
        control_g2['dDup']-control_g2['dD'], control_g2['alt'], 
        color='grey', s=scattersize, marker='x', label='upward turbulence', 
        )
    axset[1].scatter(
        control_g2['dDdown']-control_g2['dD'], control_g2['alt'], 
        color='grey', s=scattersize, marker='.', label='downward turbulence', 
        )
    
    axset[2].scatter(
        results_g3['dDup']-results_g3['dD'], results_g3['alt'], 
        color='blue', s=scattersize
        )
    axset[2].scatter(
        results_g3['dDdown']-results_g3['dD'], results_g3['alt'], 
        color='magenta', s=scattersize
        )
    axset[2].scatter(
        control_g3['dDup']-control_g3['dD'], control_g3['alt'], 
        color='grey', s=scattersize, marker='x', label='upward turb.', 
        )
    axset[2].scatter(
        control_g3['dDdown']-control_g3['dD'], control_g3['alt'], 
        color='grey', s=scattersize, marker='.', label='downward turb.', 
        )
    
    
    axset[1].set_yticks(axset[1].get_yticks())
    axset[1].set_yticklabels(['' for t in axset[1].get_yticks()])
    axset[2].set_yticks(axset[2].get_yticks())
    axset[2].set_yticklabels(['' for t in axset[2].get_yticks()])


    for ax in axset: 
        ax.set_xlim(-15, 22)
        ax.set_xticks(np.arange(-10, 21, 5))
        ax.set_xticklabels(['', '-5', '', '5', '', '15', ''])
        ax.vlines(0, -100, 3500, colors='black')
        
        ax.set_ylim(-50, 3300)
    

    axset[1].set_xlabel(r'$\Delta \delta D$'+u'(\u2030)', fontsize=12)
    axset[0].set_ylabel('altitude (m)', fontsize=12)
    
    axlabels = ['cg1', 'cg2', 'cg3']
    for ax, lab in zip(axset, axlabels):
        ax.text(
            0.5, 1.01, lab, fontsize=12, 
            ha='center', va='bottom', transform=ax.transAxes
            )
        
    axset[0].legend(loc='lower right', handletextpad=0.1, fontsize=8)
    fig.savefig("./dD_updraft-downdraft_deviations.png")
    ##_________________________________________________________________________
    ## Plot deviations of updrafts and downdrafts from mean
    
    
    
    ## Plot deviations of updrafts and downdrafts from mean
    ##_________________________________________________________________________
    fig = plt.figure(figsize=(6.5, 3.5))
    ax1 = fig.add_axes([0.11, 0.15, 0.275, 0.75])
    ax2 = fig.add_axes([0.41, 0.15, 0.275, 0.75])
    ax3 = fig.add_axes([0.71, 0.15, 0.275, 0.75])
    axset = [ax1, ax2, ax3]
    scattersize = 15
    axset[0].scatter(
        results_g1['dDup']-results_g1['dDdown'], results_g1['alt'], 
        color='black', s=scattersize
        )  
    axset[1].scatter(
        results_g2['dDup']-results_g2['dDdown'], results_g2['alt'], 
        color='black', s=scattersize
        )
    axset[2].scatter(
        results_g3['dDup']-results_g3['dDdown'], results_g3['alt'], 
        color='black', s=scattersize
        )
    
    
    axset[1].set_yticks(axset[1].get_yticks())
    axset[1].set_yticklabels(['' for t in axset[1].get_yticks()])
    axset[2].set_yticks(axset[2].get_yticks())
    axset[2].set_yticklabels(['' for t in axset[2].get_yticks()])


    for ax in axset: 
        ax.set_xlim(-15, 20)
        ax.set_xticks(np.arange(-10, 21, 5))
        ax.set_xticklabels(['', '-5', '', '5', '', '15', ''])
        ax.vlines(0, -100, 3500, colors='black')
        
        ax.set_ylim(-50, 3300)
    

    axset[1].set_xlabel(r'$\Delta \delta D(up-down) $'+u'(\u2030)', fontsize=12)
    axset[0].set_ylabel('altitude (m)', fontsize=12)
    
    axlabels = ['cg1', 'cg2', 'cg3']
    for ax, lab in zip(axset, axlabels):
        ax.text(
            0.5, 1.01, lab, fontsize=12, 
            ha='center', va='bottom', transform=ax.transAxes
            )
        
    fig.savefig("./dD_updraft-downdraft_diff.png")
    ##_________________________________________________________________________
    ## Plot deviations of updrafts and downdrafts from mean


