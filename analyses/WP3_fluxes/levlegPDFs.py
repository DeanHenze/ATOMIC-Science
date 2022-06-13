# -*- coding: utf-8 -*-
"""
Created on Sat May 14 16:49:36 2022

@author: Dean

Save PDFs of q' vs w' for cloud groupings and altitude groupings.
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
    
    
    varkeys = ["w'","q'","roll"]

    
    # Create and save files:
    if not os.path.isdir(path_savedir): os.mkdir(path_savedir)
    for pair in itertools.product(flag_altgroup, flag_cldgroup):
        
        print("Working on PDF for altitude group %i, cloud group %i" % pair)
        
        tableinfo_pair = levlegalts.loc[
            (levlegalts['altgroup']==pair[0])
            & (levlegalts['cldgroup']==pair[1])
            ]
        pdf = kde_levleggroup(
            dir_5hzdata, fnames_levlegs, 
            tableinfo_pair['ncld'], tableinfo_pair['nlevleg'], 
            varkeys
            )
        
        fname_save = "jointPDF_wprimeqprime_altgrp%i_cldgrp%i.csv" % pair
        pdf.to_csv(os.path.join(path_savedir, fname_save))



def kde_levleggroup(dir_5hzdata, fnames_levlegs, 
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
    
    
    # Column for total time of each level leg:
    tottime_df = pd.DataFrame({})
    for n, data in data_grpqc.groupby(by="ndatafile"):
        tottime_df = tottime_df.append(
            {'ndatafile':n, 'time_total': len(data.index)}, 
            ignore_index=True
            )
    data_grpqc = pd.merge(data_grpqc, tottime_df, how='left', on='ndatafile')
    
    
    # 2D grid points to get PDF at:
    dw = 0.1
    dq = 0.1
    wmin = data_grpqc["w'"].min()
    wmax = data_grpqc["w'"].max()
    qmin = data_grpqc["q'"].min()
    qmax = data_grpqc["q'"].max()
        # Extra cushion for PDF domain:
    wspan = wmax-wmin
    wmin = wmin - 0.25*wspan
    wmax = wmax + 0.25*wspan
    qspan = qmax-qmin
    qmin = qmin - 0.25*qspan
    qmax = qmax + 0.25*qspan
    
    w_1dgrid = np.arange(wmin, wmax, dw)
    q_1dgrid = np.arange(qmin, qmax, dq)
    w_2dgrid, q_2dgrid = np.meshgrid(w_1dgrid, q_1dgrid)

    # KDE bandwidth:
    neff = len(data_grpqc.index) # effective number of data points
    d = 2 # number of dims.
    bw_silv = (neff * (d + 2) / 4.)**(-1. / (d + 4)) # Silverman's method.    
    fact = 2 # multiplicative factor for bandwith (higher -> smoother PDF)
    bw = bw_silv*fact
    
    # KDE estimation and evaluation at gridpoints:
    kernel = gaussian_kde(
        data_grpqc[["w'","q'"]].values.T, 
        bw_method=bw, weights=1/data_grpqc["time_total"]
        )
    pdf = kernel(np.vstack([w_2dgrid.ravel(), q_2dgrid.ravel()]))
    pdf = pdf.reshape(w_2dgrid.shape)
    
    # Collect into a pandas df and return:
    return pd.DataFrame(
        pdf, index=q_1dgrid, columns=w_1dgrid
        )
    


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
    create_qwPDFfiles(dir_5hzdata, fnames_levlegs, path_levlegalt, path_savedir)

        
    
    

        


"""

for ncldgrp, pltc in zip(cldgrps, pltchars):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==1)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    test = kde_levleggroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, pltc
        )
"""    
    

"""
plt.figure()

tableinfo_grp = levlegalt_tab.loc[
    (levlegalt_tab['altgroup']==1)
    & (levlegalt_tab['cldgroup']==1)
    ]
test = kde_levleggroup(
    dir_5hzdata, fnames_levlegs, 
    tableinfo_grp['ncld'], tableinfo_grp['nlevleg'])
plt.contour(test.columns, test.index, test, colors='grey')
    

tableinfo_grp = levlegalt_tab.loc[
    (levlegalt_tab['altgroup']==1)
    & (levlegalt_tab['cldgroup']==2)
    ]
test = kde_levleggroup(
    dir_5hzdata, fnames_levlegs, 
    tableinfo_grp['ncld'], tableinfo_grp['nlevleg'])
plt.contour(test.columns, test.index, test, colors='blue')


tableinfo_grp = levlegalt_tab.loc[
    (levlegalt_tab['altgroup']==1)
    & (levlegalt_tab['cldgroup']==3)
    ]
test = kde_levleggroup(
    dir_5hzdata, fnames_levlegs, 
    tableinfo_grp['ncld'], tableinfo_grp['nlevleg'])
plt.contour(test.columns, test.index, test, colors='red')
"""    

    






