# -*- coding: utf-8 -*-
"""
Created on Sat May 14 16:49:36 2022

@author: Dean

PDFs of q' vs w' for cloud groupings and altitude groupings.
"""



# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns

# Local code
import iso
import cdf



# Cloud groups:
ncld_g1 = [7, 9, 11, 10, 12, 6]
ncld_g2 = [15, 3, 2, 13, 16, 14]
ncld_g3 = [1, 5, 4, 8]

ncld_g1.sort()
ncld_g2.sort()
ncld_g3.sort()



# Level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]


# Level leg altitude table:
levlegalt_tab = pd.read_csv("./p3cld_levleg_altitudes.csv")


# Group level legs by cloud group and altitude:
in_cg1 = [x in ncld_g1 for x in levlegalt_tab['ncld']]
in_cg2 = [x in ncld_g2 for x in levlegalt_tab['ncld']]
in_cg3 = [x in ncld_g3 for x in levlegalt_tab['ncld']]
levlegalt_tab['cldgroup'] = np.zeros(len(levlegalt_tab.index))
levlegalt_tab.loc[in_cg1, 'cldgroup'] = 1
levlegalt_tab.loc[in_cg2, 'cldgroup'] = 2
levlegalt_tab.loc[in_cg3, 'cldgroup'] = 3


alt_thresh1 = 300
alt_thresh2 = 800
alt_thresh3 = 1600
alt_thresh4 = 2500

altgroup = []
for x in levlegalt_tab['altmean']:
    if x < alt_thresh1:
        altgroup.append(1); continue
    elif x < alt_thresh2:
        altgroup.append(2); continue
    elif x < alt_thresh3:
        altgroup.append(3); continue
    elif x < alt_thresh4:
        altgroup.append(4); continue
    else:
        altgroup.append(5)
levlegalt_tab['altgroup'] = altgroup



def multincdata_todf(pathdir, fnames, varkeys):
    """
    Get all .nc file data for fnames located in directory pathdir. Return as 
    a single pandas df where columns correspond to varkeys (subset of keys in 
    the .nc files).
    
    Assumes that the .nc dataset is in a format convertable to a 
    pandas DataFrame.
    """
    data_alldf = pd.DataFrame({})
    for f in fnames:
        data_nc = xr.load_dataset(pathdir + f)
        data_alldf = data_alldf.append(data_nc[varkeys].to_dataframe())
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



varkeys = ["w'","q'","roll"]
 
   
def plotkde_cldaltgroup(dir_5hzdata, fnames_levlegs, 
                        ncld_list, nlevleg_list, 
                        ax, plt_kwargs={}):
    """
    Plot PDF (using KDE method) for all P-3 level leg data in a specific 
    cloud group + altidue group combination.
    """
    # All level leg data for cloud + altitude group:
    fnames_grp = fnamesubset(fnames_levlegs, ncld_list, nlevleg_list)    
    data_grp = multincdata_todf(dir_5hzdata, fnames_grp, varkeys)
    
    if len(data_grp.index) !=0:
        
        # Remove data where the roll was greater than 5 degrees:
        roll_crit = 5
        highroll = abs(data_grp['roll']) > roll_crit
        data_grpqc = data_grp.loc[~highroll]
        
        # KDE plot:
        dw = 0.04
        dq = 0.1
        wmin = data_grpqc["w'"].min() - 4*dw
        wmax = data_grpqc["w'"].max() + 4*dw
        qmin = data_grpqc["q'"].min() - 4*dq
        qmax = data_grpqc["q'"].max() + 4*dq
        #wmin, qmin = data_grpqc[["w'","q'"]].min()
        #wmax, qmax = data_grpqc[["w'","q'"]].max()
        ww, qq = np.meshgrid(
            np.arange(wmin, wmax, dw), np.arange(qmin, qmax, dq))
        neff = len(data_grpqc.index) # effective number of data points
        d = 2 # number of dims.
        fact = 1.4 # multiplicative factor for bandwith (higher -> more smooth)
        bw_silv = (neff * (d + 2) / 4.)**(-1. / (d + 4))
        kernel = gaussian_kde(data_grpqc[["w'","q'"]].values.T, bw_method=bw_silv*fact)
        pdf = kernel(np.vstack([ww.ravel(), qq.ravel()]))
        pdf = pdf.reshape(ww.shape)
        ccdf_levs = [1-0.001, 1-0.01] + list(1 - np.arange(0.05, 1, 0.2))
        pvals_cdf = cdf.pvals_ccdflevs(pdf, ccdf_levs, dx=dw, dy=dq)

        
        #plt.contour(ww, qq, pdf, levels=pvals_cdf, **plt_kwargs)
        sns.kdeplot(
            data=data_grpqc, x="w'", y="q'", ax=ax, 
            levels=[0.001, 0.01] + list(np.arange(0.05, 1, 0.2)), 
            bw_adjust=1.25, 
            **plt_kwargs
            )
        
        
    
pltcolors = ["grey", "blue", "red"]
cldgrps = [1,2,3]
#pltchars = [
#    {'cmap':'binary', 'fill':True, 'extend':'max', 'linewidths':1}, 
#    {'color':'blue', 'fill':False, 'linewidths':1}, 
#    {'color':'red', 'fill':False, 'linewidths':1}
#    ]
pltchars = [
    {'colors':'grey', 'fill':True, 'extend':'max', 'linewidths':1}, 
    {'colors':'blue', 'fill':False, 'linewidths':1}, 
    {'colors':'red', 'fill':False, 'linewidths':1}
    ]


plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltchars):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==1)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, pltc
        )
    
    
plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltchars):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==2)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, pltc
        )
    
    
plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltcolors):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==2)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, {'color':pltc}
        )
    
    
plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltcolors):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==3)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, {'color':pltc}
        )
    
    
plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltcolors):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==4)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, {'color':pltc}
        )
    

plt.figure()
ax = plt.axes()
for ncldgrp, pltc in zip(cldgrps, pltcolors):
    tableinfo_grp = levlegalt_tab.loc[
        (levlegalt_tab['altgroup']==5)
        & (levlegalt_tab['cldgroup']==ncldgrp)
        ]
    plotkde_cldaltgroup(
        dir_5hzdata, fnames_levlegs, 
        tableinfo_grp['ncld'], tableinfo_grp['nlevleg'], 
        ax, {'color':pltc}
        )
    
    
"""
tableinfo_grp1 = levlegalt_tab.loc[
    (levlegalt_tab['altgroup']==1)
    & (levlegalt_tab['cldgroup']==1)
    ]
fnames_g1 = fnamesubset(fnames_levlegs, tableinfo_grp1['ncld'], tableinfo_grp1['nlevleg'])    
data_group1 = multincdata_todf(dir_5hzdata, fnames_g1, varkeys)
# Remove data where the roll was greater than 5 degrees:
roll_crit = 5
highroll = abs(data_group1['roll']) > roll_crit
data_group1qc = data_group1.loc[~highroll]

plt.figure()
ax = plt.axes()
sns.kdeplot(
    data=data_group1qc, x="w'", y="q'", ax=ax, 
    levels=[0.01] + list(np.arange(0.05, 1, 0.1)), color='grey'
    )
"""


def plots(ncld):
    
    ncld_str = str(ncld).zfill(2)
    fnames_levlegs_cld = [f for f in fnames_levlegs if "_cld%s"%ncld_str in f]
   
    for f in fnames_levlegs_cld:
        
        data = xr.load_dataset(dir_5hzdata + f)
        data_df = data.to_dataframe() # Easier / faster to work with pandas.
        
        # Add dD column:
        data_df["dD"] = iso.qD_dD_convert('qD2dD', data_df["q"], data_df["qD"])
        
        # Remove data where the roll was greater than 5 degrees:
        roll_crit = 5
        highroll = abs(data_df['roll']) > roll_crit
        data_dfqc = data_df.loc[~highroll]
        
        # Split into rows of upward vs downward velocities:
        dataup = data_dfqc.loc[data_dfqc["w'"]>0]
        datadown = data_dfqc.loc[data_dfqc["w'"]<0]
        
        # Plots:
        plt.figure()
        ax = plt.axes()
        plt.scatter(data_dfqc["w"], data_dfqc["q"])
        sns.kdeplot(
            data=data_dfqc, x="w", y="q", ax=ax, 
            levels=[0.01] + list(np.arange(0.05, 1, 0.1)), color='red'
            )
        
        #plt.figure()
        #plt.scatter(data_dfqc["w'"], data_dfqc["q'"], s=1)
        
        #plt.figure()
        #plt.hist(dataup['q'], bins=20, histtype='step')
        #plt.hist(datadown['q'], bins=20, histtype='step')
        
        #plt.figure()
        #plt.hist(dataup["w'"], bins=20, histtype='step')
        #plt.hist(datadown["w'"], bins=20, histtype='step')
    
        print("percentage upward = %0.2f" % (100*len(dataup)/len(data_dfqc)))
        print("percentage downward = %0.2f" % (100*len(datadown)/len(data_dfqc)))
        
        
        #plt.figure()
        #ax1 = plt.subplot(1,2,1)
        #qdD_pdf(dataup, 'blue', ax1)
        #qdD_pdf(datadown, 'red', ax1)
        
        
        # Same as above execpt for data above / below the upper / lower quartile:
        #q1, q3 = np.nanquantile(data_dfqc["w'"], [0.25, 0.75])
        #dataup = data_dfqc.loc[data_dfqc["w'"]>q3]
        #datadown = data_dfqc.loc[data_dfqc["w'"]<q1]
        
        #plt.figure()
        #plt.hist(dataup["dD"], bins=20, histtype='step')
        #plt.hist(datadown["dD"], bins=20, histtype='step')
        
        ax.text(
            0.95, 0.95, "z = %i m" % round(data_dfqc["alt"].mean()), 
            ha='right', va='top', transform=ax.transAxes, fontsize=14
            )
        
        
        
def qdD_pdf(data, color, ax, scatter=False):
        
    
        if scatter: ax.scatter(data["q"], data["dD"], s=1, c=color)
        kernel = gaussian_kde([data["q"].values, data["dD"].values])
        qq, dDdD = np.meshgrid(
            np.linspace(np.min(data["q"]), np.max(data["q"])),
            np.linspace(np.min(data["dD"]), np.max(data["dD"]))
            )
        prob = np.reshape(kernel([qq.ravel(), dDdD.ravel()]).T, qq.shape)
        ax.contour(qq, dDdD, prob, colors=color)
        



#plots(6)
#plots(13)


