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
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as linsegcmap

# Local code
import cdf



## I/O paths / fnames
##_____________________________________________________________________________
# Level leg data filenames:
dir_5hzdata = "./levlegdata_5hz/"
fnames_levlegs = [f for f in os.listdir(dir_5hzdata) if f.endswith(".nc")]

# Level leg altitude table:
path_levlegalt = "./p3cld_levleg_altitudes.csv"

# Save directory:
path_wqPDFs = "./levelleg_PDFs"
fnames_wqPDFs = [f for f in os.listdir(path_wqPDFs) if "wprimeqprime" in f]
##_____________________________________________________________________________
## I/O paths / fnames



def cmap_mod(levs, cmap_name, vmin=0, vmax=1.0):
    """
    Map uneven level spacings to even colormap spacings and optional truncation.
    levs: contour levels (array) in increasing order.
    cmap (str or matplotlib.pyplot.cmap): matplotlib colormap to use.
    """
    cmap = plt.get_cmap(cmap_name)
    normedlevs = (levs-levs[0])/(levs[-1]-levs[0]) # Norm to 0-1 scale.
    #colors = cmap(np.linspace(0,1,len(levs))) # Cmap colors at even spacings.
    colors = cmap(np.linspace(vmin,vmax,len(levs))) # Cmap colors at even spacings.
    # Map:
    return linsegcmap.from_list(cmap_name, list(zip(normedlevs,colors)))



altgroups = [1,2,3,4,5]
altgrp_bnds = [300, 800, 1500, 2500]
cldgroups = [1,2,3]



contourargs = {
    1: {'colors':'red', 'linewidths':0.75}, 
    2: {'colors':'dimgrey', 'linewidths':0.75}, 
    3: {'colors':'blue', 'linewidths':0.75}
    }


fig = plt.figure(figsize=(6.5, 8))
ax_width = 0.3
ax_height = 0.25
ax1 = fig.add_axes([0.15, 3*0.075 + 2*ax_height, ax_width, ax_height])
ax2 = fig.add_axes([2*0.15 + ax_width, 3*0.075 + 2*ax_height, ax_width, ax_height])
ax3 = fig.add_axes([0.15, 2*0.075 + ax_height, ax_width, ax_height])
ax4 = fig.add_axes([2*0.15 + ax_width, 2*0.075 + ax_height, ax_width, ax_height])
ax5 = fig.add_axes([(1-ax_width)/2, 0.075, ax_width, ax_height])
axset = [ax1, ax2, ax3, ax4, ax5]


for altgrp, ax in zip(altgroups, axset):
    
    for cldgrp in cldgroups:
        
        f = [f for f in fnames_wqPDFs 
             if "_altgrp%i_cldgrp%i" % tuple([altgrp, cldgrp]) in f]
        f = f[0]
        pdf = pd.read_csv(os.path.join(path_wqPDFs, f), header=0, index_col=0)
        if len(pdf)==0: continue
        cdf_levs = [0.0025, 0.01] + list(np.arange(0.05, 0.8, 0.2))
        w = pdf.columns.values.astype(float)
        q = pdf.index.values.astype(float)
        dw = w[1] - w[0]
        dq = q[1] - q[0]
        pvals_cdf = cdf.pvals_cdflevs(pdf.values, cdf_levs, dx=dw, dy=dq)
        
        
        ax.contour(w, q, pdf, levels=pvals_cdf, **contourargs[cldgrp])
        
        if cldgrp==2:
            cmap_cldgrp2 = cmap_mod(pvals_cdf, 'Greys', vmin=0.1, vmax=0.8)
            ax.contourf(
                w, q, pdf, 
                levels=pvals_cdf, cmap=cmap_cldgrp2, 
                extend='max'
                )
        
        
## Axes limits, labels, ...
for ax in axset: 
    ax.set_xlim(-5, 5)
    ax.set_ylim(-6, 6.5)
    
axset[2].set_ylabel("w' (m/s)", fontsize=14)
axset[4].set_xlabel("q' (g/kg)", fontsize=14)

axset[0].text(0.5, 1., "z < %i m" % altgrp_bnds[0], 
                fontsize=11, ha='center', va='bottom', transform=axset[0].transAxes)
axset[1].text(0.5, 1., "%i < z < %i m" % tuple([altgrp_bnds[0], altgrp_bnds[1]]), 
                fontsize=11, ha='center', va='bottom', transform=axset[1].transAxes)
axset[2].text(0.5, 1., "%i < z < %i m" % tuple([altgrp_bnds[1], altgrp_bnds[2]]), 
                fontsize=11, ha='center', va='bottom', transform=axset[2].transAxes)
axset[3].text(0.5, 1., "%i < z < %i m" % tuple([altgrp_bnds[2], altgrp_bnds[3]]), 
                fontsize=11, ha='center', va='bottom', transform=axset[3].transAxes)
axset[4].text(0.5, 1., "%i m < z" % altgrp_bnds[3], 
                fontsize=11, ha='center', va='bottom', transform=axset[4].transAxes)


fig.savefig("./fig_levlegqwPDFs.png")

