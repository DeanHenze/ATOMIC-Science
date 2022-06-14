# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 10:58:07 2022

@author: Dean
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



def weird_integrate(x, y, z, dx, dy):
    """
    x: 1d (size n)
    y: 1d (size m)
    z: 2d (shape mxn)
    """
    xgt0 = x > 0
    i0 = np.where(xgt0[1:] ^ xgt0[:-1])[0].item() # index of cross over from neg to pos.
    ipos = 1
    ineg = 0
    
    xx, yy = np.meshgrid(x, y)
    xyz = xx*yy*z

    x_cs = []
    cumsum = []
    cumsum_pos = []
    cumsum_neg = []
    
    # First term:
    cumsum_pos.append(np.sum(xyz[:, i0+ipos])*dx*dy)
    cumsum_neg.append(np.sum(xyz[:, i0-ineg])*dx*dy)
    cumsum.append(cumsum_pos[-1] + cumsum_neg[-1])
    ipos += 1
    ineg += 1
    x_cs.append(x[i0+ipos])

    # Remaining terms:        
    for i in range(int(np.floor(len(x)/2))-3):
        cumsum_pos.append(np.sum(xyz[:, i0+ipos])*dx*dy)
        cumsum_neg.append(np.sum(xyz[:, i0-ineg])*dx*dy)
        cumsum.append(cumsum[-1] + cumsum_pos[-1] + cumsum_neg[-1])
        
        x_cs.append(x[i0+ipos])
        ipos += 1
        ineg += 1
        
    return x_cs, cumsum_neg, cumsum_pos, cumsum



altgroups = [1,2,3,4,5]
altgrp_bnds = [300, 800, 1500, 2500]
cldgroups = [1,2,3]



fig = plt.figure(figsize=(6.5, 8))
ax_width = 0.3
ax_height = 0.25
ax1 = fig.add_axes([0.15, 3*0.075 + 2*ax_height, ax_width, ax_height])
ax2 = fig.add_axes([2*0.15 + ax_width, 3*0.075 + 2*ax_height, ax_width, ax_height])
ax3 = fig.add_axes([0.15, 2*0.075 + ax_height, ax_width, ax_height])
ax4 = fig.add_axes([2*0.15 + ax_width, 2*0.075 + ax_height, ax_width, ax_height])
ax5 = fig.add_axes([(1-ax_width)/2, 0.075, ax_width, ax_height])
axset = [ax1, ax2, ax3, ax4, ax5]


fig1 = plt.figure()
ax1 = fig1.add_axes([0.15, 0.15, 0.8, 0.8])
fig2 = plt.figure()
ax2 = fig2.add_axes([0.15, 0.15, 0.8, 0.8])


for naltgrp in altgroups:
    f = [f for f in fnames_wqPDFs 
         if "_altgrp%i_cldgrp%i" % tuple([naltgrp, 3]) in f]
    f = f[0]
    pdf = pd.read_csv(os.path.join(path_wqPDFs, f), header=0, index_col=0)
    if len(pdf)==0: continue
    
    x = pdf.columns.astype(float)
    y = pdf.index.values
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    z = pdf.values
    x_cs, cumsum_neg, cumsum_pos, cumsum = weird_integrate(x, y, z, dx, dy)
    
    ax1.plot(x_cs, cumsum)
    ax2.plot(
        np.append(-1*np.flip(x_cs), x_cs), 
        np.append(np.flip(cumsum_neg), cumsum_pos)
        )





"""
for altgrp, ax in zip(altgroups, axset):
    
    for cldgrp in cldgroups:
        
        f = [f for f in fnames_wqPDFs 
             if "_altgrp%i_cldgrp%i" % tuple([altgrp, cldgrp]) in f]
        f = f[0]
        pdf = pd.read_csv(os.path.join(path_wqPDFs, f), header=0, index_col=0)
        if len(pdf)==0: continue
"""   
        
    
    
    
    
    
    
    
    
    
    