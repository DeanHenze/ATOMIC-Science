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



# Built in
import os
import itertools

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



pltchars = [
    {'colors':'red', 'fill':False, 'linewidths':1}, 
    {'colors':'grey', 'fill':True, 'extend':'max', 'linewidths':1}, 
    {'colors':'blue', 'fill':False, 'linewidths':1}
    ]

altgroups = [1,2,3,4,5]
cldgroups = [1,2,3]
fig, axset = plt.subplots(2, 3, figsize=(6.5, 5))


for altgrp, ax in zip(altgroups, axset.flatten()):
    
    for cldgrp in cldgroups:
        
        f = [f for f in fnames_wqPDFs 
             if "_altgrp%i_cldgrp%i" % tuple([altgrp, cldgrp]) in f]
        f = f[0]
        #f = fnames_wqPDFs[i+9]
        pdf = pd.read_csv(os.path.join(path_wqPDFs, f), header=0, index_col=0)
        if len(pdf)==0: continue
        cdf_levs = [0.001, 0.01] + list(np.arange(0.05, 1, 0.2))
        w = pdf.columns.values.astype(float)
        q = pdf.index.values.astype(float)
        dw = w[1] - w[0]
        dq = q[1] - q[0]
        pvals_cdf = cdf.pvals_cdflevs(pdf.values, cdf_levs, dx=dw, dy=dq)
        
        ax.contour(w, q, pdf, levels=pvals_cdf, **pltchars[cldgrp-1])


"""
for i in [0,1,2]:
    f = fnames_wqPDFs[i+9]
    pdf = pd.read_csv(os.path.join(path_wqPDFs, f), header=0, index_col=0)
    if len(pdf)==0: continue
    cdf_levs = [0.001, 0.01] + list(np.arange(0.05, 1, 0.2))
    w = pdf.columns.values.astype(float)
    q = pdf.index.values.astype(float)
    dw = w[1] - w[0]
    dq = q[1] - q[0]
    pvals_cdf = cdf.pvals_cdflevs(pdf.values, cdf_levs, dx=dw, dy=dq)
    

    
    plt.contour(w, q, pdf, levels=pvals_cdf, **pltchars[i])
"""



