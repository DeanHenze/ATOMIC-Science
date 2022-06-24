# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 16:17:19 2022

@author: Dean

Methods section figure of power spectra for w', T', q', and qD'.
"""


# Built in
import os

# Third party
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import specsmoother


## I/O
path_psdir = "./cospectra_levlegs/"
ncld='04'
fnames_pscld = [f for f in os.listdir(path_psdir) if '_cld%s' % ncld in f]


fig, axset = plt.subplots(2, 2, figsize=(6, 5))
fs = 5 # sampling frequency
varsplot = ["w'", "q'", "T'", "qD'"]


for levleg in [1, 3, 4]:
    f_levleg = [f for f in fnames_pscld if '_levleg%i' % levleg in f][0]
    ps_levleg = xr.load_dataset(path_psdir+f_levleg)

    for ax, v in zip(axset.flatten(), varsplot):
        f_smooth, P_smooth = specsmoother.spec_smoother(ps_levleg[v+v], fs, Rdi=0.08)
        f_smooth, P_smooth = f_smooth[1:], P_smooth[1:] # Remove 0 frequency.
        altmean_label = str(int(np.round(ps_levleg['alt_mean']/10)*10)) + " m"
        ax.plot(f_smooth, f_smooth*P_smooth, label=altmean_label)
        


for ax, v in zip(axset.flatten(), varsplot): 
    ax.loglog()
    ax.text(
        0.95, 0.95, v, fontsize=12, 
        ha='right', va='top', transform=ax.transAxes
        )


axset[1,0].set_xlabel('frequency (Hz)', fontsize=12)
axset[1,1].set_xlabel('frequency (Hz)', fontsize=12)
axset[0,0].set_ylabel(r'f*PSD (Hz*V$^2$)', fontsize=12)
axset[1,0].set_ylabel(r'f*PSD (Hz*V$^2$)', fontsize=12)

axset[0,0].legend(loc='upper left', fontsize=10, handlelength=1)


fig.savefig("./fig_powerspectra.png")





