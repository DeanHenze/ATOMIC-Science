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


## I/O paths and filenames
path_psdir = "./cospectra_levlegs/"
ncld='04'
fnames_pscld = [f for f in os.listdir(path_psdir) if '_cld%s' % ncld in f]
## I/O paths and filenames


## Plot power spectra
fig = plt.figure(figsize=(5, 3.5))
ax1 = fig.add_axes([0.125, 0.6, 0.35, 0.35])
ax2 = fig.add_axes([0.125+0.35+0.1, 0.6, 0.35, 0.35])
ax3 = fig.add_axes([0.125, 0.125, 0.35, 0.35])
ax4 = fig.add_axes([0.125+0.35+0.1, 0.125, 0.35, 0.35])
axset = np.array([ax1, ax2, ax3, ax4]).reshape([2,2])
#fig, axset = plt.subplots(2, 2, figsize=(6, 5))
fs = 5 # sampling frequency
varsplot = ["w'", "q'", "T'", "qD'"]

for levleg in [1, 3, 4]:
    f_levleg = [f for f in fnames_pscld if '_levleg%i' % levleg in f][0]
    ps_levleg = xr.load_dataset(path_psdir+f_levleg)

    for ax, v in zip(axset.flatten(), varsplot):
        f_smooth, P_smooth = specsmoother.spec_smoother(ps_levleg[v+v], fs, Rdi=0.08)
        #f_smooth, P_smooth = specsmoother.spec_smoother(ps_levleg[v+"w'"], fs, Rdi=0.08)
        f_smooth, P_smooth = f_smooth[1:], P_smooth[1:] # Remove 0 frequency.
        altmean_label = str(int(np.round(ps_levleg['alt_mean']/10)*10)) + " m"
        ax.plot(f_smooth, f_smooth*P_smooth, label=altmean_label)
## Plot power spectra
        


## -5/3 power law compute and plot with spectra
powerlaw = lambda x, p, A: A*x**p
A_plaw = lambda x, y, p: y/x**p
    # Use f, f*P pairs taken empirically from the spectra to get param A for the 
    # power law:
f_emp = np.array([0.3, 0.1, 0.1, 0.1])
fP_emp = np.array([0.05, 0.01, 0.002, 1.8*10**-10])
A = A_plaw(f_emp, fP_emp, -2./3)
    # Compute power and plot:
f = np.arange(0.02, 3, 0.1)
for Ai, ax in zip(A, axset.flatten()):
    fP = powerlaw(f, -2./3, Ai)
    ax.plot(f, fP, 'k-', label='-2/3 slope')
## -5/3 power law compute and plot with spectra



## Axes scale, limits, labels, ....
for ax, v in zip(axset.flatten(), varsplot): 
    ax.loglog()
    ax.text(
        0.95, 0.95, v, fontsize=12, 
        ha='right', va='top', transform=ax.transAxes
        )


#axset[1,0].set_xlabel('frequency (Hz)', fontsize=12)
#axset[1,1].set_xlabel('frequency (Hz)', fontsize=12)
#axset[0,0].set_ylabel(r'f*PSD (V$^2$)', fontsize=12)
#axset[1,0].set_ylabel(r'f*PSD (V$^2$)', fontsize=12)
fig.text(0.5, 0.05, 'f (Hz)', ha='center', va='center', fontsize=14)
fig.text(0.025, 0.5, r'f*PSD (V$^2$)', ha='center', va='center', rotation=90, fontsize=14)

axset[0,1].legend(loc='lower left', fontsize=8, handlelength=1)
## Axes scale, limits, labels, ....

fig.savefig("./fig_powerspectra.png")





