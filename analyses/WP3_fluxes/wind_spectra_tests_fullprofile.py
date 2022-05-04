# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:12:26 2022

@author: Dean


Testing full profile of fluxes computed from spectral analysis for one 
stacked level leg.
"""



# Built in:
import os

# Third party:
import numpy as np
from numpy.polynomial.polynomial import polyfit
import pandas as pd
import xarray as xr
from scipy.fft import rfft, rfftfreq
from scipy import signal
from scipy.integrate import trapz
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D



varkeys = ["T'", "q'","qD'"]

def cov(x, y):
    """
    Returns covariance of x and y. x and y must be same length.
    """
    xp = x-np.nanmean(x)
    yp = y-np.nanmean(y)
    return np.nansum(xp*yp)/len(x)


#fluxes = pd.DataFrame(
#    np.zeros([7,6]).T, 
#    columns=[vlab+"w'" for vlab in varkeys]
#    )
 
def fluxes_singlelev(data, varkeys=["T'", "q'","qD'","w'"]):
    """
    ...
    
    Inputs
    ------
    data: xarray.Dataset.
        Data for a single P-3 level leg time series at 5 Hz. 

    Returns
    -------
    fluxes: dictionary 
        Fluxes for temperature "T'", water "q'", and HDO "qD'".
    alt: scalar, float
        Mean altitude of leg.
    """
    fluxes = {}   
        
    for vk in varkeys:
        
        x1 = data[vk]
        y1 = data["w'"]
        fs = 5     # Sampling frequency
        
        # Power spectrum:
        dt_window=3*60 # time interval of window for spectral decomp, seconds
        f, P = signal.csd(    # Cross-spectrum
            x1, y1, fs=fs, nperseg=int(dt_window*fs), 
            window='hamming', detrend='linear'
            )
        P = P.real  # Co-spectrum
        f, P = f[1:], P[1:] # Remove 0 frequency.
        
        # Flux:
        cov_ps_5Hz = trapz(P, x=f) # From integral of power spectrum.
        fluxes[vk] = cov_ps_5Hz
        
    # Mean altitude of leg:
    altmean = data['alt'].mean().values
        
    return fluxes, altmean



def ptemp(T, P, P0=1013):
    """
    Calculate potential temperature give temperature T (K), and pressure P (hPa).
    Both T and P are either scalars or numpy arrays.
    """
    k = 0.286 # ratio of dry gas constant over dry air specific heat (pressure=const)
    return T*(P0/P)**k


def P_est(alt):
    """
    Estimate pressure P (hPa) from given altitude. Assume P decreases by 
    100 hPa per km.
    """
    Psfc = 1013 # assumed surface pressure
    return Psfc - 100*alt
    


def flux_sh(flux_wT, T, P):
    """
    Compute sensible heat flux. Density is computed from ideal gas law. 
    Pressure is estimated as 1013 hPa surface and 100 hPa decrease per km.

    Parameters
    ----------
    flux_wT : TYPE
        DESCRIPTION.
        
    alt: float
        Altitude in km.

    Returns
    -------
    None.

    """
    # !!!!! Replace this with constant corrected for water vapor 
    # !!!!! concentration:
    Rd = 287.1 # specific gas constant for dry air, [J/K/kg].
    cp_d = 1004.7 # specific heat of dry air at constant pressure [J/K/kg].

    rho = 100*P/(Rd*T)
    
    return flux_wT*rho*cp_d



def flux_lh(flux_wq, T, P):
    """
    
    """
    # !!!!! Replace this with constant corrected for water vapor 
    # !!!!! concentration:
    Rd = 287.1 # specific gas constant for dry air, [J/K/kg].
    Lv_T0 = 2.501*10**6 # latent heat of vaporization [J/kg] at T0=273.15.

    rho = 100*P/(Rd*T)
    
    return flux_wq*rho*Lv_T0


    
def flux_profiles(date, ncldmod):
    """
    Profiles.
    """
    
    ## All 5Hz filenames for this cloud module:
    fnames = os.listdir("./levlegdata_5hz/")
    fnamehead = "WP3_5hz_%i_cld%i_" % tuple([date, ncldmod])
    fnames_cldmod = [f for f in fnames if f.startswith(fnamehead)]
    fnames_cldmod = ["./levlegdata_5hz/"+f for f in fnames_cldmod]
    
    
    ## Get fluxes for each leg of a cloud module for flight on 1/31:
    fluxes = []
    alt = []
    T = []
    q = []
    qD = []
    #levleg_str = ['levleg'+str(i) for i in range(1,7)]
    #for lls in levleg_str:
    for f in fnames_cldmod:
        #levlegdata = xr.load_dataset('levlegdata_5hz/P-3_20200131_%s_5Hz.nc' % lls)
        levlegdata = xr.load_dataset(f)
        flx, z = fluxes_singlelev(levlegdata, varkeys=["T'", "q'","qD'","w'"])
        fluxes.append(flx)
        alt.append(z)
        T.append(levlegdata["T"].mean().values)
        q.append(levlegdata["q"].mean().values)
        qD.append(levlegdata["qD"].mean().values)
        levlegdata.close()
        
    # Merge all dictionaries:
    fluxes = {k: [d[k] for d in fluxes] for k in fluxes[0]}
    fluxes = {k: np.array(fluxes[k]) for k in fluxes}
    
    alt = np.array(alt)
    T = np.array(T) - 3 # -3 fudge term since temperatures seem too high.
    q = np.array(q)
    qD = np.array(qD)
    
    
    P = P_est(alt/1000) # input altitude in km.
    theta = ptemp(T, P)
    
    fig, axes = plt.subplots(1, 4, figsize=(10, 6))
    twinaxes = [ax.twiny() for ax in axes]
    axes[0].scatter(
        flux_sh(fluxes["T'"], T, P), alt
        )
    twinaxes[0].plot(theta, alt)
    
    axes[1].scatter(flux_lh(fluxes["q'"]/1000, T, P), alt) # Divide by 1000 for kg/kg.
    twinaxes[1].plot(q, alt)
    axes[3].scatter(fluxes["w'"], alt)
    
    axes[0].set_xlabel("SH flux (W/m2)", fontsize=14)
    twinaxes[0].set_xlabel(r"$\theta$ (K)", fontsize=14)
    axes[1].set_xlabel("LH flux (W/m2)", fontsize=14)
    twinaxes[1].set_xlabel("q (g/kg)", fontsize=14)
    axes[2].set_xlabel(r"$\delta D_{flux}$ (permil)", fontsize=14)
    twinaxes[2].set_xlabel(r"$\delta D$ (permil)", fontsize=14)
    axes[3].set_xlabel(r"w'w' (m/s)$^2$", fontsize=14)
    
    
    RD_SMOW = 155.76*(10**-6) # HDO/H2O ratio of standard ocean water.
    dD_fluxes = ((fluxes["qD'"]/fluxes["q'"])/RD_SMOW - 1)*1000
    dD_profile = ((qD/q)/RD_SMOW - 1)*1000
    axes[2].scatter(dD_fluxes, alt)
    twinaxes[2].plot(dD_profile, alt)
    twinaxes[2].set_xlim(-200, 0)
    axes[2].set_xlim(-200, 0)
    axes[3].set_xlim(0, 0.55)



#date = 20200124
date = 20200205
ncldmod=2
flux_profiles(date, ncldmod)


"""
from numpy.polynomial.polynomial import polyfit
slopes = []
fig, axes_qdD = plt.subplots(2, 3, figsize=(8,8))
for lls, ax in zip(levleg_str, axes_qdD.flatten()):
    levlegdata = xr.load_dataset('levlegdata_5hz/P-3_20200131_%s_5Hz.nc' % lls)
    
    picdata = levlegdata[["q","qD"]].to_dataframe()
    picdata.dropna(how='any', axis=0, inplace=True)
    data_dD = ((picdata["qD"]/picdata["q"])/RD_SMOW - 1)*1000
    ax.scatter(picdata["q"], picdata["q"]*data_dD)
    results = polyfit(picdata["q"], picdata["q"]*data_dD, 1)

    #dD_levleg = ((levlegdata["qD"]/levlegdata["q"])/RD_SMOW - 1)*1000
    #ax.scatter(levlegdata["q"], levlegdata["q"]*dD_levleg)
    #levlegdata.close()
    #results = polyfit(levlegdata["q"], levlegdata["q"]*dD_levleg, 1)
    slopes.append(results[1])
    
    
    
dDresults_df = pd.DataFrame({'flux_ps':dD_fluxes,'flux_slope':slopes}, 
                            index=alt)
print(dDresults_df)
"""



