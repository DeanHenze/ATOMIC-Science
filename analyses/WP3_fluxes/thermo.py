# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:18:05 2022

@author: Dean

Collection of functions to compute atmospheric / thermodynamic quantities.
"""



import numpy as np



def P_est(alt):
    """
    Estimate pressure P (hPa) from given altitude (in km). Assume P decreases by 
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
    Latent heat flux. Density is computed from ideal gas law. 
    """
    # !!!!! Replace this with constant corrected for water vapor 
    # !!!!! concentration:
    Rd = 287.1 # specific gas constant for dry air, [J/K/kg].
    Lv_T0 = 2.501*10**6 # latent heat of vaporization [J/kg] at T0=273.15.

    rho = 100*P/(Rd*T)
    
    return flux_wq*rho*Lv_T0
