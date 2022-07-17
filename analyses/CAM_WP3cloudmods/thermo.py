# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:18:05 2022

@author: Dean

Collection of functions to compute atmospheric / thermodynamic quantities.
"""



import numpy as np


## Physical constants and parameters:
#------------------------------------------------------------------------------
g = 9.8 # gravitational acceleration [m/s2].
cp_d = 1004.7 # specific heat of dry air at constant pressure [J/K/kg].
cp_v = 1846.1 # specific heat of water vapor at constant pressure [J/K/kg].
cl = 4218 # specific heat of liquid water [J/K/kg].
P0 = 1000 # Reference surface pressure [mbar].
k = 0.286 # ratio of dry gas constant over dry air specific heat (pressure=const)
es0 = 6.11 # reference saturation vapor pressure [mbar].
T0 = 273.15 # reference temp for es0.    
Rd = 287.1 # specific gas constant for dry air, [J/K/kg].
Rv = 461.5 # specific gas constant for water vapor, [J/K/kg].
Lv_T0 = 2.501*10**6 # latent heat of vaporization [J/kg] at T0=273.15.
zs_2017 = 700 # African continent mean surface land elevation (just a guess) [m].
gamma = 0.01 # Dry adiabatic lapse rate.
R_D_SMOW = 155.76*(10**-6) # HDO/H2O ratio of standard ocean water.
eps=0.622 # Ratio dry air gas constant over water vapor gas constant.
#------------------------------------------------------------------------------



def air_density(T, P):
    """
    Air density from temperature (K) and pressure (hPa), assuming a dry 
    air specific gas constant.
    """
    return (100*P)/(Rd*T)



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



def ptemp(T, P):
    """
    Calculate potential temperature give temperature T (K), and pressure P (hPa).
    Both T and P are either scalars or numpy arrays.
    """
    return T*(P0/P)**k



def vtemp(T, q):
    """
    Calculate virtual temperature from temperature T (K), and mass mixing ratio 
    q (mass unitless). T and q can be either scalars or numpy arrays.  
    """
    return T*(1 + (q/eps))/(1+q)



def vptemp(T, q, P):
    """
    Calculate virtual potential temperature from temperature T (K), mass mixing ratio 
    q (mass unitless), and pressure P (hPa). T, q, and P can be either scalars or 
    numpy arrays.  
    """
    return ptemp(vtemp(T, q), P)



def satvars(T,P):
    """
    Calculate saturation vapor pressure, specific humidity, and mass mixing ratio,
    given temperature and pressure. Returns a dictionary of vapor pressure in mbar, 
    unitless mass specific humidity, and unitless mixing ratio.
    
    Inputs: T, P: numpy scalars or arrays of floats for temperature (K) and pressure 
    (mbar).
    """
    esat = es0*np.exp((Lv_T0/Rv)*(1/T0 - 1/T)) # Vapor pressure.
    qsat = eps*esat/(P-esat+eps*esat) # Specific humidity. 
    wsat = eps*esat/(P-esat) # Mass mixing ratio.
    return {'esat':esat,'qsat':qsat,'wsat':wsat}



def Fb_convert(Fb_Wm2, T, P, q):
    """
    Convert buoyancy flux from W/m2 to m2/s3. All inputs are either scalars 
    or numpy arrays with the same shape.
    
    Fb_Wm2: Buoyancy flux in units of W/m2.
    T, P: Temperature (K) and pressure (hPa)
    q: water mixing ratio (kg/kg)
    """
    rho = air_density(T, P)
    thetav = vptemp(T, q, P) # Virtual potential temperature
    return Fb_Wm2*(g/thetav)/(cp_d*rho)







