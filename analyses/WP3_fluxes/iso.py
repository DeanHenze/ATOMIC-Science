import numpy as np



### Physical constants and parameters:
#------------------------------------------------------------------------------
g = 9.8 # gravitational acceleration [m/s2].
cp = 1005 # specific heat of dry air at constant pressure [J/K/kg].
Ps_ref = 1013 # Reference surface pressure [mbar].
    # This temperature is needed to calculate a P(z). It is the assumed 
    # temperature at Ps_ref. For this I've just assumed that in 
    # the SEA, average surface temperature and pressure are 295K and 1013
    # mbar.
Ts_ref=295
k = 0.286 # ratio of dry gas constant over dry air specific heat (pressure=const)
es0 = 6.11 # reference saturation vapor pressure [mbar].
T0 = 273 # reference temp for es0.    
Rv = 461 # specific gas constant for water vapor, [J/K/kg].
L = 2.5*10**6 # latent heat of vaporization [J/kg].
zs_2017 = 700 # African continent mean surface land elevation (just a guess) [m].
gamma = 0.01 # Dry adiabatic lapse rate.
R_D_SMOW = 155.76*(10**-6) # HDO/H2O ratio of standard ocean water.
R_18_SMOW = 2005.2*(10**-6) # 18O/16O ratio of standard ocean water.
#------------------------------------------------------------------------------



### Method to find the array index of the element which is closest to some 
### asked for value:
###____________________________________________________________________________
def i_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
###____________________________________________________________________________



###____________________________________________________________________________
"""
Convert an input humidity from units of ppmv to g/kg or vise versa. 
Inputs:
    units: either 'gkg' or 'ppmv'; the units returned by the fxn.
    q (float) is the humidity.
"""
###____________________________________________________________________________
def ppmv_gkg_convert(units, q):
    eps=0.622
    if units=='gkg': return q*(eps/1000)
    if units=='ppmv': return q*(1000/eps)
###____________________________________________________________________________
    
    
    
###____________________________________________________________________________
"""
Convert between HDO mixing ratio, qD, and HDO delta notation, dD [permil].

Inputs:
    conversion: Either 'qD2dD' to convert from HDO mixing ratio to dD, or 
        'dD2qD' for the reverse.
    a1: Numpy array, any shape. Holds values for humidity used in the conversion
        formulas. 
    a2: Numpy array, same shape as a1. If conversion is set to 'qD2dD' then a2 
        should be assigned qD values (same units as a1). If conversion is set to 'dD2qD', a2 is
        assigned dD values.
"""
###____________________________________________________________________________
def qD_dD_convert(conversion, a1, a2):
    if conversion=='qD2dD':
        RD = a2/a1 # Deuterium ratio qD/q.
        return 1000*(RD/R_D_SMOW - 1) # Formula for dD.    
    
    elif conversion=='dD2qD':
        RD = R_D_SMOW*(a2/1000 + 1) # Deuterium ratio computed from dD.
        return RD*a1 # qD = RD*q.
###____________________________________________________________________________

    

#______________________________________________________________________________
"""
Calc equilibrium fractionation coefficient for either D/H or 18O/16O. The 
coefficient is computed from the passed temperature(s)
 
Inputs:
    - isotopologue (str): Either 'D' or '18O'.
    - T (scalar or array-like): Temperature(s) in kelvin.
"""
#______________________________________________________________________________
def alpha_e(isotopologue, T):
    T=np.array(T)
    
    if isotopologue=='D':
        return np.exp( 1158.8*(T**3/10**12) - 1620.1*(T**2/10**9)
                      + 794.84*(T/10**6) - 161.04/10**3 + 2.9992*(10**6/T**3) )
    
    if isotopologue=='18O':
        return np.exp( -7.685/10**3 + 6.7123/T - 1.6664*(10**3/T**2)
                      + 0.35041*(10**6/T**3) )
#______________________________________________________________________________



#______________________________________________________________________________
"""
Return dD value for vapor in equilibrium with standard ocean water of a 
specified temperature (T) and dD of the liquid (dD_l). 
T: scalar or numpy array. 
dD: scalar or numpy array.
"""
#______________________________________________________________________________
def dD_eqvapor(T, dD_l = 0, isotopologue='D'):
    ae = alpha_e(isotopologue, T) # Equilibrium fractionation factors.
    return 1000*( (dD_l + 1000)/(1000*ae) - 1 )
    #return 1000*(1./ae - 1)
#______________________________________________________________________________



#______________________________________________________________________________
"""
Calculate theoretical MBL dD and d18O values from the CGE with MJ79 closure.

The Craig-Gordon equation (CGE) is:
    Re = a_k.[(Roc/a_e) - h.Ra]/(1-h)
    
where Re, Roc, and Ra are the isotopeic ratios of evaporation, ocean, and the
atmosphere in which the ocean is evaporating, respectively; a_e and a_k are 
the equilibrium and kinetic fractionations; h is the relative humidity of the 
atmsophere normalized to the sea surface temperature.

I think the Merlivat and Jouzel 1979 (MJ79) closure to the CGE sets Ra=Re. 
From this I derived:
    Re = (1/a_e).[Roc/(q + h.(a_k-1))]
 
a_k for deuterium will be taken from Merlivat and Jouzel as 1.0053 for smooth and 
1.0035 for rough regimes. Temporarily I cannot find a kinetic factor for 18O, 
so I am deriving it from the MJ79 results that (1) a_k = 1 - k, and 
(2) k_D/k_18O = 0.88. The 'regime' input is a string, set to either 'smooth' 
or 'rough', which determines what values of kinetic fractionation factors to 
use. It's default is 'smooth'.

Inputs:
    sst (float or array): sea surface temperature, used to calc a_e.
    
    h (float or array, same size as sst): relative humidity (as a decimal, not 
      %) of air with respect to the sst.
    
    dD_oc, d18O_oc (floats or numpy arrays same size as sst): optional values 
        for the isotope ratios of the ocean. Otherwise defaults to 0 for both 
        isotopologues.
        
    regime (string, default='smooth'): 'smooth' or 'rough' to use different 
        values of the kinetic fractionation factor.
        
Returns:
    Two element list. First is a numpy float or array of dD values, second is
    a float or array of d18O values.
"""
#______________________________________________________________________________
def delta_evap_MJ79(sst, h, dD_oc=0, d18O_oc=0, regime='smooth'):
    
# get np.arrays:
    sst=np.array(sst); h=np.array(h)
    
# Kinetic fractionation factors:
    if regime=='smooth':
        a_k_D = 1.0053 # smooth regime
    elif regime=='rough':
        a_k_D = 1.0035 # rough regime
    # Calc a_k_18O by combining (1) and (2) in the header notes:
    # Also found in Benetti et al 2015 that a_k_18 = 1.006 for smooth conditions.
    a_k_18 = 1 - ( (1-a_k_D)/0.88 )
    
# Calculate equilibrium fractionation factors from SST:
    #a_e_D = np.exp( 1158.8*(sst**3/10**12) - 1620.1*(sst**2/10**9)            \
    #                + 794.84*(sst/10**6) - 161.04/10**3 + 2.9992*(10**6/sst**3) )
    a_e_D = alpha_e('D', sst)
    a_e_18 = np.exp( -7.685/10**3 + 6.7123/sst - 1.6664*(10**3/sst**2)        \
                     + 0.35041*(10**6/sst**3) )
        
# Convert ocean isotope ratios from delta notation to regular fractions:
    R_D_oc = R_D_SMOW*(dD_oc/1000 + 1)
    R_18O_oc = R_18_SMOW*(d18O_oc/1000 + 1)
    
# Calculate MJ79 prediction of MBL isotope ratios and delta values:
    R_D = (R_D_oc/a_e_D)/(h + a_k_D*(1-h))
    dD = 1000*(R_D/R_D_SMOW - 1)
    R_18 = (R_18O_oc/a_e_18)/(h + a_k_18*(1-h))
    d18O = 1000*(R_18/R_18_SMOW - 1)

    return dD, d18O
#______________________________________________________________________________   



###____________________________________________________________________________
"""
Estimate lifting condensation level from given surface temp and humidity.
qs, Ts, and zs are mean African continent surface humidity [unitless],
temperature [Kelvin], and elevation [meters]. Return lifting condensation
level [meters] and temperature [Kelvin]. All inputs are floats.
"""
###____________________________________________________________________________
def LCL(qs, Ts, zs):
    z=np.arange(0,20000,100) # Height array [m] to generate P(z).
    P = Ps_ref*( 1 - (g*z)/(cp*Ts_ref) )**(1/k) # Pressure profile P(z).
    # Dewpoint temp profile:
    Td = ( (1/T0) - (Rv/L)*np.log(qs*P/es0) )**-1
    # Dry adiabatic lapse profile over land surface:
    Tz = Ts - gamma*(z-zs)
    # Find index of point-of-intersection between T_d and Tz:
    i_lcl = np.argmin(np.abs(Td - Tz))
    return z[i_lcl], Tz[i_lcl]
###____________________________________________________________________________
    
    

###____________________________________________________________________________
    """
    Estimate humidity of air starting at height zi [meters] and ending
    at height zf [meters], assuming a moist adiabatic (ma) profile and total removal of 
    excess water. Starting temperature of air is Ti=T(zi) [K]. All inputs are floats.
    """
###____________________________________________________________________________
def q_ma(zi, Ti, zf):
    # Temperature at zf assuming moist adiabatic lapse from the LCL upward:
    T_zf = Ti - 0.006*(zf-zi) 
    # Pressure at zf:
    Pf = Ps_ref*( 1 - (g*zf)/(cp*Ts_ref) )**(1/k)        
    # Saturation mixing ratio at T(zf):
    esat = es0*np.exp((L/Rv)*(1/T0 - 1/T_zf)) # [mbar].
    qsat = esat/Pf # Convert to mixing ratio.
    return qsat
###____________________________________________________________________________
  
   

###____________________________________________________________________________
"""
For range of humidity values q [ppmv], estimate dewpoint temperatures assuming
a typical profile P(z), dry adiabatic lapse rate up to an LCL, and a moist 
adiabatic lapse rate from the LCL upward:
    
Inputs:
-------
q: array of floats; humidity values.
Ts, qs: floats; surface temperature [K] and humidity [unitless]. 
zs: float; surface elevation [meters].
"""
###____________________________________________________________________________
def T_sat(q, Ts, qs, zs):

# Remove ppmv units of qs and q:
    q = q/10**6
    qs = qs/10**6
    
# Generate pressure profile:
    dz=100
    zf=20000
    Pz = Ps_ref*( 1 - (g*np.arange(zs,zf,dz))/(cp*Ts_ref) )**(1/k)
    
# Generate temperature profile which goes dry adiabatic up to LCL, then
# moist adiabatic upwards of LCL:
    z_lcl, T_lcl = LCL(qs, Ts, zs) # Get LCL.
    T_da = Ts - 0.01*(np.arange(zs,z_lcl,dz)-zs) # Dry adiabatic profile.
    T_ma = T_lcl - 0.006*(np.arange(z_lcl,zf,dz)-z_lcl) # Moist adabatic.
    Tz = np.append(T_da, T_ma) # Combine dry and moist profiles.
    
# For T(z) and P(z), calc saturation humidity:
    esat = es0*np.exp((L/Rv)*(1/T0 - 1/Tz)) # [mbar].
    qsat = esat/Pz # Convert to mixing ratio.  
    
# For each value of q, find the height at which q ~ qsat, and get temperature
# at that height:
    T_return = np.zeros(len(q)) # Temperature values to return.
    for i in range(len(q)):
        i_sat = np.argmin(np.abs(qsat-q[i])) # Index of closest saturation humidity.
        T_return[i] = Tz[i_sat]
        
    return T_return
###____________________________________________________________________________



###____________________________________________________________________________
"""
dD vs q curve for Rayleigh model

Inputs:
    q0, d0: floats; initial humidiy [ppmv] and delta-deuterium [permil] of air parcel.
    units: Either 'ppmv' or 'gkg'. The units to plot humidity in. 
    q: array-like, floats; Domain of q over which to calc dD.
    ax: matplotlib.pyplot axes instance; axes to plot on. Default is
        None, in which case q and dD arrays are returned. Else q,dD are plotted
        on the specified axes and no arrays are returned.     
    get_xtra: bool, default false. If true, return temperature T used to calc
        fractionation at each step.
"""
###____________________________________________________________________________
def Rayleigh_curve_dD(q0, d0, q, ax=None, units=None, get_xtra=False):
    
# Estimate dewpoint temperatures at each q, for assumed values of surface humidity,
# temperature, and elevation:
    Ts=298; qs=25000; zs=0
    T = T_sat(q, Ts, qs, zs)
# Equilibrium fractionation factor at temperatures T:
    a_e_D = alpha_e('D', T)
# Get R0 from d0:
    R0 = (d0/1000 + 1)*R_D_SMOW
# Calc R = qD/q from Rayleigh equation:
    R = R0*(q/q0)**(a_e_D - 1)
# Convert R to dD:
    dD = 1000*(R/R_D_SMOW - 1)
# Units of q, only need to do something if gkg is chosen:
    if units=='gkg':q = q*(0.622/1000)
# Plot or return arrays:
    if ax==None:
        if get_xtra: return dD, T
        return dD
    else:
        ax.plot(q,dD,'k-', linewidth=2.5)  
###____________________________________________________________________________



###____________________________________________________________________________
"""
Computs dD for moist adiabatic process (e.g. no condensate removal). 
Equation taken from Noone 2012.

q0, dD0: floats.
    Starting values of q (ppmv) and dD (permil) for model.
q: float or numpy array.
    Compute dD for these values of q (ppmv).
"""
###____________________________________________________________________________
def moist_adiabatic(q0, dD0, q):

    # Estimate dewpoint temperatures at each q, for assumed values of surface humidity,
    # temperature, and elevation:
    Ts=298; qs=25000; zs=0
    T = T_sat(q, Ts, qs, zs)
    # Equilibrium fractionation factor at temperatures T:
    ae_D = alpha_e('D', T)    

    return dD0 + (ae_D - 1)*(q - q0)/(q - ae_D*(q - q0))
###____________________________________________________________________________

        
        
###____________________________________________________________________________
"""
dD vs q mixing curve between two given mixing end members.
    
equations for mixing are:
    q_mix = f*q1 + (1-f)*q2 => f = (q_mix - q2)/(q1 - q2)
    R_mix = (f*q1_D + (1-f)*q2_D)/q_mix
where subscripts 1 denotes quantities for the moist mixing end member and 2 denotes
quantities for the dry end member. q1_D and q2_D are deuterium mixing ratios.

Inputs:
    q1, dD1: floats; the humidity [any units] and delta deuterium [permil] values of the moist mixing 
        end member. Note dD1 is the delta value, not the mixing ratio which shows
        up in the above equation.
    q2, dD2: floats; analogous quantities for the dry end member.
    q: array-like, floats; Domain of q over which to calc dD.
    ax: matplotlib.pyplot axes instance; axes to plot on. Default is
        None, in which case fxn returns an array of calc'd dD values from the 
        mixing model. Else q,dD are plotted on the specified axes and no arrays 
        are returned.  
    linewidth: integer, default 4. linewidth argument to pass to matplotlib if plotting.
    return_f: bool, default=False. If set to true, will additionally return
        values of the mixing fraction f used along the mixing curve.
"""
###____________________________________________________________________________
def mix_curve(q1, dD1, q2, dD2, q, ax=None, linewidth=4, return_f=False):
    
    R_D_SMOW = 155.76*(10**-6) # HDO/H2O ratio of standard ocean water.
    
    f = (q - q2)/(q1 - q2) # Calc mixing fraction for each q. 
    
    # Deuterium mixing ratio end members from the delta values:
    q1_D = (dD1/10**3 + 1)*R_D_SMOW*q1
    q2_D = (dD2/10**3 + 1)*R_D_SMOW*q2
    
    # Calculate dD values for each f:
    R_D_mix = (f*q1_D + (1-f)*q2_D)/q
    dD_mix = 10**3*(R_D_mix/R_D_SMOW - 1)
    
    # Return or plot:
    if ax==None:
        if return_f: return dD_mix, f
        return dD_mix
    else:
        ax.plot(q, dD_mix, color='k', linestyle='dashed', alpha=1, linewidth=linewidth)  
###____________________________________________________________________________