# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:47:26 2022

@author: Dean

Figure for vertical profiles of <w'w'>, sensible and latent heat fluxes, 
and dD of the flux. 
"""



import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

# Local code
import profileplotter
import rangescaler
import prfrestructure
import thermo_deSzoeke
import oversampler



## Data directories and filenames.
##_____________________________________________________________________________    
path_keyaltstable = "../WP3_cloudmodule_char/cldmod_keyaltitudes.csv"

dir_flux = "./fluxes_levlegs/"
fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]

path_rhbflux = ("../../data/RHB/metflux/EUREC4A_ATOMIC_RonBrown_10min"
                "_nav_met_sea_flux_20200109-20200212_v1.3.nc"
                )
##_____________________________________________________________________________
## Data directories and filenames.



def get_fluxprofiles(ncld_list, keyalts_table, dir_flux, 
                     scale_altkeys=[]):
    """
    Inputs
    ------
    keyalts_table: pandas.DataFrame.
        If scaling altitude by a subset of columns in this table, those 
        columns cannot have NAN's.
        
    scale_altkeys: iterable of str's.
        Optional. Keys in keyalts_table to scale altitude by. If not empty,
        altitudes will be mapped as:
            0 <-- 0
            1 <-- 1st key in scale_alt
            2 <-- 2nd key in scale_alt
            3 ...
    """

    keyalts_table = keyalts_table.set_index('ncld')


    # All flux filenames:
    fnames_flux = [f for f in os.listdir(dir_flux) if f.endswith('.nc')]


    # Collect flux profiles into a dictionary of pandas.DataFrame's.
    # Each DataFrame will contain all profiles for a specific var.
    fluxvarkeys = ["u'u'_bar", "v'v'_bar", "w'w'_bar", "TKE", "TKE_h",
                   "flux_sh", "flux_lh", "flux_b", "dD_flux", 
                   "q'w'_bar", "T'w'_bar", "anisotropy_ratio", "ratio_Fb_sensible"
                   ]
    fluxprfs_dict = {}
    for k in fluxvarkeys: fluxprfs_dict[k] = pd.DataFrame({})
    
    
    for n in ncld_list:
        
        print("Working on cloud module %i." % n)
        
        # Key altitude info for this module:
        keyalts = keyalts_table.loc[n]
        
        # Load flux data:
        fname_flux = [f for f in fnames_flux if "_cld%s_" % str(n).zfill(2) in f]
        if len(fname_flux) == 0: continue
        fname_flux = fname_flux[0]        
        fluxes = xr.load_dataset(dir_flux + fname_flux)
    
        # Optional scale altitude.
        # 0=sea level.
        if len(scale_altkeys) != 0:
            alt_scalepoints = [0] + [keyalts[k] for k in scale_altkeys]
            alt = rangescaler.piecewise_linscale(
                fluxes['alt'].values, 
                alt_scalepoints, np.arange(len(alt_scalepoints))
                )
        else:
            alt = fluxes['alt'].values
            
        # Merge flux profiles of current module into resp. DataFrame's:
        for k in fluxvarkeys:
            
            fluxvar_pd = pd.DataFrame(  # easier to work in pandas.
                {k+"_cld%i"%n: fluxes[k].values}, 
                index=alt
                )
                
            fluxprfs_dict[k] = pd.merge(
                fluxprfs_dict[k], fluxvar_pd, 
                how='outer', left_index=True, right_index=True
                )   
            
            
    return fluxprfs_dict




def plot_RHBmeanfluxes(path_rhbflux, ax_sh, ax_lh, ax_bflux):
    """
    Plot mean, std of RHB surface fluxes during IOP portion with P-3 flights. 
    Plot SHF and LHF on ax_sh and ax_lh repspectively.
    """    
    rhb_metflux = xr.load_dataset(path_rhbflux)
    
    # Remove measurements for dates before the P-3 started sampling:
    dts_wp3start = np.datetime64("2020-01-17")
    rhb_metflux = rhb_metflux.where(rhb_metflux['time'] > dts_wp3start, drop=True)

    # Compute buoyancy flux:
    flux_b = thermo_deSzoeke.buoyancy_flux(
        rhb_metflux['hs_bulk'], rhb_metflux['hl_bulk'], 
        T=rhb_metflux['tair'], p=rhb_metflux["psealevel_ship"]*100, 
        q=rhb_metflux["qair"]/1000
        )
    rhb_metflux['bflux_bulk'] = flux_b[0]

    # Means and std's:
    rhbfluxkeys = ["hs_bulk", "hl_bulk", "bflux_bulk", "ustar", "gust"]
    flux_P3iopmean = rhb_metflux[rhbfluxkeys].mean(dim='obs')
    flux_P3iopstd = rhb_metflux[rhbfluxkeys].std(dim='obs')
    
    # Error bar plot:
    rhb_altsamp = 20/500 # height of RHB instruments, roughly normalized by z_ML.
    ax_sh.errorbar(
        -1*flux_P3iopmean['hs_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hs_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )
    ax_lh.errorbar(
        -1*flux_P3iopmean['hl_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['hl_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )
    ax_bflux.errorbar(
        -1*flux_P3iopmean['bflux_bulk'], rhb_altsamp, 
        xerr = flux_P3iopstd['bflux_bulk'], 
        marker='s', markersize=8, mfc='orange', mec='black', ecolor='black', 
        zorder=99
        )



def xaxes_cleanup(axset_wind, axset_scalar):
    """
    Set x-axis labels, limits, etc. for sets of wind and flux axes.
    """

    ## Limits, labels, etc for wind / turbulence axes:    
    for ax in axset_wind[0:3]: 
        ax.set_xlim(0, 0.8)
        ax.set_xticks(np.arange(0, 1, 0.2))
    axset_wind[3].set_xlim(0, 1.2)  
    axset_wind[3].set_xticks(np.arange(0, 1.2, 0.2))
    axset_wind[0].set_xlabel(r"u'u' (m/s)$^2$", fontsize=12)
    axset_wind[1].set_xlabel(r"v'v' (m/s)$^2$", fontsize=12)
    axset_wind[2].set_xlabel(r"w'w' (m/s)$^2$", fontsize=12)
    axset_wind[3].set_xlabel(r"TKE (m/s)$^2$", fontsize=12)

    ## Limits, labels, etc for flux axes:    
    axset_scalar[0].set_xlim(-30, 35)
    axset_scalar[1].set_xlim(-100, 650)
    axset_scalar[2].set_xlim(-5e-4, 1.7e-3)
    axset_scalar[3].set_xlim(-160, -30)
    axset_scalar[0].set_xlabel(r"$F_{sh}$ (W/m$^2$)", fontsize=12)
    axset_scalar[1].set_xlabel(r"$F_{lh}$ (W/m$^2$)", fontsize=12)
    axset_scalar[2].set_xlabel(r"$F_{b}$ ($m^2/s^3$)", fontsize=12)
    axset_scalar[3].set_xlabel(r"$\delta D_{flux}$ (permil)", fontsize=12)



def yaxes_cleanup(axset_wind, axset_scalar, yticks, yticklabels):
    """
    Set y-axis labels, limits, etc. for sets of wind and flux axes.
    """
    ## Limits, labels, etc for wind / turbulence axes:    
    for ax in axset_wind:
        ax.set_ylim(min(yticks), max(yticks))
        ax.set_yticks(yticks)
    axset_wind[0].set_yticklabels(yticklabels)
    for ax in axset_wind[1:]: ax.set_yticklabels(["" for t in yticks])
    axset_wind[0].set_ylabel('altitude (m)', fontsize=12)
    
    ## Limits, labels, etc for flux axes:    
    for ax in axset_scalar:
        ax.set_ylim(min(yticks), max(yticks))
        ax.set_yticks(yticks)
    axset_scalar[0].set_yticklabels(yticklabels)
    for ax in axset_scalar[1:]: ax.set_yticklabels(["" for t in yticks])
    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)



def fig_scatter():
    """
    Create and save figures for turbulence and flux profiles where altitude 
    is scaled by LCL and trade inversion bottom.
    """
    
    keyalts_table = pd.read_csv(path_keyaltstable) # key altitudes for each cld mod.
    
    # Cloud module number groupings:
    ncld_g1 = [1, 5, 4, 8]
    ncld_g2 = [7, 9, 11, 10, 6]
    ncld_g3 = [12, 15, 3, 2, 13, 16, 14]  
    #ncld_g1 = [1, 5, 4]
    #ncld_g2 = [8, 7, 9, 11, 10, 6]
    #ncld_g3 = [15, 3, 2, 13, 16, 14, 12]  
    
    ncld_g1.sort()
    ncld_g2.sort()
    ncld_g3.sort()
    ncld_groups = [ncld_g1, ncld_g2, ncld_g3]
    
    #scale_altkeys = ["z_lcl"] # scale altitude by these quantities.
    #altbinwidth=0.5
    scale_altkeys = [] # scale altitude by these quantities.
    altbinwidth = 250

    # Get flux profiles grouped; groups are elements of a list:
    fluxprfs_grouped = []
    for cg in ncld_groups:
        fluxprfs_grouped.append(
            get_fluxprofiles(       # Returns flux profiles as a dictionary of 
                cg, keyalts_table,  # pandas.DataFrame's. 
                dir_flux, scale_altkeys=scale_altkeys
                )
            )
        
        
    ## Plot turbulence and flux profiles
    ##_________________________________________________________________________
    def scatterwithmean(data_prfs, varkey, ax, color, altbinwidth, gridalt, n):
        """
        """
        for key_prf in data_prfs.columns:
            prf = data_prfs[key_prf]
            ax.scatter(prf, data_prfs.index, color=c, s=8, alpha=0.4)
        

        # Get mean and std prfiles using oversampling:
        df = data_prfs.stack().reset_index()
        if 'level_0' in df.columns:
            ovs = oversampler.oversample_1d(df[0], df['level_0'], gridalt, ffact=0.6, return_stdev='yes')
        elif 'level_1' in df.columns:
            ovs = oversampler.oversample_1d(df[0], df['altleg'], gridalt, ffact=0.6, return_stdev='yes')
        
        # Plot mean and std of mean:
        ax.plot(ovs['mean'], gridalt, color=color)
        ax.fill_betweenx(
            gridalt, 
            ovs['mean']-ovs['stdev']/3**0.5,  # ~3 observations per mean
            ovs['mean']+ovs['stdev']/3**0.5, 
            color=color, alpha=0.2
            )
        
        return ovs['mean']
    
        

    fig_scalar = plt.figure(figsize=(6.5, 3))
    axset_scalar = (
        fig_scalar.add_axes([0.1, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.325, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.55, 0.2, 0.2, 0.75]),
        fig_scalar.add_axes([0.775, 0.2, 0.2, 0.75]),
        )
    fig_wind = plt.figure(figsize=(6.5, 3))
    axset_wind = (
        fig_wind.add_axes([0.1, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.325, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.55, 0.2, 0.2, 0.75]),
        fig_wind.add_axes([0.775, 0.2, 0.2, 0.75]),
        )
    
    
    colors = ['red', 'grey', 'blue']
    cldgrps_n = [1,2,3]
    tkeh_meanprfs = []
    anisotropy_ratio_meanprfs = []
    ratio_Fbsensible_meanprfs = []
    ww_meanprfs = []
    tke_meanprfs = []
    shf_meanprfs = []
    lhf_meanprfs = []
    bf_meanprfs = []
    dDf_meanprfs = []
    gridalt = np.arange(200, 3400, 250)

    
    save_dir = "./mean_profiles/" # Save mean profiles here
    if not os.path.isdir(save_dir): os.mkdir(save_dir)
    for fpgroup, c, n in zip(fluxprfs_grouped, colors, cldgrps_n):

        tkeh_meanprfs.append(
            scatterwithmean(fpgroup["TKE_h"], "TKE_h", axset_wind[0], c, altbinwidth, gridalt, n))
        ww_meanprfs.append(
            scatterwithmean(fpgroup["w'w'_bar"], "w'w'_bar", axset_wind[1], c, altbinwidth, gridalt, n))
        tke_meanprfs.append(
            scatterwithmean(fpgroup["TKE"], "TKE", axset_wind[2], c, altbinwidth, gridalt, n))
        anisotropy_ratio_meanprfs.append(
            scatterwithmean(fpgroup["anisotropy_ratio"], "anisotropy_ratio", axset_wind[3], c, altbinwidth, gridalt, n))
        shf_meanprfs.append(
            scatterwithmean(fpgroup["flux_sh"], "flux_sh", axset_scalar[0], c, altbinwidth, gridalt, n))
        lhf_meanprfs.append(
            scatterwithmean(fpgroup["flux_lh"], "flux_lh", axset_scalar[1], c, altbinwidth, gridalt, n))
        bf_meanprfs.append(
            scatterwithmean(fpgroup["flux_b"], "flux_b", axset_scalar[2], c, altbinwidth, gridalt, n))
        #ratio_Fbsensible_meanprfs.append(
        #    scatterwithmean(fpgroup["ratio_Fb_sensible"], "ratio_Fb_sensible", axset_scalar[3], c, altbinwidth, gridalt, n))
        dDf_meanprfs.append(
            scatterwithmean(fpgroup["dD_flux"], "dD_flux", axset_scalar[3], c, altbinwidth, gridalt, n))


    # RHB surface flux means, stds for the P-3 sampling time period:
    plot_RHBmeanfluxes(path_rhbflux, axset_scalar[0], axset_scalar[1], axset_scalar[2])
    
    
    # Ref line at x=0 for buoyancy flux:
    axset_scalar[2].vlines(
        0, -100, 3300, 
        color='black', linewidth=1.5, linestyle='dashed', zorder=100
        )
    
    # Rectangle patches at altitude ranges of median cloud tops
    axset_wind[1].plot((0.55, 0.55), (1400,2000), linewidth=3, color='red')
    axset_wind[1].plot((0.5, 0.5), (1600,2000), linewidth=3, color='grey')
    axset_wind[1].plot((0.55, 0.55), (2500,2800), linewidth=3, color='blue')
    ##_________________________________________________________________________    
    
    """
    ## Plot w' skewness
    ##_________________________________________________________________________
    # Load and merge data tables:
    wmom = pd.read_csv("./WP3_wmoments_levlegs.csv")
    keyalts_table = pd.read_csv("../WP3_cloudmodule_char/cldmod_keyaltitudes.csv")
    wmomdata = wmom.merge(
        keyalts_table, 
        left_on="ncld", right_on="ncld", how="left"
        )
    
    
    wmomprfs = []
    for ncld_list in [ncld_g1, ncld_g2, ncld_g3]:
        wmomprfs.append(
            prfrestructure.restruct_profiles(
                wmomdata, 'wskew', ncld_list, 
                keyalts_table, scale_altkeys=[]
                )
            )
        
    # Plot skewness profiles:
    wskew_meanprfs = []
    for prfs, c, n in zip(wmomprfs, colors, cldgrps_n):  
            wskew_meanprfs.append(
                scatterwithmean(prfs, 'wskew', axset_wind[3], c, altbinwidth, gridalt, n))
            
    # Ref line at x=0:
    axset_wind[3].vlines(
        0, -100, 3300, 
        color='black', linewidth=1.5, linestyle='dashed', zorder=100
        )
    ##_________________________________________________________________________
    """
    
    ## Save mean profiles as .csv files:
    ##_________________________________________________________________________
    varprfs = [tkeh_meanprfs, ww_meanprfs, tke_meanprfs, shf_meanprfs, 
               lhf_meanprfs, bf_meanprfs, dDf_meanprfs, 
               anisotropy_ratio_meanprfs, ratio_Fbsensible_meanprfs]
    varkeys = ["TKE_h", "wp2_bar", "TKE", "flux_sh", 
               "flux_lh", "flux_b", "dD_flux",
               "anisotropy_ratio", "ratio_Fbsensible"] 
    for prfs, varkey in zip(varprfs, varkeys): 
        cgkeys = ["cg%i" %n for n in cldgrps_n]
        df_save = pd.DataFrame(dict(zip(cgkeys, prfs)))
        df_save['alt'] = gridalt
        df_save.to_csv(save_dir + "meanprf_%s_WP3.csv" % varkey, index=False)
    ##_________________________________________________________________________
        


    axset_wind[0].set_yticks(axset_wind[0].get_yticks())
    axset_wind[0].set_yticklabels(['' for t in axset_wind[1].get_yticks()])
    axset_wind[0].set_ylim(-100, 3300)
    axset_wind[0].set_xlabel(r"TKE$_h$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[0].set_xlim(0, 0.8)
    axset_wind[0].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    axset_wind[0].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)

    axset_wind[1].set_yticks(axset_wind[0].get_yticks())
    axset_wind[1].set_ylim(-100, 3300) 
    axset_wind[1].set_xlabel(r"$\bar{w'w'}$ (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[1].set_xlim(0, 0.6)
    axset_wind[1].set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    axset_wind[1].set_xticklabels(['0', '', '0.2', '', '0.4', '', '0.6'], fontsize=9)
    
    axset_wind[2].set_yticks(axset_wind[0].get_yticks())
    axset_wind[2].set_yticklabels(['' for t in axset_wind[1].get_yticks()])
    axset_wind[2].set_ylim(-100, 3300)
    axset_wind[2].set_xlabel(r"TKE (m$^2$ s$^{-2}$)", fontsize=12)
    axset_wind[2].set_xlim(0, 1.)
    axset_wind[2].set_xticks([0, 0.2, 0.4, 0.6, 0.8])
    axset_wind[2].set_xticklabels(['0', '0.2', '0.4', '0.6', '0.8'], fontsize=9)

    #axset_wind[3].set_yticks(axset_wind[0].get_yticks())
    #axset_wind[3].set_ylim(-100, 3300)
    ##axset_wind[3].set_xlabel(r"$\bar{w'w'w'}$", fontsize=12)
    #axset_wind[3].set_xlabel(r"$S_w$", fontsize=12)
    #axset_wind[3].set_xlim(-1.3, 3.1)
    #axset_wind[3].set_xticks([-1.2, -0.6, 0, 0.6, 1.2, 1.8, 2.4, 3])
    #axset_wind[3].set_xticklabels(['-1.2', '', '0', '', '1.2', '', '2.4', ''], fontsize=9)
    axset_wind[3].set_yticks(axset_wind[0].get_yticks())
    axset_wind[3].set_ylim(-100, 3300)
    axset_wind[3].set_xlabel("anisotropy \nratio", fontsize=12)
    axset_wind[3].set_xlim(0, 1.8)
    axset_wind[3].set_xticks([0, 0.5, 1, 1.5])
    axset_wind[3].set_xticklabels(['0', '0.5', '1.0', '1.5'], fontsize=9)



    axset_scalar[0].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[0].set_xlim(-25, 25) 
    axset_scalar[0].set_ylim(-100, 3300) 
    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)
    axset_scalar[0].set_xlabel(r"SHF (W m$^{-2}$)", fontsize=12)
    axset_scalar[0].set_xticks([-20, -10, 0, 10, 20])
    axset_scalar[0].set_xticklabels(['-20', '', '0', '', '20'], fontsize=9)

    axset_scalar[1].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[1].set_ylim(-100, 3300) 
    axset_scalar[1].set_xlabel(r"LHF (W m$^{-2}$)", fontsize=12)
    axset_scalar[1].set_xlim(-60, 490)
    axset_scalar[1].set_xticks([0, 100, 200, 300, 400])
    axset_scalar[1].set_xticklabels(['0', '100', '200', '300', '400'], fontsize=9)

    axset_scalar[2].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[2].set_xlim(-0.0007, 0.00165) 
    axset_scalar[2].set_ylim(-100, 3300) 
    axset_scalar[2].set_xlabel(r"F$_b$ ($m^2 s^{-3} 10^{-3}$)", fontsize=12)
    axset_scalar[2].set_xticks(np.array([-0.5, -0.25, 0, 0.25, 0.5, 0.75, 1., 1.25, 1.5])*1e-3)
    axset_scalar[2].set_xticklabels(['-0.5', '', '0', '', '0.5', '', '1.0', '', '1.5'], fontsize=9)

    
    axset_scalar[3].set_yticks(axset_wind[0].get_yticks())
    axset_scalar[3].set_ylim(-100, 3300) 
    axset_scalar[3].set_xlim(-250, 0) 
    axset_scalar[3].set_xlabel(r"$\delta D_{flux}$"+u'(\u2030)', fontsize=12)
    axset_scalar[3].set_xticks([-250, -200, -150, -100, -50, 0])
    axset_scalar[3].set_xticklabels(['', '-200', '', '-100', '', '0'], fontsize=9)

    #axset_scalar[3].set_yticks(axset_wind[0].get_yticks())
    #axset_scalar[3].set_ylim(-100, 3300) 
    #axset_scalar[3].set_xlim(-2.5, 2.5) 
    #axset_scalar[3].set_xlabel(r"$F_{b,s}/F_{b}$", fontsize=12)
    
    yticklabels = axset_wind[0].get_yticks().astype(int).astype(str)
    axset_wind[0].set_yticklabels(yticklabels, fontsize=9)
    axset_wind[0].set_ylabel('altitude (m)', fontsize=12)
    for ax in axset_wind[1:]:
        ax.set_yticklabels(['' for t in ax.get_yticks()])

    axset_scalar[0].set_ylabel('altitude (m)', fontsize=12)
    axset_scalar[0].set_yticklabels(yticklabels, fontsize=9)
    for ax in axset_scalar[1:]:
        ax.set_yticklabels(['' for t in ax.get_yticks()])

    
    fig_wind.savefig("./fig_wind+turb_profiles.png")
    fig_scalar.savefig("./fig_scalarflux_profiles.png")



if __name__=="__main__":
    fig_scatter()