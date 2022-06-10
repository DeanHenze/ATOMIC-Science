# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 09:55:52 2022

@author: Dean
"""

import numpy as np
import xarray as xr



def pvals_cdflevs(pdf, cdfvals, dx=None, dy=None):
    """
    For a 1d or 2d probability distribution, find probability values corresponding 
    to requested cumulative probabilities. 
    
    For a requested cumulative prob of, for example, 0.75, find the 
    probability level p where if a contour line were drawn at p, the region 
    outside the contour would include 75% of the samples in the distribution.
    
    Inputs:
        pdf (1d or 2d array): values of a normalized probability distribution.
        
        cdfvals (1d list/tuple): complementary cumulative probabilities.
            If computed pvalues are going to be used with the matplotlib contour 
            fxn, ccdfvals should be monotonically decreasing.
        
        dx, dy (floats, default=none): The bin widths of the dimensions over which 
            the pdf was taken. For 1d distributions, only dx is used. For 2d, both 
            dx and dy are needed.
    """
    
    # Flatten pdf to 1d list if it is 2d:
    if pdf.ndim==2: 
        probflat = [item for sublist in pdf for item in sublist]
    
    # Sort the 1d list by increasing values and cast as np.array:
    psort = np.array(sorted(probflat))
        
    # Construct a cumulative distribution from this sorted list:
        # initialize cdf; set 1st element to 0:
    cdf = np.zeros(len(psort))
        # cdf[0] = psort[0] 
        # Compute cdf:
    for i in range(1,len(psort)): cdf[i] = cdf[i-1] + psort[i]*dx*dy
    
    # Need to remove elements from psort and cdf corresponding to any redundant 
    # leading 0's in cdf. These 0's mean that there are multiple points in the 
    # pdf with prob 0. They are not needed and the below algorithm requires cdf
    # to be monotonically increasing:
    i0s = np.where(cdf==0)[0]
    if len(i0s)!=0:
        psort = psort[i0s[-1]:] # Keep from index of last 0 onward.
        cdf = cdf[i0s[-1]:]

    # Find probability values corresponding to 1 minus each input cdf value. For 
    # example, if the contour level encompassing 90% of the data is desired, find 
    # the pvalue corresponding to cdf=0.1:    
    pvals = [] # Put all computed prob values here.
        # Put pdf into an xarray with cdf coords for easy interpolation to desired 
        # cdfvals:
    xa = xr.DataArray(psort, coords={'cdfvals':cdf}, dims=['cdfvals'])
    for v in cdfvals:
        pvals.append(xa.interp(cdfvals=(v)).values)
        
    return pvals



def pvals_ccdflevs(pdf, ccdfvals, dx=None, dy=None):
    """
    Similar to the 'pvals_cdflevs' fxn but for the complementary cumulative 
    probability.
    
    In this case, a requested complementary cumulative prob of, for example, 
    0.75, would locate the probability level p where if a contour line were 
    drawn at p, the region enclosed by the contour would include 75% of the 
    samples in the distribution.
    
    
    Inputs:
        pdf (1d or 2d array): values of a normalized probability distribution.
        
        ccdfvals (1d list/tuple): complementary cumulative probabilities.
            If computed pvalues are going to be used with the matplotlib contour 
            fxn, ccdfvals should be monotonically decreasing.
        
        dx, dy (floats, default=none): The bin widths of the dimensions over which 
            the pdf was taken. For 1d distributions, only dx is used. For 2d, both 
            dx and dy are needed.
    """
    
    cdfvals = 1 - np.array(ccdfvals)
    pvals = pvals_cdflevs(pdf, cdfvals, dx=dx, dy=dy)    
    return pvals