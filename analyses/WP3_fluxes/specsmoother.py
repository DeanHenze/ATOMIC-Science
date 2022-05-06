# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 10:15:02 2022

@author: Dean

Smooths a spectral density function, applying increased smoothing at higher 
frequencies such that smoothing is proportional to fractional increases in 
frequency (e.g. smoothing increases logarithmically).

Original code in MATLAB courtesy of Chris Fairall. The code as been 
translated to Python here.
"""


import numpy as np


def spec_smoother(spec, frq_sample, Rdi=0.118):
    """
    Smooths spectra.
   
    Inputs:
    -------
    spec: 1D numpy array-like.
        The spectrum to be smoothed.
    frq_sample: float
        The sampling rate.
    Rdi: float.
        Larger number in Rdi results in less points in the smoothed spectrum.
        Some examples:
            Rdi=.2; 31 points output from 4097 input
            Rdi=.15; 41 points output from 4097 input
            Rdi=.118; 50 points output from 4097 input    
        
    Returns:
    --------
    spec_smooth: 
        The smoothed spectrum
    frq_smooth 
    The frequency of smoothed points
    """
    
    Mj = 0 #Mj = 1, matlab starts at 1 indexing
    Di = 0
    spec_len = len(spec)
    ii = 0
    jj=0 #jj=1, matlab starts at 1 indexing
    Dfx = frq_sample/(spec_len*2)        # frq_sample==sampling rate
    
    
    spec_smooth = [] # Smoothed spectra appended here.
    frq_smooth = [] # Frequencies for each smoothed spectra point appended here.
    
    while ii < spec_len:
        if ii==(spec_len-1): break
       
        # Frequency window:
        Di = round(Rdi*(Mj+Di));
        frq_smooth.append((Mj+Di/2)*Dfx) # Corresponding frequency of the window.
        
        
        # Get average spectal value over frequency window:
        sum_specwindow = 0
        for ii in range(Mj, Mj+Di+1):
            if ii == (spec_len-1): break
            sum_specwindow += spec[ii+1]
        spec_smooth.append(sum_specwindow/(Di+1))

        
        Mj = Mj+Di+1
        jj = jj+1
    
    
    spec_smooth = np.array(spec_smooth[0:len(frq_smooth)-1])
    frq_smooth = np.array(frq_smooth[0:len(frq_smooth)-1])
    
    return frq_smooth, spec_smooth
    