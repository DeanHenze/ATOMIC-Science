# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:47:26 2022

@author: Dean
"""



import numpy as np
#from sklearn.preprocessing import minmax_scale
import matplotlib.pyplot as plt



def piecewise_linscale(x, x_bnds, xscaled_bnds):
    """
    Linear scaling of piecewise segments. Each segment is scaled to specified 
    bounds.
    
    Inputs
    ------
    x: iterable.
        Each element of x is a piecewise segment. Even if only one piece, 
        it should be put in an iterable.
    
    xscale_bnds: iterable of 2-tuples.
        The bounds to scale each segment by.
    """
    
    # Break x into piecewise segments:
    x_belowbnds = None
    i_belowbnds = np.where(x<x_bnds[0][0])[0]
    if len(i_belowbnds) != 0:
        x_belowbnds = x[i_belowbnds]
        
    x_abovebnds = None
    i_abovebnds = np.where(x>=x_bnds[-1][1])[0]
    if len(i_abovebnds) != 0:
        x_abovebnds = x[i_abovebnds]
        
    xpiece = []
    for bnds in x_bnds:
        b1 = bnds[0]
        b2 = bnds[1]
        xpiece.append(x[np.where((x>=b1) & (x<b2))])
        
        
    # Scaled segments:
    def linscale(xp, bnds, bnds_scaled):
        b1 = bnds[0]; b2 = bnds[1]
        bs1 = bnds_scaled[0]; bs2 = bnds_scaled[1]
        x_std = (xp - b1) / (b2 - b1)        
        return x_std * (bs2 - bs1) + bs1
        
    xscaled = np.array([])
    for xp, bnds, bnds_scaled in zip(xpiece, x_bnds, xscaled_bnds):
        xpiece_scaled = linscale(xp, bnds, bnds_scaled)
        xscaled = np.append(xscaled, xpiece_scaled)
    if x_belowbnds is not None:
        x_belowbnds_scaled = linscale(x_belowbnds, x_bnds[0], xscaled_bnds[0])
        xscaled = np.append(x_belowbnds_scaled, xscaled)
    if x_abovebnds is not None:
        x_abovebnds_scaled = linscale(x_abovebnds, x_bnds[-1], xscaled_bnds[-1])
        xscaled = np.append(xscaled, x_abovebnds_scaled)

        
    return xscaled
    

x = np.arange(0,100,1)
#x_bnds = [(0,20), (20,40), (40,60), (60,100)]
x_bnds = [(20,40), (40,60)]
xscaled_bnds = [(2,5), (5,10)]

xs = piecewise_linscale(x, x_bnds, xscaled_bnds)

plt.figure()
plt.scatter(x, xs)



















