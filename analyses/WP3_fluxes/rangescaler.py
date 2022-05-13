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
    #if len(x) != len(xscale_bnds):
    #    print("x and xscale_bnds are not the same length")
    #    return
    
    
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
    xpiece_scaled = []
    for xp, bnds, bnds_scaled in zip(xpiece, x_bnds, xscaled_bnds):
        b1 = bnds[0]; b2 = bnds[1]
        bs1 = bnds_scaled[0]; bs2 = bnds_scaled[1]
        x_std = (xp - b1) / (b2 - b1)
        xpiece_scaled.append(x_std * (bs2 - bs1) + bs1)

        
    # Reshape:
    xscaled = np.array([])
    for xp in xpiece_scaled:
        xscaled = np.append(xscaled, xp)
        
    return xscaled
    

x = np.arange(0,100,1)
x_bnds = [(0,20), (20,40), (40,60), (60,100)]
xscaled_bnds = [(1,2), (2,5), (5,10), (10,11)]

xs = piecewise_linscale(x, x_bnds, xscaled_bnds)

plt.figure()
plt.scatter(x, xs)



















