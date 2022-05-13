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
    if len(x_bnds) != len(xscaled_bnds):
        print("Length of x_bnds and xscaled_bnds must match.")
        return None
    if len(x_bnds)<2:
        print("Must specify at least 2 bounds.")
        return None
    
    
    # Break x into piecewise segments:
    xpiece_below = None # Any values below the lowest bound go here.
    i_belowbnds = np.where(x<x_bnds[0])[0]
    if len(i_belowbnds) != 0:
        xpiece_below = x[i_belowbnds]
        
    x_pieceabove = None # Any values above the highest bound go here.
    i_abovebnds = np.where(x>=x_bnds[-1])[0]
    if len(i_abovebnds) != 0:
        xpiece_above = x[i_abovebnds]
        
    xpiece = [] # All segments within bounds here.
    for i in range(len(x_bnds)-1):
        b1 = x_bnds[i]
        b2 = x_bnds[i+1]
        xpiece.append(x[np.where((x>=b1) & (x<b2))])
        
        
    # Scale segments:
    def linscale(xp, b1, b2, bs1, bs2):
        # xp = values to scale.
        # b1, b2 = unscaled bounds; bs1, bs2 = scaled bounds
        x_std = (xp - b1) / (b2 - b1)        
        return x_std * (bs2 - bs1) + bs1
        
    xscaled = np.array([])
        
    for i in range(len(x_bnds)-1): # Segments within bounds.
        xpiece_scaled = linscale(
            xpiece[i], 
            x_bnds[i], x_bnds[i+1], xscaled_bnds[i], xscaled_bnds[i+1]
            )
        xscaled = np.append(xscaled, xpiece_scaled)
    
    if xpiece_below is not None:
        xpiece_below_scaled = linscale(
            xpiece_below, 
            x_bnds[0], x_bnds[1], xscaled_bnds[0], xscaled_bnds[1]
            )
        xscaled = np.append(xpiece_below_scaled, xscaled)
    
    if xpiece_above is not None:
        xpiece_above_scaled = linscale(
            xpiece_above, 
            x_bnds[-2], x_bnds[-1], xscaled_bnds[-2], xscaled_bnds[-1]
            )
        xscaled = np.append(xscaled, xpiece_above_scaled)

        
    return xscaled
    


def quicktest():
    x = np.arange(0,100,1)
    #x_bnds = (0,20,40,60,100)
    #xscaled_bnds = (1,2,5,10,11)
    x_bnds = (20,40,60)
    xscaled_bnds = (2,5,10)
    
    xs = piecewise_linscale(x, x_bnds, xscaled_bnds)
    
    plt.figure()
    plt.scatter(x, xs)
    plt.xlabel("unscaled x", fontsize=14)
    plt.ylabel("scaled x", fontsize=14)


















