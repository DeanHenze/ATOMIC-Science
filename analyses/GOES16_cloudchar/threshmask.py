# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 09:03:48 2022

@author: Dean
"""


import numpy as np


def threshmask(data, low_thresh=None, high_thresh=None):
    """
    Returns boolean array-like of values higher than low_thresh and lower than 
    high_thresh.
    
    data: array-like.
    
    low_thresh, high_thresh: scalars.
        Values relevant to 'data'.
    """
    if low_thresh is None: low_thresh = -np.inf
    if high_thresh is None: high_thresh = np.inf
    
    mask1 = data > low_thresh
    mask2 = data < high_thresh
    return mask1 & mask2