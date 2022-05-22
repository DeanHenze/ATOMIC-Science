# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 09:17:59 2022

@author: Dean
"""

import numpy as np



def correlation(ds1, ds2, cmax_plus, cmax_minus): 
    """
    Time-lagged cross correlation calc. ds1 shifted with respect to ds2.
    Inputs: data sets 1 and 2, and the maximum shift (integer number of array elements),
    between the two to condiser for the calc.
    """
    
    #print('\t# Starting correlation fxn.')
    
    nc = cmax_plus + cmax_minus + 1
    
    ### Prelim stuff:
    a1 = np.array(ds1); a2 = np.array(ds2)
    l1 = len(a1); l2 = len(a2)
        # Cut short if the arrays do not have the same length:
    if l1 != l2:
        print('Data sets do not have the same number of elements')
        return [np.nan for i in range(nc)],[np.nan for i in range(nc)], \
               [np.nan for i in range(nc)]
        # cut short if the array lengths are less than the user-specified
        # maximum offset:
    if (l1/2) < np.max([cmax_plus, cmax_minus]):
        print('Arrays are too short for correlation calc at given maximum '
              'offset.')
        return [np.nan for i in range(nc)],[np.nan for i in range(nc)], \
               [np.nan for i in range(nc)]
    
    
    ### Calculate correlations for a1 shifted right wrt a2, from 0 to c_max:
    ###
        # corrR will hold the correlation calculations for right-shifts. The 
        # index of corrR corresponds to the amount by which a1 has been shifted
        # (include the 0 shift with the right shifts):
    corrR = np.zeros(cmax_plus+1) # the +1 is to include the case of 0 shift
        # initialize array N_R, which will hold the number of non-nan samples  
        # used to calc the correlation for each right-shift:
    N_R = np.zeros(cmax_plus+1) # the +1 is to include the case of 0 shift    
        # loop through correlation calc for each shift:
    for j in range(cmax_plus+1):
        # subsections of a1 and a2 which overlap after shift:
        s1 = a1[0:l1-j]; s2 = a2[j:l2]
        # only use subset of s1 and s2 where both of them have non-Nan and 
        # non-INF values:
        ss1 = s1[np.isfinite(s1*s2)]
        ss2 = s2[np.isfinite(s1*s2)]
        # define n as the length of s1 and s2 after all reductions in size:
        n = len(ss1)
        # calc means and standard deviations:
        mu1 = np.nanmean(ss1); std1 = np.nanstd(ss1)
        mu2 = np.nanmean(ss2); std2 = np.nanstd(ss2)
        # calculate the correlation, and assign it to the appropriate index
        # of corrR:
        corrR[j] = np.nansum((ss1-mu1)*(ss2-mu2))/(n*std1*std2)
        # assign n to N_R:
        N_R[j] = n
        
    ### Calculate correlations for a1 shifted left wrt a2, from -c_max to -1.
    ### This is accomplished by using the same algorithm as above but now 
    ### shifting a2 right with respect to a1:
    ###
        # corrL will hold the correlation calculations for left-shifts. The 
        # 1st element of corrL contains the calc for -c_max, and the last 
        # element contains the calc for -1:
    corrL = np.zeros(cmax_minus)
        # initialize array N_L, which will hold the number of non-nan samples  
        # used to calc the correlation for each left-shift:
    N_L = np.zeros(cmax_minus) # the +1 is to include the case of 0 shift 
        # Loop through correlation calc for each shift:
    for j in range(0, cmax_minus):
        # Subsection of a1 and a2. j runs from 0 to c_max-1 rather than 1 to 
        # c_max, so add 1 to j:
        s1 = a1[j+1:l1]; s2 = a2[0:l2-(j+1)] 
        # Only use subset of s1 and s2 where both of them have non-Nan and 
        # non_INF values:
        ss1 = s1[np.isfinite(s1*s2)]
        ss2 = s2[np.isfinite(s1*s2)]
        # Define n as the length of s1 and s2 after all reductions in size:
        n = len(ss1)
        # Calc means and standard deviations:
        mu1 = np.nanmean(ss1); std1 = np.nanstd(ss1)
        mu2 = np.nanmean(ss2); std2 = np.nanstd(ss2)
        # Calculate the correlation, and assign it to the appropriate index
        # of corrL (the extra -1 in the indexing brackets is to make python's
        # 0-indexing scheme work out):
        corrL[(cmax_minus-1)-j] = np.nansum((ss1-mu1)*(ss2-mu2))/(n*std1*std2) 
        # assign n to N_L:
        N_L[(cmax_minus-1)-j] = n
    
    # Assign to a var, c, all the shifts from -c_max to c_max performed above:
    c = np.arange(-cmax_minus, cmax_plus+1, 1)
    # Return the shifts, correlation calc for each shift, and number of non-
    # nan samples used in the corr calc for each shift:
    #print('\t# Completed correlation fxn.')
    return c, np.append(corrL,corrR), np.append(N_L,N_R)
