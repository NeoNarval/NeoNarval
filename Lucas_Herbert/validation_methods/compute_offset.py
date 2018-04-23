#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python modules imports :
import numpy as np
import matplotlib.pyplot as plt


""" 
This script is used to deal with the offset which can make the computing of the spikes more difficult. We can, thanks to the following function, compute the offset and then deal with the spectrum minus the computed offset.
"""


"""
The following function searchs for the minima of a spectrum and then computes its offset by interpolating between those minima. There are some particular cases of blended spikes which we have to consider, that's why we use the average value of intensity to discriminate those cases.
Inputs :
- lambdas : list of wavelengths 
- intensities : list of the associated intensities (spectrum)
Output :
- spectrum : list of the normalized intensities after computation of the offset.
"""

def normalize_offset(lambdas, intensities):
    
    offset = [0 ]*len(intensities)
    i_minima = [] # intensities of the minima
    l_minima = [] # wavelengths of the minima
    minima_indices = [] # indices of the minima
    
    # The average is used to avoid considering the blended spikes as offset : when two spikes are too close, the minimum between those maxima is so high in intensity that if we cut it, we make a big mistake in the offset computation. We choose to avoid considering those case by just ignoring this kind of minimum
    average = (np.sum(intensities))/(len(intensities))
    
    for i in range ( 2 , len(lambdas)-1 ):
        
        # We can easily find the minima with the local mathematical definition
        if ( intensities[i-1] > intensities[i] and intensities[i] < intensities[i+1] and intensities[i] <= 3*average ) :  
            
            minima_indices.append(i)
            i_minima.append(intensities[i])
            l_minima.append(lambdas[i])
    
    # As we said, we use the average to avoid adding certains minima which are definitively not part of the offset


    # Now that we have our minima, we need to fill the list with interpolated values between the minima in order to reach the lenght of the intensities list...
    
    for k in range(0,len(l_minima)-1):

        delta_lambda  = l_minima[k+1] - l_minima[k]
        if (delta_lambda != 0):
            coeff_dir = (1/delta_lambda)*(i_minima[k+1]-i_minima[k])
        else : 
            coeff_dir = 0
        for i in range(minima_indices[k],minima_indices[k+1]):
            offset[i] = i_minima[k]+coeff_dir*(lambdas[i]-l_minima[k] )

    
    # Computing the flat part of the intensities at the beginning of each order :
    l=0
    while( intensities[l] == intensities[l+1] ):
        offset[l] = intensities[l]
        l += 1
    offset[l] = intensities[l]
    l +=1
    offset[l] = intensities[l] # Don't forget the last 2 points to avoid casualties
    # Then we have to fill between those points dealing with the flat part at the beginning and the first minimum : we are gonna do a simple interpolation
    if ((l_minima[0] - lambdas[l])!=0):
        local_coeff_begining = (i_minima[0] - intensities[l])/(l_minima[0] - lambdas[l])
    else : 
        local_coeff_begining = 0
    for i in range(l, minima_indices[0]):
        offset[i] = intensities[l]+local_coeff_begining*(lambdas[i]-lambdas[l])
    
    
    # Computing the flat part at the end of each order :
    l= len(intensities)-1
    while( intensities[l] == intensities[l-1] ):
        offset[l] = intensities[l]
        l -= 1    
    offset[l] = intensities[l] 
    l -= 1  
    offset[l] = intensities[l] # Don't forget the last 2 points to avoid casualties 
    # Then the same idea is applied : we inteprolate between the last minimum and this wavelengths which is the begining of the flat part at the end of the order.
    if ((lambdas[l] - l_minima[-1])!=0):    
        local_coeff_end = (intensities[l] - i_minima[-1])/(lambdas[l] - l_minima[-1])
    else : 
        local_coeff_end = 0
    for i in range(minima_indices[-1],l):
        offset[i] = i_minima[-1]+local_coeff_end*(lambdas[i]-l_minima[-1])     
    
    
    # Computation of the new spectrum without the offset :
    spectrum = [intensities[i] - offset[i] for i in range(len(intensities)) ]

    
    
    return(spectrum) 
    
    



