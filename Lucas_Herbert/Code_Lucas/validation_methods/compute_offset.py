import numpy as np
import matplotlib.pyplot as plt
""" 
This script is used to deal with the offset which can make the computing of the spikes more difficult. We can, thanks to the following function, compute the offset and then deal with the spectrum minus the computed offset.
The following function searchs for the minima of a spectrum and then computes its offset by interpolating between those minima. There are some particular cases of blended spikes which we have to consider, that's why we use the average value of intensity to discriminate thsoe cases.
"""

def normalize_offset(lambdas, intensities):
    
    offset = [0*len(intensities)]
    i_minima = []
    l_minima = []
    minima_indices = []
    
    # The average is used to avoid considering the blended spikes as offset : when two spikes are too close, the minimum between those maxima is so high in intensity that if we cut it, we make a big mistake in the offset computation. We choose to avoid considering those case by just ignoring this kind of minimum
    average = (np.sum(intensities))/(len(intensities))
    
    for i in range ( 2 , len(lambdas)-1 ):
        
        if ( intensities[i-1] >= intensities[i] and intensities[i] <= intensities[i+1] and intensities[i] <= 3*average ) :  
            
            minima_indices.append(i)
            i_minima.append(intensities[i])
            l_minima.append(lambdas[i])
    
    # As we said, we use the average to avoid adding certains minima which are definitively not part of the offset


    # Now that we have our minima, we need to fill the list with interpolated values between the minima in order to reach the lenght of the intensities list...
    
    for k in range(0,len(l_minima)-1):

        delta_lambda  = l_minima[k+1] - l_minima[k]
        coeff_dir = (1/delta_lambda)*(i_minima[k+1]-i_minima[k])
        for i in range(minima_indices[k],minima_indices[k+1]):
            offset.insert(i,  i_minima[k]+coeff_dir*(lambdas[i]-l_minima[k] ) )
        
    for i in range(0,minima_indices[0]-1):
        offset.insert(i,intensities[i])
        
    for i in range(minima_indices[-1],len(intensities)):
         offset.insert(i,intensities[i]) 
        
    spectrum = [intensities[i] - offset[i] for i in range(len(intensities)) ]
    
    return(spectrum) 
    
    



