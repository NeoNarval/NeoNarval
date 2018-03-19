import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

""" 
In this script, we will search for all the spikes of a normalized spectrum and then fit a gaussian to each of the interesting spikes in order to find its precise centre and width.
"""




""" 
First function : we need to find the lists of wavelengths and intensities corresponding to each spike of the spectrum in order to go further.
"""

def find_spikes_data(lambdas,intensities):
    
    spikes_data = []
    
    intensities_maxima = []
    lambdas_maxima = []
    maxima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] and intensities[i] >= 0  and intensities[i] >= 0.5) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
            maxima_indices.append(i)
    # Now we need to find the slit around each maximum
    
    for j in range(len(lambdas_maxima)):
        
        local_spike_data = []
        
        # left minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index-1] ):
            local_spike_data.append([lambdas[index],intensities[index]])
            index -= 1
        local_spike_data.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        # right minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index+1] ):
            local_spike_data.append([lambdas[index],intensities[index]])
            index += 1
        local_spike_data.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        local_spike_data.sort() # We sort the list according to the lambdas order
        local_spike_lambdas = [ local_spike_data[i][0] for i in range(len(local_spike_data)) ]
        local_spike_intensities = [ local_spike_data[i][1] for i in range(len(local_spike_data)) ]
        spikes_data.append([local_spike_lambdas,local_spike_intensities])
    
    return(spikes_data)
    
    
for i in range(len(data)):
    plt.plot(data[i][0],data[i][1])
    plt.show()