#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Python's global modules imports :
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

# Validation methods imports :
import validation_methods.gaussian_fit as gfit
import validation_methods.compute_order as cporder

""" 
In this script, we will search for all the spikes of a normalized spectrum and then fit a gaussian to each of the interesting spikes in order to find its precise centre and width.
"""

file_path = '/home/stagiaire/Documents/Données utiles/th_calibered.fits' # current path computed
"""
Open the ".fits" file which we want to reduce in a table.
"""    
    
def read_fits(file_path):
    
    data_file = pyfits.open(file_path)
    
    # We search for our data in the .fits file : the wawelenghts and intensities of the thorium argon spectrum, it is in the first extension of the .fits file
    
    data = data_file[1].data
    
    # Now we have the right table, let's builds the two lists of values which interest us
    
    lambdas = data['wavelength_lane1']
    intensities = data['intensity_lane1']
    
    return(lambdas,intensities)
    
lambdas, intensities = read_fits(file_path)[0], read_fits(file_path)[1]


""" 
First function : we need to find the lists of wavelengths and intensities corresponding to each spike of the spectrum in order to go further. The following function takes as an input two lists : the wavelengths and their intensities. The output is the list of the data of each spike : each element of this list is a list of two lists : the lambdas and the intensities of the spike. 
"""

def find_spikes_data(lambdas,indices,intensities):
    
    spikes_data = []
    
    intensities_maxima = []
    lambdas_maxima = []
    maxima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] and intensities[i] >= 0  and intensities[i] >= 3) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
            maxima_indices.append(i)
    # Now we need to find the spike around each maximum
    
    for j in range(len(lambdas_maxima)):
        
        local_spike_data = []
        
        # left minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index-1] ):
            local_spike_data.append([lambdas[index],indices[index],intensities[index]])
            index -= 1
        local_spike_data.append([lambdas[index],indices[index],intensities[index]]) # don't forget the last pointindices_of_maxima = []
        
        # right minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index+1] ):
            local_spike_data.append([lambdas[index],indices[index],intensities[index]])
            index += 1
        local_spike_data.append([lambdas[index],indices[index],intensities[index]]) # don't forget the last point
        
        local_spike_data.sort() # We sort the list according to the lambdas order
        local_spike_lambdas = [ local_spike_data[i][0] for i in range(len(local_spike_data)) ]
        local_spike_indices = [ local_spike_data[i][1] for i in range(len(local_spike_data)) ]
        local_spike_intensities = [ local_spike_data[i][2] for i in range(len(local_spike_data)) ]
        spikes_data.append([local_spike_lambdas, local_spike_indices,local_spike_intensities])
    
    return(spikes_data)
    

"""
For each spike we have found, we need to localize more precisly its position to know its exact wavelengths. That's why we are gonna fit a gaussian on each spike thanks to another algorithm using lmfit and the least squares method. The centre of the gaussian will be the precise wavelength of the spike, and we can also have an access to informations like its width, all those data helping us to discriminate among the spikes. We will for example be able to filter the spikes which don't have a right fitting or which have a too large width, because they are not physically realistic enough.
"""


"""
The following function will compute the centers of the spikes of an given order and return data about the fits. 
"""

def fit_spikes_order(n) :
    
    spikes_fits_data = []
    # To verify the job has been done correctly, we can plot the different things we do.The original data are the wavelengths and intensities we have in the begining, plotted in black. 
    orders = cporder.search_orders(lambdas, intensities)
    order_lambdas = orders[n][0]
    order_indices = orders[n][1]
    order_intensities = cporder.compute_order(n)
    plt.figure(1)
    plt.plot(order_lambdas, order_intensities, color='black')
    plt.xlabel("Wavelengths(Angstrom)")
    plt.figure(3)
    plt.plot(order_indices, order_intensities, color='black')
    plt.xlabel("Indices")
    plt.show()
    # Then we find the differents spikes in those data and we can also plot them to show what we have considered as a spike or not. The only usefull lane in this function for now is the spikes_data computation.
    spikes_data = find_spikes_data(order_lambdas,order_indices, order_intensities)
    for i in range(len(spikes_data)):
        plt.figure(1)
        x = spikes_data[i][0]
        indices = spikes_data[i][1] 
        y = spikes_data[i][2]
        plt.plot(x,y, color='blue')
        plt.figure(3)
        plt.plot(indices,y, color='blue')
    
        # When we have the spikes data, we need to compute each spike with a gaussian fit.
        
        spike_fit_data = []
        plt.figure(1)
        spike_fit_data.insert(0,gfit.fit_the_spike(x,y))
        plt.figure(3)
        spike_fit_data.insert(1,gfit.fit_the_spike(indices,y))
        spikes_fits_data.insert(i,spike_fit_data)
        
    plt.show() 
      
    return(spikes_fits_data)












