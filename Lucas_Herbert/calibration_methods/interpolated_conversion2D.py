#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python's modules imports
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

# Other module's imports
import validation_methods.tests_fit as fit
import calibration_methods.interpolated_conversion as cnv 
import validation_methods.matching as mtch
import validation_methods.compute_order as cporder
import calibration_methods.itered_calibration as itrd

# Definition of global values
global order_len # length of an order
order_len = 4612
global order_total # number of orders
order_total = 36


"""
This module will use 2D polynomial fits instead of 1D. Basically does the same job as interpolated_conversion.
"""

"""
We can use our previous results to establish the correspondance between the indices scale and the wavelengths scale. Thus, we need to read those data, which is the role of the following function.
Input :
- path : string giving the path to the file containing the drs matching data.
Output : the pkl data. 
"""

path = 'matching_data_0.1_ultra_wide_selection_0.03Itreshold_Angstrom.pkl'

def matching_reader(path):
    
    file = open(path,'r')
    data = pickle.load(file)
    lambdas = data['Spikes_wavelengths']
    errors = data['Matching_gaps_between_drs_and_atlas']
    indices = data['Spikes_indices']
    
    
    # We will record those data in a particular pickle to be able to have an access later
    file.close()
    return(data)

data = matching_reader(path)    
spikes_lambdas = data['Spikes_wavelengths']
spikes_indices = data['Spikes_indices']
    
    
"""
This function will reconstitute each order's indices and associated wavelengths, knowing the fact that there is exactly  values (and thus also indices) per order.
Input : 
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- n int, the order's number- A pkl file recorded at the given path.
OutPut : 
- selected_order_indices : list containing the indices of the order n bewteen 0 and 4611.
- selected_order_wavelengths : list containing the wavelengths corresponding to those indices.
"""

def select_order2D(spikes_lambdas,spikes_indices,n):
    
    # Indices of the order which will be computed
    index_begining = order_len*n
    index_end = order_len*(n+1)
    
    selected_order_indices = []
    selected_order_wavelengths = []
    
    for i in range(len(spikes_indices)):
        
        if ( index_begining <= spikes_indices[i] <= index_end ):
            
            selected_order_indices.append(spikes_indices[i] )
            selected_order_wavelengths.append(spikes_lambdas[i])
    
    selected_order_wavelengths.sort()
    
    # The indices of the order can be between n*order_len and (n+1)*order_len but we now want it between 0 and order_len-1  so let's come back to this form :
    selected_order_indices2D = []
    
    for ind in selected_order_indices :
        index = ind % order_len
        selected_order_indices2D.append(index)
    
    selected_order_indices2D.sort()
    
    
    plt.figure(18)
    plt.plot(selected_order_indices2D,selected_order_wavelengths,'y+')
    plt.show()
    
    return(selected_order_indices2D,selected_order_wavelengths)

"""
Now that we have found the indices and the wavelengths for each order, we need to find the mapping law between indices and wavelengths. We can do it with an interpolation between the matched values.
The following function will find the law : lambda = f(indice) by interpolating between several wavelengths, for a given order.
Input :
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- n : int, number of the order.
Output :
- coeffs : list of float representing the interpolation polynom's coefficients from the higher degree to the smaller.
"""

def polynomial_interpolation_legendre2D(spikes_lambdas,spikes_indices,n):
    
    # Optical order
    m = 21 + n
    
    # Order of the polynomial fit :
    p = 10
    
    # Selection of the rigth order to fit in indices and lambdas
    indices = select_order2D(spikes_lambdas,spikes_indices,n)[0]
    lambdas = select_order2D(spikes_lambdas,spikes_indices,n)[1]
    # Printing some details about the computation
    print("Order number :"+str(m))
    print("Number of used spikes : "+str(len(indices)))
    absc = [i for i in range(order_len) ]

    try :    
    
        coeffs = np.polynomial.legendre.legfit(indices,lambdas,p)
        new_lambdas = np.polynomial.legendre.legval(absc,coeffs)
        
    except :
        new_lambdas = lambdas
        print("Error in the interpolation : not enough spikes to interpolate!")
    # Computation of the fit

    plt.figure(18)
    plt.title("Interpolation results")
    plt.plot(indices,lambdas,'+',color='black')
    plt.plot(absc,new_lambdas,color='red')
    plt.show()
    
    return(new_lambdas)
    
def plotallconversion():

    for i in range(34):
        polynomial_interpolation_legendre2D(spikes_lambdas,spikes_indices,i)
    
"""
The following function will compute a conversion using the function above, and then compute the matching or a given order. It will plot the results so that we can compare the effiency of the 2D fit versus the itered conversion. 
Inputs : 
- order : the number of the order to compute
Outputs : 
- None
"""

def order_conversion2D(order):
    
    clean_plt()
    order_lambdas = polynomial_interpolation_legendre2D(spikes_lambdas,spikes_indices,order)
    lambdas_file = open("temporary_file_for_2Dconversion",'w')
    pickle.dump(order_lambdas,lambdas_file)
    lambdas_file.close()
    
    matching_results = mtch.order_matching("temporary_file_for_2Dconversion",0.1,order,0.1)
    
    return(None)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
