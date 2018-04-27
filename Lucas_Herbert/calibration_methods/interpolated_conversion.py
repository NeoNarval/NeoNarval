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



"""
The module calibration_methods called interpolated_conversion uses the results obtained by the matching.py script from validation_methods. It takes the lists of wavelengths and indices of the matched spikes and use it to compute a polynom which fits the wavelengths given the indices for each order. The role of this script is to find the law giving the wavelength in Angstroms as a function of the indice (int). Thus we will be able to iterate the conversion as explained in the itered_calibration script.
"""


"""
A little function which opens the ".fits" file which we want to reduce in a table.
Input :
- file_path : name of the file.
Outputs :
- lambdas : list of wavelengths of the file
- intensities : list of the associated intensities ( = spectrum )
"""    
    
def read_fits(file_path):
    
    data_file = pyfits.open(file_path)
    
    # We search for our data in the .fits file : the wawelenghts and intensities of the thorium argon spectrum, it is in the first extension of the .fits file
    
    data = data_file[1].data
    
    # Now we have the right table, let's builds the two lists of values which interest us
    
    lambdas = data['wavelength_lane1']
    intensities = data['intensity_lane1']
        
    return(lambdas,intensities)
    
    
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


    file.close()
    return(data)

data = matching_reader(path)    
spikes_lambdas = data['Spikes_wavelengths']
spikes_indices = data['Spikes_indices']


"""
This function will reconstitute each order's indices and associated wavelengths, knowing the fact that there is exactly 4612 values (and thus also indices) per order.
Input : 
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- n int, the order's number- A pkl file recorded at the given path.
OutPut : 
- selected_order_indices : list containing the indices of the order n bewteen 0 and 4611.
- selected_order_wavelengths : list containing the wavelengths corresponding to those indices.
"""

def select_order(spikes_lambdas,spikes_indices,n):
    
    # Indices of the order which will be computed
    index_begining = 4612.0*n
    index_end = 4612.0*(n+1)
    
    selected_order_indices = []
    selected_order_wavelengths = []
    
    for i in range(len(spikes_indices)):
        
        if ( index_begining <= spikes_indices[i] <= index_end ):
            
            selected_order_indices.append(spikes_indices[i] )
            selected_order_wavelengths.append(spikes_lambdas[i])
    
    selected_order_indices.sort()
    selected_order_wavelengths.sort()
    
    plt.figure(18)
    plt.plot(selected_order_indices,selected_order_wavelengths,'yo')
    plt.show()
    
    return(selected_order_indices,selected_order_wavelengths)


    
"""
This function is the definition of a polynom.
Inputs :
- x : list of floats 
- coeffs : list of float representing the polynom's coefficients from the higher degree to the smaller.
Output :
- List of the float values f(x[i]) where f is the polynom's function and the x[i]  are x's values.
"""

def polynom(x,coeffs):
        
    return( [ np.sum([coeffs[len(coeffs)-i-1]*xi**(i) for i in range(len(coeffs))]) for xi in x ] )
    
    
    
"""
Now that we have found the indices and the wavelengths for each order, we need to find the mapping law between indices and wavelengths. We can do it with an interpolation between the matched values.
The following function will find the law : lambda = f(indice) by inteprolating between several wavelengths, for a given order.
Input :
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- n : int, number of the order.
Output :
- coeffs : list of float representing the interolation polynom's coefficients from the higher degree to the smaller.
"""

def polynomial_interpolation(spikes_lambdas,spikes_indices,n):
    
    
    # Order of the polynomial fit :
    p = 5
    
    # Selection of the rigth order to fit in indices and lambdas
    indices = select_order(spikes_lambdas,spikes_indices,n)[0]
    lambdas = select_order(spikes_lambdas,spikes_indices,n)[1]
    # Printing some details about the computation
    print("Order number :"+str(n))
    print("Number of used spikes : "+str(len(indices)))
    try :    
        coeffs = np.polyfit(indices,lambdas,p)
    except :
        coeffs = p*[0]
        print("Error in the interpolation : not enough spikes to inteprolate!")
    # Computation of the fit
    p = np.poly1d(coeffs)
    polyfitted = p(indices)
    plt.figure(18)
    plt.title("Interpolation results")
    plt.plot(indices,lambdas,'o',color='black')
    #plt.plot(indices,polyfitted,color='purple')
    plt.show()
    return(coeffs)



"""
The following function will convert a list from arturos to Angstroms.
Input :
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- arturos_list : list of int or float representing the list of the wavelengths in the arturo scale.
- order : int, number of the order to convert.
Output :
- order_lambdas_Angstroms : list of floats, conversion of arturos_list in Angstroms
"""



def arturo_convert_angstroms_order(spikes_lambdas,spikes_indices,order):
    
    print(" ")
    print(" _________________________ Converting  _________________________ ")
    print(" ")
    
    arturos_list = [i for i in range(0,166032)]
    lambdas_Angstroms = []
    
    order_lambdas_arturos = [] 
    order_lambdas_Angstroms = []
    
    # Selection of the arturos to convert :

    order_lambdas_arturos = arturos_list[4612*order:4612*(order+1)]

    # Convertion using the inteprolation :
    order_coeffs = polynomial_interpolation(spikes_lambdas,spikes_indices,order)
    print("Polynomial coefficients",order_coeffs)
    pol = np.poly1d(order_coeffs)
    order_lambdas_Angstroms = pol(order_lambdas_arturos)
    
    return(order_lambdas_Angstroms)

"""
This function converts all orders, using arturo_convert_angstroms_order.
Input :
- spikes_lambdas : list of the spikes wavelengths
- spikes_indices : list of the corresponding indices
- arturos_list : list of int or float representing the list of the wavelengths in the arturo scale.
Output :
- lambdas_Angstroms : list of floats, conversion of arturos_list in Angstroms for all orders.
"""
def arturo_convert_angstroms(spikes_lambdas,spikes_indices):
    
    lambdas_Angstroms = []
    
    for order in range(36):
        order_lambdas_Angstroms = arturo_convert_angstroms_order(spikes_lambdas,spikes_indices,order)
        for i in range(len(order_lambdas_Angstroms)):
            lambdas_Angstroms.insert(order+(i),order_lambdas_Angstroms[i])
            
            
    # We have to replace the two lasts order's wavelengths with the ThAr_calibered wevelengths because the interpolation is impossible due to the lack of spikes in the red wavelengths. That will avoid us many troubles later, it is not a good solution but it is the only one we have found yet.
    lambdas_ThAr_calibered = pickle.load(open("ThAr_calibered_original_lambdas.pkl",'r'))
    lambdas_Angstroms[0:2*4612] = lambdas_ThAr_calibered[0:2*4612]
    # Comparing to the original drs wavelengths with a little plot : 
    plt.figure(16)
    plt.title("Comparison : black = calibered lambdas, red = interpolated lambdas")
    plt.plot(lambdas_ThAr_calibered,color='black')
    plt.plot(lambdas_Angstroms,color='red')
    plt.show()
    

    return(lambdas_Angstroms)


"""
A little function which records the results of the previous functions in a pkl file.
Input :
- path : string, name of the file and path to the file
Output :
None.
"""


def convert_to_pickle(path):
    
    lambdas_converted = arturo_convert_angstroms(spikes_lambdas,spikes_indices,arturos_list)
    data_file = open(path,'w')
    data_pickle = pickle.dump(lambdas_converted,data_file)
    data_file.close()
    
    






















