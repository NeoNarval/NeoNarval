#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python modules imports :
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pickle
# Local methods imports :
import validation_methods.compute_offset as cpoffset

"""
All this script is written to open the file we need to compute,and then to compute each order.
"""
file_path = 'validation_methods/Validation_files/th_calibered.fits' # current path computed




    
def read_fits(file_path):
    
    """
    A little function which opens the ".fits" file which we want to reduce in a table.
    Input :
    - file_path : name of the file.
    Outputs :
    - lambdas : list of wavelengths of the file
    - intensities : list of the associated intensities ( = spectrum )
    """    
    data_file = pyfits.open(file_path)
    
    # We search for our data in the .fits file : the wawelenghts and intensities of the thorium argon spectrum, it is in the first extension of the .fits file
    
    data = data_file[1].data
    
    # Now we have the right table, let's builds the two lists of values which interest us
    
    lambdas = data['wavelength_lane1']
    intensities = data['intensity_lane1']
        
    return(lambdas,intensities)
    
    

def fits_to_pkl():
        
    """
    The following function will record the lambdas and intensities from the original ThAr_calibered.fits) in a pickle so that we can use it in an easier way than beofre (with the fits file).
    Inputs :
    None
    Output : 
    None
    """    

    lambdas = read_fits(file_path)[0]
    intensities = read_fits(file_path)[1]
    
    lambdas_file = open("validation_methods/Validation_files/ThAr_calibered_lambdas.pkl",'w')
    intensities_file = open("validation_methods/Validation_files/ThAr_calibered_original_intensitites.pkl",'w')
    
    lambdas_pkl = pickle.dump(lambdas,lambdas_file)
    intensities_pkl = pickle.dump(intensities,intensities_file)
    
    lambdas_file.close()
    intensities_file.close()
    
    print("ThAr_calibered_original_data recorded!")
    return(None)




def search_orders(lambdas_pkl_path,intensities_pkl_path):
    
    """
    This function uses two path to lists as an input. The first path is leading to the list of wavelengths and the second to the list of associated intensities. Since the orders are not isolated from each other, we need to build several lists, one for each order, representing the wavelengths and intensities of each order. This function does the job by returning the wavelengths and intensities of each order. 
    Inputs: 
    - lambdas_pkl_path : path of the file containing the list of wavelengths
    - intensities : path of the file containing the list of the associated intensities ( = spectrum )
    Outputs :
    - orders : list of lists containing each order wavelengths and intensities
    """
    # ThAr_calibered_lambdas = read_fits(file_path)[0]
    lambdas = pickle.load(open(lambdas_pkl_path,'r'))
    # We have to complete the two lasts order of the wavelengths list with the ThAr list because when we convert a new wavelengths list we don't have values in those two order due to the lack of spikes in the red wavelengths. You can refer to the comments of the module called "calibration_methods" and "matching.py" to understand.
    
    intensities = pickle.load(open(intensities_pkl_path,'r'))
    
    orders = []
    i=0
    
    while ( i <= len(lambdas)-2 ) :
        current_order_lambdas = []
        current_order_indices = []
        current_order_intensities = []
        while ( i<= len(lambdas)-2 and lambdas[i] - lambdas[i+1] <= 50 ):
            current_order_lambdas.append(lambdas[i])
            current_order_intensities.append(intensities[i])
            i+=1
        else :
            i+=1
            
        # We add the last values that the "whil" loop has forgotten which are at the end of each order but sill part of the order
        current_order_lambdas.append(lambdas[i-1])
        current_order_intensities.append(intensities[i-1])    
        
        current_order = [current_order_lambdas, current_order_indices, current_order_intensities]
        orders.append(current_order)
    
    orders.sort()
    
    indice = 0
    for i in range(len(orders)):
        for j in range(len(orders[i][0])):
            orders[i][1].insert(j,indice)
            indice += 1
            
        
    return(orders)




def compute_order(n) :
    
    
    """
    Once an order has been defined thanks to search_orders, we can compute its offset and return the precise normalized spectrum we want to study. That's the role of this function. You can notice that we will always the "ThAr_calibered_lambdas.pkl" file because we don't care about which wavelengths we choose.
    Input :
    - n : int, number of the order to compute
    Output :
    - order_normalized_spectrum : list of the normalized intensities of the computed order.
    """
    orders = search_orders("validation_methods/Validation_files//ThAr_calibered_lambdas.pkl","validation_methods/Validation_files/ThAr_calibered_original_intensitites.pkl")
    order_lambdas = orders[n][0]
    order_indices = orders[n][1]
    order_intensities = orders[n][2]
    
    order_normalized_spectrum = cpoffset.normalize_offset(order_lambdas,order_intensities)
    
    # plt.figure(1)
    # plt.plot(order_lambdas,order_intensities,color='brown')
    # plt.figure(3)
    # plt.plot(order_indices,order_intensities,color='brown')
    # plt.show()
    return(order_normalized_spectrum)
    
    

def plot_order(n):
    
    
    """
    We can also plot an order after its computation.
    Input :
    - n : int, number of the order to plot.
    Output :
    None.
    """
    orders = search_orders("validation_methods/Validation_files/ThAr_calibered_lambdas.pkl","validation_methods/Validation_files/ThAr_calibered_original_intensitites.pkl")
    order_lambdas = orders[n][0]
    order_intensities = orders[n][2]
    order_indices = orders[n][1]
    computed_order_intensities = compute_order(n)
    # plt.figure(1)
    #plt.plot(order_lambdas,order_intensities, color='black')
    # plt.plot(order_lambdas, computed_order_intensities, color='red')
    # plt.xlabel("Wavelengths(Angstrom)")
    
    plt.figure(4)
    plt.plot(order_indices,order_intensities, color='black')
    plt.figure(5)
    plt.plot(order_indices, computed_order_intensities, color='purple')
    
    
    plt.show()

    