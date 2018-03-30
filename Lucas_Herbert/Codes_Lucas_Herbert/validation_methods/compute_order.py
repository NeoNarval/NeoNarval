#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python modules imports :
import pyfits
import numpy as np
import matplotlib.pyplot as plt
# Local methods imports :
import validation_methods.compute_offset as cpoffset

"""
All this script is written to open the file we need to compute,and then to compute each order.
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
This function uses two lists as an input. The first is the list of wavelengths and the second is the list of associated intensities. Since the orders are not isolated from each other, we need to build several lists, one for each order, representing the wavelengths and intensities of each order. This function does the job by returning the wavelengths and intensities of each order. 
"""


def search_orders(lambdas,intensities):
    
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
        else  :
            i+=1
                
        current_order = [current_order_lambdas, current_order_indices, current_order_intensities]
        orders.append(current_order)
    
    orders.sort()
    
    indice = 0
    for i in range(len(orders)):
        for j in range(len(orders[i][0])):
            orders[i][1].insert(j,indice)
            indice += 1
            
        
    return(orders)


    
"""
Once an order has been defined thanks to search_orders, we can compute its offset and return the precise normalized spectrum we want to study. That's the role of this function.
"""


    
def compute_order(n) :
    
    orders = search_orders(lambdas, intensities)
    order_lambdas = orders[n][0]
    order_indices = orders[n][1]
    order_intensities = orders[n][2]
    
    order_normalized_spectrum = cpoffset.normalize_offset(order_lambdas,order_intensities)
    
    plt.figure(1)
    plt.plot(order_lambdas,order_intensities,color='brown')
    plt.figure(3)
    plt.plot(order_indices,order_intensities,color='brown')
    plt.show()
    return(order_normalized_spectrum)
    

"""
We can also print an order after its computation
"""
def plot_order(n):
    
    orders = search_orders(lambdas,intensities)
    order_lambdas = orders[n][0]
    order_intensities = orders[n][2]
    order_indices = orders[n][1]
    computed_order_intensities = compute_order(n)
    plt.figure(1)
    plt.plot(order_lambdas,order_intensities, color='black')
    plt.plot(order_lambdas, computed_order_intensities, color='red')
    plt.xlabel("Wavelengths(Angstrom)")
    
    plt.figure(2)
    plt.plot(order_indices,order_intensities, color='black')
    plt.plot(order_indices, computed_order_intensities, color='red')
    plt.show()
    