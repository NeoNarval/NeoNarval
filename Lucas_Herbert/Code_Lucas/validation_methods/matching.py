#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python's modules imports
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
import validation_methods.compute_spikes as cpspikes
import validation_methods.read_ThAr_atlas as rdAtlas


"""
This script's objective is to use the data reduction's computing (offset computing, maxima and spikes computing and their gaussian fit, etc) to compare the output to the ThAr Atlas. The main goal is to see which wavelengths are matching and the precision of those matchings.
"""





""" 
This function will compare the two lists used as an input. It will return the proportion of values in the first list that find an equivalent with the given precision in the second list, and also a list with the matching values and their gap
"""
def data_matching(x,y,precision):
    
    matching_data = []
    
    for i in range(len(x)):
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        data = [x[i],best_matching_value,gap_min]            
        if (best_matching_value != None ):
            matching_data.insert(i,data)
    
    matching_rate = float(len(matching_data))/float(len(x))
    print("Matched, Total",len(matching_data),len(x),"Matching rate",matching_rate)
    return(matching_data,matching_rate)
    
    
"""
This functions returns the matching rate and the matchings gaps for one particulr order. The output data are the list of the matched wevelengths (couple of wavelengths in each list which are close enough to be be considered as the same one in reality) and of the corresponding gaps (distance between the two matched wavelengths) and finally the matching proportion among the order's different spikes wavelengths.
"""

def order_matching(n,precision):
    
    spikes_fits_data = cpspikes.fit_spikes_order(n)
    spikes_centers = []
    spikes_widths = []
    for i in range(len(spikes_fits_data)):
        spikes_centers.insert(i,spikes_fits_data[i][0])
        spikes_widths.insert(i,spikes_fits_data[i][1])
    plt.figure(2)
    plt.bar(spikes_centers,spikes_widths,color='red')
    atlas_lambdas = rdAtlas.read_ThAr_Atlas()
    data = data_matching(spikes_centers,atlas_lambdas,precision)
    gaps = []
    lambdas_gaps = []
    matching_data = data[0]
    
    for i in range(len(matching_data)):
        gaps.insert(i,matching_data[i][2])
        lambdas_gaps.insert(i,matching_data[i][0])
    plt.figure(2)
    plt.bar(lambdas_gaps,gaps, color='brown')
    plt.title("Error = f(Lambda)")
    plt.show() 
    average_error = float(np.sum(gaps)/len(matching_data))
    average_width = float(np.sum(spikes_widths)/len(spikes_widths))
    print("Average width",average_width)
    print("Average error",average_error)
    return(data[1],average_error,lambdas_gaps,gaps)
    


"""
This function calls order_matching for each order and print the results.
"""

def global_matching(precision):
    gaps = []
    lambdas_gaps = []
    for n in range(36):
        data = order_matching(n,precision)
        for i in range(len(data[2])):
            lambdas_gaps.append(data[2][i])
            gaps.append(data [3][i])
    return(lambdas_gaps,gaps)        






















