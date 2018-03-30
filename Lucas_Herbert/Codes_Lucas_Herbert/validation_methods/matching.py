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

# Validation methods imports :
import validation_methods.compute_spikes as cpspikes
import validation_methods.read_ThAr_atlas as rdAtlas
import validation_methods.tests_fit as tfit


"""
This script's objective is to use the data reduction's computing (offset computing, maxima and spikes computing and their gaussian fit, etc) to compare the output to the ThAr Atlas. The main goal is to see which wavelengths are matching and the precision of those matchings.
"""





""" 
This function will compare the two lists used as an input. It will return the proportion of values in the first list that find an equivalent with the given precision in the second list, and also a list with the matching values and their gap.
"""

def lambdas_matching(x,indices,y,precision):
    
    matching_data = []
    
    for i in range(len(x)):
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        data = [x[i],best_matching_value,gap_min,indices[i]]            
        if (best_matching_value != None ):
            matching_data.insert(i,data)
    
    if len(x) != 0 :
        matching_rate = float(len(matching_data))/float(len(x))
    else : 
        matching_rate = 0
    print("Matched, Total",len(matching_data),len(x),"Matching rate",matching_rate)
    return(matching_data,matching_rate)
    
    
"""
This functions returns the matching rate and the matchings gaps for one particulr order. The output data are the list of the matched wevelengths (couple of wavelengths in each list which are close enough to be be considered as the same one in reality) and of the corresponding gaps (distance between the two matched wavelengths) and finally the matching proportion among the order's different spikes wavelengths.
"""

def order_matching(n,precision):
    
    
    # First of all, we need to fit a gaussian on each spike, we will use the function "fit_the_spike"
    spikes_fits_data = cpspikes.fit_spikes_order(n)
    
    # Some of those fits are not good enough to be used : we have a factor, the "chi-square", which describes how good the fit is. It's a normalized indicator of the average error between the fitted value and the original value (g_fit(xi) - yi). Because of those bad fits, indicated by a higher chi-square, we are gonna filter the fits accordind to their chi-square. The critical chi-square is determined empirically by looking at the "wrong" fits and choosing above which chi-square the fits is wrong.
    critical_chi_square = 0.5
    # The same idea can be applied with the width : if the width of the fitted spike is more than a critical value, we can forget this spike because it's not a real spike (which would be thiner)
    critical_width = 0.1
    
    filtered_spikes_fits_data = []
    
    for i in range(len(spikes_fits_data)):
        
        chi_square = read_chi_square(spikes_fits_data[i][0][2])
        width = spikes_fits_data[i][0][1]
        
        if ( chi_square < critical_chi_square and  width < critical_width ) :
            filtered_spikes_fits_data.append(spikes_fits_data[i])
    
    filtered_spikes_fits_data.sort()

    spikes_centers_lambdas = []
    spikes_centers_indices = []
    spikes_widths = []
    spikes_fits_chi_square = []
    for i in range(len(filtered_spikes_fits_data)):
        spikes_centers_lambdas.insert(i,filtered_spikes_fits_data[i][0][0])
        spikes_centers_indices.insert(i,filtered_spikes_fits_data[i][1][0])
        spikes_widths.insert(i,filtered_spikes_fits_data[i][0][1])
        spikes_fits_chi_square.insert(i,read_chi_square(filtered_spikes_fits_data[i][0][2]))
    plt.figure(2)
    plt.plot(spikes_centers_lambdas,spikes_widths,'.',color='black')
    plt.plot(spikes_centers_lambdas,spikes_fits_chi_square,'.',color='purple')
    plt.figure(4)
    plt.plot(spikes_centers_indices,spikes_widths,'.',color='black')
    plt.plot(spikes_centers_indices,spikes_fits_chi_square,'.',color='purple')
    atlas_lambdas = rdAtlas.read_ThAr_Atlas()
    
    data = lambdas_matching(spikes_centers_lambdas,spikes_centers_indices,atlas_lambdas,precision)
    
    gaps = []
    matched_lambdas = []
    matched_lambdas_indices = []
    matching_data = data[0]
    
    for i in range(len(matching_data)):
        gaps.insert(i,matching_data[i][2])
        matched_lambdas.insert(i,matching_data[i][0])
        matched_lambdas_indices.insert(i,matching_data[i][3])
        
    plt.figure(2)
    plt.plot(matched_lambdas,gaps,'.',color='red')
    plt.title("Error = f(Lambda)")
    plt.xlabel("Wavelengths(Angstrom)")
    plt.show() 
    
    plt.figure(4)
    plt.plot(matched_lambdas_indices,gaps,'.',color='red')
    plt.title("Error = f(Indice)")
    plt.xlabel("Indices")
    plt.show() 
    average_error = float(np.sum(gaps)/len(matching_data))
    average_width = float(np.sum(spikes_widths)/len(spikes_widths))
    print("Order",n) 
    print("Average width",average_width)
    print("Average error",average_error)

    return(data[1],average_error,matched_lambdas,matched_lambdas_indices,gaps)
    
""" 
For each gaussian fit, there is a report which is a string countaining a lot of informations about the fit. The chi_square is a part of those informations which says if the fit is good or not, indicating the error between the fit and the fitted data. We need to get the chi_square from each report in order to discriminate the spikes fits. That's what this function does.

In order to get the chi_square from the report, here is a little function which read the report and finds the chi_square.
"""


def read_chi_square(report) :

    data = report[175:205]
    i = 0
    while ( data[i] != '=' ):
        i += 1
    i += 2
    chi_square_string = ''
    while ( data[i] != ' ' ):
        chi_square_string += data[i]
        i += 1
    chi_square = float(chi_square_string)        
    
    return(chi_square)

"""
This function calls order_matching for each order and print the results.
"""

def global_matching(precision):
    gaps = []
    lambdas_gaps = []
    indices = []
    for n in range(36):
        data = order_matching(n,precision)
        for i in range(len(data[2])):
            lambdas_gaps.append(data[2][i])
            gaps.append(data [4][i])
            indices.append(data[3][i])
    return(lambdas_gaps,gaps,indices)        




"""
A little function which records the matching data in a pickle.
"""

#path = '/home/stagiaire/Documents/Codes_Lucas_Herbert/matching_data_0.001_Angstrom.pkl'

def matching_pickle(path):
    
    data_file = open(path,'w')
    data = global_matching(1)
    
    matching_data = {}
    matching_data['Spikes_wavelengths'] = data[0]
    matching_data['Matching_gaps_between_drs_and_atlas'] = data[1]
    matching_data['Spikes_indices'] = data[2]
    matching_data_pickle = pickle.dump(matching_data,data_file)
    print("Recorded")
    data_file.close()
    
"""
A little function to classify the data of the given pkl. It takes a pkl file with wavelengths and the associated errors and regroup the wavelengths according to their errors.
"""

path = 'Données utiles/matching_data_1_Angstrom.pkl'

def matching_reader(path):
    
    file = open(path,'r')
    data = pickle.load(file)
    lambdas = data['Spikes_wavelengths']
    errors = data['Matching_gaps_between_drs_and_atlas']
    indices = data['Spikes_indices']
    plt.figure(5)
    plt.plot(lambdas,errors,'.',color='blue')
    
    # plt.figure(6)
    # plt.plot(indices,errors,'.',color='red')
    # plt.show()
    file.close()
    return(data)
    
    
'''
The functions above this message have been used as a test. It has not been effective so we let it here as a comment in order not to loose the code but we won't use it.

#"""
# A little function which takes the list of lambdas and errors from the order_matching and computes the affine function which fits those data and returns the a and b of the "ax+b".
# """
# 
# def affine_order_matching(n,precision):
#     
#     data = order_matching(n,precision)
#     
#     x = data[2]
#     y = data[3]
#     
#     return(tfit.fit_affine(x,y))
# 
# """
# A little function which takes the list of lambdas and errors from the global_matching and computes the affine function which fits those data and returns the a and b of the "ax+b".
# """
# 
# def affine_global_matching(precision) :
#     
#     data = global_matching(precision)
#     
#     x = data[0]
#     y = data[1]
#     
#     return(tfit.fit_affine(x,y))
#     
# """
# A little function which takes the list of lambdas and errors from the order_matching and computes the polynomial (degree 2) function which fits those data and returns the a, b  and c of the "ax²+bx+c".
# """
# 
# def pol2_order_matching(n,precision):
#     
#     data = order_matching(n,precision)
#     
#     x = data[2]
#     y = data[3]
#     
#     return(tfit.fit_pol2(x,y))
#     /home/stagiaire/Documents/Codes_Lucas_Herbert/matching_data_0.001_Angstrom.pkl'
# """
# Same function but global
# """
# 
# 
# def pol2_global_matching(precision):
#     
#     data = global_matching(precision)
#     
#     x = data[0]
#     y = data[1]
#     
#     return(tfit.fit_pol2(x,y))   
#     

'''


    
"""
Now that we are able to identify which wavelengths have a satisfying matching (under 1 Angstrom), we want to demonstrate that actually the majority is around 0.005A, which would correpond to the actual accurency of the calibration alogirithm. We want to print a graph where we can see the proportion of wavelengths corresponding to different groups of errors : [0,0.001],[0.001,003];[0.003,0.006],[0.006,0.01],etc.
That's the role of this function. 
"""
path = '/home/stagiaire/Documents/Codes_Lucas_Herbert/matching_data_1_Angstrom.pkl'

def errors_distribution(path):
    
        data = matching_reader(path)
        lambdas = data['Spikes_wavelengths']
        errors = data['Matching_gaps_between_drs_and_atlas']
        plt.figure(6)
        plt.hist(errors, bins=1000)
        plt.xlabel('Error in Angstrom')
        plt.ylabel('Proportion of the matched wavelengths')
        plt.show()
            
            
    
    
    
    
    