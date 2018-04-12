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
This function will compare the two lists used as an input. It will return the proportion of values in the first list that find an equivalent with the given precision in the second list, and also a list with the matching values and their gaps.
Inputs :
- x : list of the values of the first list we want to compare (drs wavelengths list actually).
- indices : list of the associated indices (we want to know which x matches with which y and what is the index of each of those x because we will use them later).
- y : second list of values (atlas wavelengths list actually).
- precision : float representing the maximum possible difference between one value in x and one value in y to consider them as matching values.
Outputs :
- matching_data : list of lists containing the data of each matched values : the two matching values in x and y and their difference (called "gap") and the associated index for the x value.
- matching_rate : float representing the proportion of the x's values which find a matching value in y.
"""

def lambdas_matching(x,indices,y,precision):
    
    matching_data = []  # initiating the list we will use to register the matching wavelengths and their differences
    
    for i in range(len(x)):
        
        
        # initiating the minimum gap and best matching value research
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        # Now that we have the minimum gap and the matching value, we can register this data
        data = [x[i],best_matching_value,gap_min,indices[i]]            
        if (best_matching_value != None ):
            matching_data.insert(i,data)
    
    if len(x) != 0 :
        matching_rate = float(len(matching_data))/float(len(x))
    else : 
        matching_rate = 0
    print("Matched, Total : "+str(len(matching_data))+"/"+str(len(x))+" => Matching rate = "+str(matching_rate))
    return(matching_data,matching_rate)
    
    
"""
This functions returns the matching rate and the matching data for one particular order, between this order's computing spikes wavelengths and the ThAr Atlas wavelengths, using lambdas_matching.
Before computing the matchings, it filters the spikes wavelengths according to the accuracy of their fit and the thickness of their width in order to keep only the right values (those which are accurate enough to be compared without doing a big mistake). This function will use the read_chi_square function which is written and described below to filter the chi_square values.
Inputs : 
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use) 
- n : int, number of the order
- precision : float representing the maximum possible difference between one value in x and one value in y to consider them as matching values.
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
Outputs :
- data[1] : this is the matching_rate from the lambdas_matching of the order's spikes wavelengths and the Atlas wavelengths.
- average_error : float which is the average gap between the matched values.
- matched_lambdas : list of the wavelengths of the spikes which have found a corresponding Altas wavelength closer than the input precision. 
EDIT : we will try to replace the list of wavelengths matching by the Atlas wavelengths which have matched, but ONLY if precision <= 0.01, in order to have a better impatch on conversion with the iterated process of the calibrarion_methods module.
- matched_lambdas_indices : list of the indices of the matched_lambdas (those are actually floats because those are the centers of the indices gaussian fits), which will be used in the interpolated_coversion algorithm.
- gaps : list of the diffrences between the matched wavelengths
"""

def order_matching(lambdas_path,max_detection,n,precision):
    print(" ")
    print(" _________________________ Matching  _________________________ ")
    print(" ")
    print("Order : "+str(n)) 
    # First of all, we need to fit a gaussian on each spike, we will use the function "fit_the_spike"
    spikes_fits_data = cpspikes.fit_spikes_order(lambdas_path,max_detection,n)
    
    # Remember, the spikes_fits_data is a list of spike_fit_data, thats is to say a list of lists giving informations about a spike's fit. In each list in the spikes_fits_data list, we find the data concerning the two fits (one in wavelengths in Angstroms and the other in indices).
    

    
    # Some of those fits are not good enough to be used : we have a factor, the "chi-square", which describes how good the fit is. It's a normalized indicator of the average error between the fitted value and the original value (g_fit(xi) - yi). Because of those bad fits, indicated by a higher chi-square, we are gonna filter the fits accordind to their chi-square. The critical chi-square is determined empirically by looking at the "wrong" fits and choosing above which chi-square the fits is wrong.
    critical_chi_square = 0.5
    # The same idea can be applied with the width : if the width of the fitted spike is more than a critical value, we can forget this spike because it's not a real spike (which would be thiner)
    critical_width = 0.1
    
    
    
    
    filtered_spikes_fits_data = []  # initiating the usefull data list
    
    
    # this loop will select the spikes with good chi_squares and widths
    for i in range(len(spikes_fits_data)):
        report = spikes_fits_data[i][0][2]
        if (report != "Computation failed for this spike : default data = naive fit"):
            chi_square = read_chi_square(report)
        else :
            chi_square = 100 
        width = spikes_fits_data[i][0][1]
        
        if ( chi_square < critical_chi_square and  width < critical_width ) :
            filtered_spikes_fits_data.append(spikes_fits_data[i])
    
    filtered_spikes_fits_data.sort()
    
    # initiating some more lists which will countain the specific data for each spike's fit : its wavelength_center, its indice_center, its width and its chi_square. 
    spikes_centers_lambdas = []
    spikes_centers_indices = []
    spikes_widths = []
    spikes_fits_chi_square = []
    
    for i in range(len(filtered_spikes_fits_data)):
        # let's register all the data of each spike's fit in the filtered list, all in wavelenghts (lambdas)
        spikes_centers_lambdas.insert(i,filtered_spikes_fits_data[i][0][0]) 
        spikes_widths.insert(i,filtered_spikes_fits_data[i][0][1])
        spikes_fits_chi_square.insert(i,read_chi_square(filtered_spikes_fits_data[i][0][2]))
        # in indices, only the center is interesting for us
        spikes_centers_indices.insert(i,filtered_spikes_fits_data[i][1][0])
        
    # We can plot all of those data    
    # plt.figure(2)
    # plt.plot(spikes_centers_lambdas,spikes_widths,'.',color='black')
    # plt.plot(spikes_centers_lambdas,spikes_fits_chi_square,'.',color='purple')
    # plt.figure(4)
    # plt.plot(spikes_centers_indices,spikes_widths,'.',color='black')
    # plt.plot(spikes_centers_indices,spikes_fits_chi_square,'.',color='purple')
    
    
    # Now that we have all what we need to compare to the atlas, let's get the atlas spikes wavelengths.
    atlas_lambdas = rdAtlas.read_ThAr_Atlas()
    
    # Let's compare our data with the atlas data.
    data = lambdas_matching(spikes_centers_lambdas,spikes_centers_indices,atlas_lambdas,precision)
    
    
    # We need to initialize to last lists to countain the comparison's results.
    gaps = []
    matched_lambdas_drs = []
    matched_lambdas_atlas = []
    matched_lambdas_indices = []
    matching_data = data[0]
    
    for i in range(len(matching_data)):
        gaps.insert(i,matching_data[i][2])
        matched_lambdas_drs.insert(i,matching_data[i][0])
        matched_lambdas_atlas.insert(i,matching_data[i][1])
        matched_lambdas_indices.insert(i,matching_data[i][3])
        
        
    # We can plot the comparison's result for the in different scales (wavelengths and indices)
    plt.figure(2)
    plt.plot(matched_lambdas_atlas,gaps,'yo')
    plt.plot(matched_lambdas_drs,gaps,'.',color='red')
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
    #print("Average width",average_width)
    print("Average error : "+str(average_error))
    
    
    # Modification of the returned wavelengths : if the precision is high enough, that can be understood as if there is no doubt about the matching, we can replace the drs wavelengths by the atlas wavlengths in the output of this function because it will be used to compute a new conversion which can benefit from this idea. It is only a try for the moment.
    
    matched_lambdas = []
    
    for i in range(len(matched_lambdas_drs)) :
        if ( gaps[i] <= 0.01 ) :
            matched_lambdas.insert(i,matched_lambdas_atlas[i])
        else :
            matched_lambdas.insert(i,matched_lambdas_drs[i])
    
    
    return(data[1],average_error,matched_lambdas,matched_lambdas_indices,gaps)
    
    
    
""" 
For each gaussian fit, there is a report which is a string countaining a lot of informations about the fit. The chi_square is a part of those informations which says if the fit is good or not, indicating the error between the fit and the fitted data. We need to get the chi_square from each report in order to discriminate the spikes fits. That's what this function does.
In order to get the chi_square from the report, here is a little function which read the report and finds the chi_square.
Input : 
- report : string returned by the fitting algorithm giving all the informations about the fit, thus also giving the chi_square.
Output :
- chi_square : float giving information about the accuracy of the fit. The smaller it is, the more accurate is our fit.
"""


def read_chi_square(report) :
    
    data = report[172:205]
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
Inputs : 
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use) 
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : float representing the maximum possible difference between one value in x and one value in y to consider them as matching values.
Outputs :
- lambdas_gaps : list of the matching wavelengths
- gaps : list of the associated gaps
- indices : list of the associated indices
"""

def global_matching(lambdas_path,max_detection,precision):
    gaps = []
    lambdas_gaps = []
    indices = []
    for n in range(36):
        data = order_matching(lambdas_path,max_detection,n,precision)
        for i in range(len(data[2])):
            lambdas_gaps.append(data[2][i])
            gaps.append(data [4][i])
            indices.append(data[3][i])
    return(lambdas_gaps,gaps,indices)        


"""
A little function which records the matching data in a pickle, for a given order.
Inputs :
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use)
- order : int, order which will be computed 
- record_path : string indicating the file name.
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : float representing the maximum possible difference between one value in x and one value in y to consider them as matching values.
Output :
None.
"""

path = 'matching_data_0.1A_selection_0.01Itreshold.pkl'

def order_matching_pickle(lambdas_path,order,record_path,max_detection,precision):
    
    data_file = open(record_path,'w') # opening the file
    data = order_matching(lambdas_path,max_detection,order,precision)
    
    matching_data = {} # we will use a dictionary to register our data
    matching_data['Spikes_wavelengths'] = data[2]  # precising the name of the columns
    matching_data['Matching_gaps_between_drs_and_atlas'] = data[4]
    matching_data['Spikes_indices'] = data[3]
    matching_data_pickle = pickle.dump(matching_data,data_file) # storing of the data
    #print("Matching pickle recorded")
    data_file.close()
 

"""
A little function which records the matching data in a pickle.
Inputs :
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use) 
- record_path : string indicating the file name.
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : float representing the maximum possible difference between one value in x and one value in y to consider them as matching values.
Output :
None.
"""


def matching_pickle(lambdas_path,record_path,max_detection,precision):
    
    data_file = open(record_path,'w') # opening the file
    data = global_matching(lambdas_path,max_detection,precision)
    
    matching_data = {} # we will use a dictionary to register our data
    matching_data['Spikes_wavelengths'] = data[0]  # precising the name of the columns
    matching_data['Matching_gaps_between_drs_and_atlas'] = data[1]
    matching_data['Spikes_indices'] = data[2]
    matching_data_pickle = pickle.dump(matching_data,data_file) # storing of the data
    #print("Macthing pickle recorded")
    data_file.close()
    
"""
A little function which reads the matching data from a pkl file.
Input : 
- path : string indicating the name of the file
Output : 
- data : data contained by the file
"""


def matching_reader(path):
    
    file = open(path,'r')  # opening the file
    data = pickle.load(file)  
    lambdas = data['Spikes_wavelengths']
    errors = data['Matching_gaps_between_drs_and_atlas']
    indices = data['Spikes_indices']
    file.close()
    return(data)
    

    
"""
Now that we are able to identify which wavelengths have a satisfying matching (under "precision" Angstrom), we want to demonstrate that actually the majority is around 0.005A, which would correpond to the actual accurency of the calibration algorithm. We want to plot a graph where we can see the proportion of wavelengths corresponding to different groups of errors : [0,0.001],[0.001,003];[0.003,0.006],[0.006,0.01],etc.
That's the role of this function. 
Input :
- path : name of the pkl file where the matching data are stored.
Output : 
None.
"""


def errors_distribution(path):
    
        data = matching_reader(path)
        lambdas = data['Spikes_wavelengths']
        errors = data['Matching_gaps_between_drs_and_atlas']
        
        standard_deviation = np.std(errors)
        plt.close(25)
        print("Standard deviation : "+str(standard_deviation))
        plt.figure(25)
        plt.hist(errors, bins=100,color='red')
        plt.xlabel('Error in Angstrom')
        plt.ylabel('Proportion of the matched wavelengths')
        plt.show()
            
            
    
    
    
    
    