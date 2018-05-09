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

import calibration_methods.interpolated_conversion as cnv 
import validation_methods.matching as mtch
import validation_methods.compute_order as cporder

# Definition of global values
global order_len # length of an order
order_len = 4612
global order_total # number of orders
order_total = 36

"""
This python module's role is to use the previous functions in validation_methods and calibration_methods in order to improve the calibration. It starts by using an already converted wavelengths list (converted by the cross correlation method used by Arthur Fezard in Thar3.py) to find which were the right spikes in this list thanks to the comparison with the ThAr Atlas. When those matching wavelengths have been found, it interpolates between those values to find each order's conversion law. Thanks to this law, it is able to convert the original arturo's scaled wavelengths list to an Angstrom's scaled list. This one is used to find the matching spikes, which are used to interpolate the conversion law, which is used to find another wavelengths list, etc. By iterating the process, we hope we can find the better possible conversion. The algorithm will stop the iteration when the errors between matched wavelengths and converted wavelengths are satisfying. We can use the standard deviation of the errors between matched and converted wavelengths as an indicator of the precision of our conversion.
"""



    
"""
The following function will use the validation_methods module's scripts to generate the list of matching wavelengths we need to achieve the interpolation, for a given order. It will also plot some usefull informations. Then it will use this matching wavelengths list to compute a new conversion and record the result (which is a new wavelengths list from which we can iterate the process).
Inputs :
- order : int, number of the order
/!\ If order = 34 or 35, this function won't be able to do the job correctly due to the lack of spikes in those two orders, so it will just return the list of wavelengths computed previously in ThAr_calibered
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
- record_path : it is the path of the pkl file created to record our results
Outputs :
- lambdas_Angstrom : list of the new converted_wavelengths
- initial_standard_deviation : the value which tells us about the accuracy of our conversion
- initial_mean : average error between the matched wavelengths (lambda_drs - lambda_atlas)
"""


def compute_order_first_conversion(order,max_detection,precision,converted_lambdas_file_path):
    
    # First of all, we will compute the list of matched wavelengths and indices, using the matching algorithms and the inputs. The wavelengths used for the first computation will always be the original wavelengths recorded in ThAr_calibered_lambdas.pkl.
    
    mtch.order_matching_pickle("ThAr_calibered_lambdas.pkl",order,'temporary_file_for_conversion.pkl',max_detection,precision)
    
    data = mtch.matching_reader('temporary_file_for_conversion.pkl')
    
    spikes_lambdas = data['Spikes_wavelengths']
    spikes_errors = data['Matching_gaps_between_drs_and_atlas']
    spikes_indices = data['Spikes_indices']
    
    # The computation of the initial_standard_deviation is really important because it stands for the overall accuracy of our programm. The main goal of this all script is to minimize this value, so we need to keep this "initial_standard_deviation" in memory because we will compare it to the "new_standard_deviation" after the new conversion and matching computation. From this comparison we will say if the new conversion is better, and when it has reached a stable value, stop the iteration and finally say : we have the better possible conversion for those parameters.
    initial_mean = np.mean(spikes_errors)
    initial_standard_deviation = np.std(spikes_errors)
    print("Initial mean : "+str(initial_mean))
    print("Initial standard deviation : "+str(initial_standard_deviation))
    
    # We can print the errors distribution to have an idea of the precision of our matching.
    plt.figure(15)
    plt.hist(spikes_errors, bins=100,color='red')
    plt.xlabel('Error in Angstrom')
    plt.ylabel('Proportion of the matched wavelengths')
    plt.show()
    
    # Now we have our list of wavelengths and indices, so we can compute the interpolation.
    conversion_results = cnv.arturo_convert_angstroms_order(spikes_lambdas,spikes_indices,order)
    lambdas_Angstrom = conversion_results[0]
    
    # Recording the results 
    conversion_result_file = open(converted_lambdas_file_path,'w')
    data_pickle = pickle.dump(lambdas_Angstrom,conversion_result_file)
    conversion_result_file.close()
    
    if (order == 34 or order == 35 ):

        lambdas_Angstrom = cporder.search_orders("ThAr_calibered_lambdas.pkl","ThAr_calibered_original_intensitites.pkl")[order][0]
    
    # Comparing to the original drs wavelengths with a little plot : 
    lambdas_ThAr_calibered = cporder.search_orders("ThAr_calibered_lambdas.pkl","ThAr_calibered_original_intensitites.pkl")[order][0]
    indices_list = [i for i in range(0,order_total*order_len)]
    order_lambdas_indices = indices_list[order_len*order:order_len*(order+1)]
    plt.figure(18)
    plt.title("Comparison : black = calibered lambdas, red = interpolated lambdas")
    plt.plot(order_lambdas_indices,lambdas_ThAr_calibered,'b+-')
    plt.plot(order_lambdas_indices,lambdas_Angstrom,'r--')
    plt.show()
    
    return(lambdas_Angstrom,initial_standard_deviation,initial_mean)


"""
The following function will use the validation_methods module's scripts to generate the list of matching wavelengths we need to achieve the interpolation, for all orders. It will also plot some usefull informations. Then it will use this matching wavelengths list to compute a new conversion and record the result (which is a new wavelengths list from which we can iterate the process).
Inputs :
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
- record_path : it is the path of the pkl file created to record our results
Outputs :
- lambdas_Angstrom : list of the new converted_wavelengths
- initial_stds : the list of the values (standard_deviaitons) which tells us about the accuracy of our conversion
"""


def compute_first_conversion(max_detection, precision, converted_lambdas_file_path):
    
    # Creation of the lists containing the global informations that we will compute
    lambdas_Angstrom = []
    initial_stds = []
    
    # For each order, we will record the converted wavelengths in a different file with the number of the order at the end
    for order in range(36):
        
        # Defining the record file's name
        record_file = converted_lambdas_file_path + "_order_"+str(max_detection)+"_"+str(precision)
        
        # Computing the results 
        order_results = compute_order_first_conversion(order,max_detection,precision,record_file)
        
        order_converted_lambdas = order_results[0]
        order_initial_std = order_results[1]
        
        lambdas_Angstrom.insert(order,order_converted_lambdas)
        initial_stds.insert(order,order_initial_std)
        
        # Recording the results 
        pickle.dump(initial_stds,"initial_stds""_order_"+str(max_detection)+"_"+str(precision))
    
    return(lambdas_Angstrom, initial_stds)



    
    
"""
This version concerns only a particular order, we'll extend it to all orders after that.
Now we have initiated the process by computing the first conversion and recording the first standard_deviation, we need to compute the validation again and see if the standard_deviation has been improved or not, and use the matching results to compute another conversion to iterate the process if need be. The following function will do this job, using a converted list of wavelengths to compute the corresponding standard_deviation. Then it will compute a new conversion and record the wevelengths in a pickle file.
Inputs :
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use) 
- order : int, number of the order which we want to compute
/!\ If order = 34 or 35, this function won't be able to do the job correctly due to the lack of spikes in those two orders, so it will just return the list of wavelengths computed previously in ThAr_calibered
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
- new_converted_lambdas_path : path to the file where the converted wavelengths will be recorded
Outputs :
- lambdas_Angstroms : new converted list of wavelengths
- new_standard_deviation : new computed standard deviation
- new_mean : new average error beteween drs and atlas in Angstroms
"""

def compute_order_conversion(lambdas_path, order, max_detection, precision, new_converted_lambdas_path):
    
    # Let's compute the matching data : spikes indices and wavelengths, erros, etc in order to get the necessary data for our polynomial inteprolation.
    
    mtch.order_matching_pickle(lambdas_path,order,"temporary_file_for_order_conversion.pkl",max_detection,precision)
    data = mtch.matching_reader("temporary_file_for_order_conversion.pkl")
    
    spikes_lambdas = data['Spikes_wavelengths']
    spikes_errors = data['Matching_gaps_between_drs_and_atlas']
    spikes_indices = data['Spikes_indices']
    
    # Now, we've all we need to start the new conversion!
    
    # Let's compute the new standard deviation.
    
    new_mean = np.mean(spikes_errors)
    new_standard_deviation = np.std(spikes_errors)
    print("Current mean : "+str(new_mean))
    print("Current standard deviation : "+str(new_standard_deviation))
    
    # We can print the errors distribution to have an idea of the precision of our matching.
    plt.figure(15)
    plt.hist(spikes_errors, bins=100,color='red')
    plt.xlabel('Error in Angstrom')
    plt.ylabel('Proportion of the matched wavelengths')
    plt.show()
    
    # We can also compute the error in m/s and record it in a pickle
    
    errors_file = open("results/Calibration errors/calibration_errors_"+str(order)+"_"+str(max_detection)+"_"+str(precision),'w')
    calibration_errors = [ 3*10**8*spikes_errors[i]/spikes_lambdas[i] for i in range(len(spikes_errors)) ]
    errors_data = [spikes_lambdas,calibration_errors]
    pickle.dump(errors_data,errors_file)
    errors_file.close()
    
    # Now we have our list of wavelengths and indices, so we can compute the interpolation.
    conversion_results = cnv.arturo_convert_angstroms_order(spikes_lambdas,spikes_indices,order)
    lambdas_Angstrom = conversion_results[0]
    order_coeffs = conversion_results [1]
    
    # Recording the results 
    conversion_result_file = open(new_converted_lambdas_path,'w')
    data_pickle = pickle.dump(lambdas_Angstrom,conversion_result_file)
    conversion_result_file.close()
    
    coeffs_file = open("results/Interpolation coefficients/Interpolation_coefficients_order_"+str(order)+"_"+str(max_detection)+"_"+str(precision),'w')
    pickle.dump(order_coeffs,coeffs_file)
    coeffs_file.close()
    
    # Dealing with the two lasts orders, which don't have any data to interpolate for now.
    
    if (order == 34 or order == 35 ):

        lambdas_Angstrom = cporder.search_orders("ThAr_calibered_lambdas.pkl","ThAr_calibered_original_intensitites.pkl")[order][0]
    
    # Comparing to the original drs wavelengths with a little plot : 
    
    lambdas_ThAr_calibered = cporder.search_orders("ThAr_calibered_lambdas.pkl","ThAr_calibered_original_intensitites.pkl")[order][0]
    indices_list = [i for i in range(0,order_total*order_len)]
    order_lambdas_indices = indices_list[order_len*order:order_len*(order+1)]
    plt.figure(18)
    plt.title("Comparison : black = calibered lambdas, red = interpolated lambdas")
    plt.plot(order_lambdas_indices,lambdas_ThAr_calibered,'b+-')
    plt.plot(order_lambdas_indices,lambdas_Angstrom,'r--')
    plt.show()
    
    return(lambdas_Angstrom,new_standard_deviation,new_mean)


"""
This version compute all orders one after the other.
Now we have initiated the process by computing the first conversion and recording the first standard_deviation, we need to compute the validation again and see if the standard_deviation has been improved or not, and use the matching results to compute another conversion to iterate the process if need be. The following function will do this job, using a converted list of wavelengths to compute the corresponding standard_deviation. Then it will compute a new conversion and record the wevelengths in a pickle file.
Inputs :
- lambdas_path : path to the file containing the wavelengths list which will be used to compute the validation scripts (since we are gonna compare different wavelengths lists to know if the conversion has been good, we need to indicate which wavelengths we want to use) 
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
- new_converted_lambdas_path : path to the file where the converted wavelengths will be recorded
Outputs :
- new_standard_deviation : new computed standard deviation
- new_lambdas : new list of converted wavlengths
"""

def compute_all_order_conversion(lambdas_path, max_detection, precision, new_converted_lambdas_path):
    
    # Creation of the lists containing the global informations that we will compute
    lambdas_Angstrom = []
    new_stds = []
    
    # For each order, we will record the converted wavelengths in a different file with the number of the order at the end.
    for order in range(36):
        
        # Defining an accurate name for the pkl file
        record_file = new_converted_lambdas_path + "_order_ " + str(order) +str(max_detection)+"_"+str(precision)+"_"+str(stop_value)
        
        # Computing the results
        order_results = compute_order_first_conversion(order,max_detection,precision,record_file)
        
        order_converted_lambdas = order_results[0]
        order_new_std = order_results[1]
        
        lambdas_Angstrom.insert(order,order_converted_lambdas)
        new_stds.insert(order,order_new_std)
        
        # Recording the results 
        pickle.dump(order_converted_lambdas,record_file)
    
    return(lambdas_Angstrom, new_stds)
    
    
"""  
Now we have everything we need to build the iteration : we can initiate it with compute_first_conversion and then use compute_conversion to compare the successive strandard deviations and decide if we iterate or not.
The following function will do this job and implement the iteration process.
Inputs :
- order : int, number of the order 
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
Outputs :
- std_evolution : list containing the successives computed standard_deviation
- mean_evolution : list containing the successives computed means
"""

def order_iterated_conversion(order,max_detection,precision,stop_value,final_converted_lambdas_path):
    print(" ")
    print(" _________________________ Initial conversion _________________________ ")
    print(" ")
    
    # Defining the detailed record file's name 
    record = "results/Converted wavelengths"+final_converted_lambdas_path+"_"+str(order)+"_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value)
    
    first_results = compute_order_first_conversion(order,max_detection,precision,"results/first_converted_lambdas.pkl")
    old_lambdas = first_results[0]
    old_standard_deviation = first_results[1]
    old_mean = first_results[2]
    std_evolution = [old_standard_deviation]
    mean_evolution = [old_mean]
    print(" ")
    print(" _________________________ First Iteration _________________________ ")
    print(" ")
    
    # Initiating the iteration loop
    
    new_results = compute_order_conversion('results/first_converted_lambdas.pkl',order,max_detection,precision,record)
    new_lambdas = new_results[0]
    new_standard_deviation = new_results[1]
    new_mean = new_results[2]
    std_evolution.append(new_standard_deviation)
    mean_evolution.append(new_mean)
    
    # First record of the results in a pkl file :
    std_file = open("results/standard_deviation_evolution_"+str(order)+"_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),'w')
    pickle.dump(std_evolution,std_file)
    std_file.close()
    
    mean_file = open("results/mean_evolution_"+str(order)+"_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),'w')
    pickle.dump(mean_evolution,mean_file)
    mean_file.close()
    
    # This variable is an iteration limit to make sure the programm ends
    stop = 2
    
    while ( (abs(old_standard_deviation - new_standard_deviation) > stop_value) and (stop != 100) ):
        
        print(" ")
        print(" _________________________ "+" New iteration "+str(stop)+" _________________________ ")
        print(" ")
        
        # Updating the old parameters with the current ones 
        
        old_standard_deviation = new_standard_deviation
        old_mean = new_mean
        old_lambdas = new_lambdas
        
        # Computing the new results
        
        new_results = compute_order_conversion(record,order,max_detection,precision,record)
        
       
        # Updating the iteration's parameter with the new computed ones
        
        new_lambdas = new_results[0]
        new_standard_deviation = new_results[1]
        new_mean = new_results[2]
        
        std_evolution.append(new_standard_deviation)
        mean_evolution.append(new_mean)
        
        # Recording the results in a pkl file :
        std_file = open("results/standard_deviation_evolution_"+str(order)+"_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),'w')
        pickle.dump(std_evolution,std_file)
        std_file.close()
        
        mean_file = open("results/mean_evolution_"+str(order)+"_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),'w')
        pickle.dump(mean_evolution,mean_file)
        mean_file.close()
        
    # Updating the emergency stop variable
        stop += 1
    
    print("Standard deviation evolution :")
    print(std_evolution)
    print("Mean evolution :")
    print(mean_evolution)
    
    print("---- Iteration completed ----")
    
    #Temporary for easier tests
    plt.close(3)
    plt.close(15)
    plt.close(4)
    
    return(std_evolution,mean_evolution)
    
    
"""  
Now we have everything we need to build the iteration : we can initiate it with compute_first_conversion and then use compute_conversion to compare the successive strandard deviations and decide if we iterate or not.
The following function will do this job and implement the iteration process.
The version below is the same as order_iterated_conversion but for every order (except 34 and 35 as because it makes no sense to match without converting after...), so it takes something like an eternity to compute, do it only if you have several hours or even days!
Inputs :
- max_detection : float representing the minimum value of the detected maxima in the compute_spikes algorithm.
- precision : when comparing the atlas to the drs data, we will consider that two close wavelengths are matching only if they are closer than precision Angstroms.
Outputs :
- std_evolutions : list containing the successives computed standard deviation evolutions lists (one per order)
- mean_evolutions : list containing the successives computed mean evolutions lists (one per order)
"""

def all_orders_iterated_conversion(max_detection,precision,stop_value,iteration_converted_lambdas):
    
    # Defining a more detailed name :
    iteration_converted_lambdas = iteration_converted_lambdas+str(max_detection)+"_"+str(precision)+"_"+str(stop_value)
    
    print(" ")
    print(" ========== STARTING THE FULL CONVERSION ========== ")
    
    print(" ")
    # We will record the different standard deviation evolutions (one per order) in a big list which is defined just below.
    std_evolutions = []
    mean_evolutions = []
    
    # We will record each lambdas_list in different file which will have a name depending on the order with a subtil "_orderN" (where N is the number of the order, obviously) at the end. The oly thing left to do in this function is to call for each of those names the function which iterates an order. Let's go :
    
    # Since orders 34 and 35 are not converted (because of the lack of spikes in the ThAr red wavelengths), we will stop at 33...
    
    for order in range(34):
        
        file_name = iteration_converted_lambdas+"_order_"+str(order)
        order_results = order_iterated_conversion(order,max_detection,precision,stop_value,file_name)
        std_evolutions.insert(order,order_results[0])
        mean_evolutions.insert(order,order_results[1])
        
        
        # Recording the results in a pkl file :
        std_evolutions_file = open("std_evolutions_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),"w")
        pickle.dump(std_evolutions,std_evolutions_file)
        std_evolutions_file.close()
        
        mean_evolutions_file = open("mean_evolutions_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),"w")
        pickle.dump(mean_evolutions,mean_evolutions_file)
        mean_evolutions_file.close()
    
    
        # Recording the results in a pkl file :
        std_evolutions_file = open("std_evolutions_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),"w")
        pickle.dump(std_evolutions,std_evolutions_file)
        std_evolutions_file.close()
        
        mean_evolutions_file = open("mean_evolutions_"+str(max_detection)+"_"+str(precision)+"_"+str(stop_value),"w")
        pickle.dump(mean_evolutions,mean_evolutions_file)
        mean_evolutions_file.close()
        
    return(std_evolutions,mean_evolutions)
    
    
    


"""
Some of the values used in the ThAr Atlas pkl can't be used for the comparison because they match with wavelengths which are blinded in the ThAr calibered spectrum. We won't use them so we will modify the ThAr pkl to avoid those problematic values. Tis function will be used to do those modifications. It will be done order per order empically, using the other algorithms. The new list (which we will record with a pkl), will be used to compute the new matching and improve the results.
"""

def ThAr_pickler():
    
    name = "ThAr_Atlas_filtered.pkl"
    file = open(name,'r')
    Th_Ar_list = pickle.load(file)
    
    print(Th_Ar_list)
    
    print("Which wavelengths should be deleted?")
    file.close()
    
    
    deleted_lambda = input()
    
    file = open("ThAr_Atlas_filtered.pkl",'w')
    if ( not Th_Ar_list.__contains__(deleted_lambda) ):
        print("Error : this wavelengths is not even in the list!")
    Th_Ar_list.remove(deleted_lambda)
    
    pickle.dump(Th_Ar_list,file)
    file.close()
    print("Modified!")
    



"""
Little script to plot everything we have computed before
"""
"""
means = []
stds = []

for i in range(0,34) :
    
    plt.figure(110)
    plt.title("Std and Mean evolutions for the order "+str(i)+ "red = Mean, blue = Std")
    
    mean_file = open("results/Mean evolution/mean_evolution_"+str(i)+"_0.1_0.1_1e-10",'r')
    mean_evolution = pickle.load(mean_file)
    plt.plot(mean_evolution,'r+')
    
    
    std_file = open("results/Std evolution/standard_deviation_evolution_"+str(i)+"_0.1_0.1_1e-10",'r')
    std_evolution = pickle.load(std_file)
    plt.plot(std_evolution, 'b+')
    
    means.insert(i,mean_evolution[-1])
    stds.insert(i,std_evolution[-1])
    
    plt.show()

plt.figure(100)
plt.title("Average error for each order after iteration convergence")
plt.plot(means,'b+')
plt.figure(101)
plt.title("Standard deviation for each order's error after iteration convergence")
plt.plot(stds,'b+')
plt.show()
"""
"""
average_lambdas = []

for i in range(34) : 
    lambda_file = open("results/Converted wavelengths/Converted wavelengthstest0.1_0.1_1e-08_order_"+str(i)+"_"+str(i)+"_0.1_0.1_1e-08",'r')
    lambdas = pickle.load(lambda_file)
    lambda_file.close()
    average_lambdas.insert(i,np.mean(lambdas))


errors_file = open("mean_evolutions_0.1_0.1_1e-08",'r')
errors = pickle.load(errors_file)
errors_file.close()
final_errors = [errors[i][-1] for i in range(len(errors))]

radial_vel = [final_errors[i]/average_lambdas[i]*3*10**8 for i in range(len(final_errors))]
# plt.figure(56)
# plt.plot(final_errors,'ro')
plt.figure(200)
plt.plot(average_lambdas,radial_vel,'ro')
plt.xlabel("Wavelength (Angtroms)")
plt.ylabel("Absolute calibration error (m/s)")
plt.title("Evolution of calibration error after computation of the errors between Atlas and DRS (erreur moyenne et lambda moyen Angtroms par ordre puis vitesse)")
plt.show()
"""

"""
plt.figure(31)
plt.title('Calibration errors ("ThAr_calibered.fits" compared to ThAr Atlas)  ')

for i in range(34):
    order_data_file = open("results/Calibration errors/calibration_errors_"+str(i)+"_0.1_0.1",'r')
    data = pickle.load(order_data_file)
    order_data_file.close()
    
    lambdas = data[0]
    errors = data[1]
    plt.plot(lambdas,errors,'r+')
plt.xlabel("Wavelength (Angstroms)")
plt.ylabel("Calibration error (m/s)")
plt.show()

plt.figure(32)
plt.title("Evolution of calibration error after computation of the errors between Atlas and DRS (moyenne sur les erreurs en m/s)")
average_lambdas = []
average_errors = []
for i in range(34):
    order_data_file = open("results/Calibration errors/calibration_errors_"+str(i)+"_0.1_0.1",'r')
    data = pickle.load(order_data_file)
    order_data_file.close()
    
    lambdas = data[0]
    errors = data[1]
    
    average_lambdas.insert(i,np.mean(lambdas))
    average_errors.insert(i,np.mean(errors))
plt.plot(average_lambdas,average_errors,'ro')
plt.xlabel("Wavelength (Angstroms)")
plt.ylabel("Calibration error (m/s)")
plt.show()
"""


"""
all_coeffs = []
for i in range(34):
    
    coeffs_file = open("results/Interpolation coefficients/Interpolation_coefficients_order_"+str(i)+"_"+str(0.1)+"_"+str(0.1),'r')
    
    coeffs  = pickle.load(coeffs_file)
    print(coeffs)
    coeffs_file.close()
    all_coeffs.insert(i,coeffs)
    plt.figure(40)
    plt.title("Degre du coefficient : "+str(5))
    plt.plot(i,coeffs[0],'r+')
    plt.figure(41)
    plt.title("Degre du coefficient : "+str(4))
    plt.plot(i,coeffs[1],'b+')
    plt.figure(42)
    plt.title("Degre du coefficient : "+str(3))
    plt.plot(i,coeffs[2],'g+')
    plt.figure(43)
    plt.title("Degre du coefficient : "+str(2))
    plt.plot(i,coeffs[3],'r.')
    plt.figure(44)
    plt.title("Degre du coefficient : "+str(1))
    plt.plot(i,coeffs[4],'b.')
    plt.figure(45)
    plt.title("Degre du coefficient : "+str(0))
    plt.plot(i,coeffs[5],'g.')
    plt.show()
"""
"""
data = read_fits("th_calibered.fits")
lambdas_artur = data[0]

lambdas_all_orders = []

for i in range(34):
    
    name = "results/Converted wavelengths/Converted wavelengthstest0.1_0.1_1e-08_order_"+str(33-i)+"_"+str(33-i)+"_0.1_0.1_1e-08"
    file = open(name,'r')
    lambdas = pickle.load(file)
    file.close()
    for j in range(len(lambdas)):
        lambdas_all_orders.append(lambdas[j])

plt.plot(lambdas_artur,color = 'black')

for i in range(2*order_len):
    lambdas_all_orders.insert(i,lambdas_artur[i])


plt.plot(lambdas_all_orders,color='red')
plt.show()

file1 = open("Wavelengths_corrected_after_iteration",'w')
pickle.dump(lambdas_all_orders,file1)
file1.close()
"""
    


























    