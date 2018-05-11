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
import calibration_methods.itered_calibration as itrd

# Definition of global values
global order_len # length of an order
order_len = 4612
global order_total # number of orders
order_total = 36

"""
We will use the poly2D_fit and the itered_calibration algorithms to compute a new calibration, using our fitted coefficients for the conversion polynom. Those coefficients come from the iterated conversions and have been fitted and computed again to try to improve their precision. 
The idea is to replace the actual conversion method, which uses the interpolated_conversion algorithms by a new one using th fitted coefficients, and use it in the itered_calibration alogrithms, so we will use previously made code and modify it a little.
"""


"""
At this point, we have computed the best possible precision order per order (the best we currently can do with our current ideas obvously). Now we will try to see if we can improve this precision of calibration by using 2D fits of the coefficients of the polynom which is used to convert wavelengths from Arturos to Angstroms. For now, we consider that each polynom used for the conversion only depends on the order, but maybe we can find a tendance, a law which gives an idea of the evolution of the polynomial coefficents according to the order. In other words, the parameter will now be the order number and the data to fit the coefficients. We will find a new law for each coefficient (each degree), giving the coefficients. We will use those to create the corresponding polynom and try to convert and match the original wavelengths with that. Then we will be able to tell if it is a way to improve the results or not.
The structure of the lists of coefficients used will be the following : a list of lists which will have for each order : [a0, a1,a2,etc] where ai is the coeff of degree P-i where P is the degree of the polynom. (higher degree first then the other, and degree 0 at the end of each list).
Inputs :
None.
Outputs :
- lists_of_fitted_coeffs : list of the fitted coeffs, order per order
"""

def polyfit_dispersion2D():
    
    N = 34 # How many orders do you want to fit (the N first orders will be fitted)
    
    # Order of the polynomial fit :
    p = 20
    
    lists_of_coeffs = []
    list_of_fitted_coeffs = []
    
    # Creating 6 lists of coefficients to fit
    list_of_coeffs0 = []
    list_of_coeffs1 = []
    list_of_coeffs2 = []
    list_of_coeffs3 = []
    list_of_coeffs4 = []
    list_of_coeffs5 = []
    
    # Importing the coefficients from the previously recorded pickles in our lists to fit
    
    for i in range(N):
    
        coeffs_file = open("results/Interpolation coefficients/Interpolation_coefficients_order_"+str(i)+"_"+str(0.1)+"_"+str(0.1),'r')
        old_coeffs  = pickle.load(coeffs_file)
        coeffs_file.close()

        list_of_coeffs0.insert(i,old_coeffs[0])
        list_of_coeffs1.insert(i,old_coeffs[1])
        list_of_coeffs2.insert(i,old_coeffs[2])
        list_of_coeffs3.insert(i,old_coeffs[3])
        list_of_coeffs4.insert(i,old_coeffs[4])
        list_of_coeffs5.insert(i,old_coeffs[5])
    
    lists_of_coeffs.insert(0,list_of_coeffs0)
    lists_of_coeffs.insert(1,list_of_coeffs1)
    lists_of_coeffs.insert(2,list_of_coeffs2)
    lists_of_coeffs.insert(3,list_of_coeffs3)
    lists_of_coeffs.insert(4,list_of_coeffs4)
    lists_of_coeffs.insert(5,list_of_coeffs5)
    
    
    indices = [i for i in range(N)]
    
    lists_of_fitted_coeffs = [ [] for order in range(N) ]
    
    for i in range(6):
        
        
        try :    
        
            coeffs = np.polyfit(indices,lists_of_coeffs[i],p)
            
        except :
            print("Polynomial fitting of the coefficients failed!")
                        
        # Computation of the fit
        pol = np.poly1d(coeffs)
        fitted_coeffs = pol(indices)
        plt.figure(60+i)
        plt.title("Interpolation results")
        plt.plot(indices,lists_of_coeffs[i],'o',color='black')
        plt.plot(indices,fitted_coeffs,'o',color='purple')
        plt.show()
        
        for order in range(N):
            lists_of_fitted_coeffs[order].insert(i,fitted_coeffs[order])
    
    record_file = open("fitted_coeffs_degree"+str(p),'w')
    pickle.dump(lists_of_fitted_coeffs,record_file)
    record_file.close()
    
    return(lists_of_fitted_coeffs)
    
    
    
    
    

"""
This function will convert one order from Arturos to Angstroms, using the fitted coeffs computed thanks to poly2D_fit.
Inputs :
- order : number of the order to convert
- lists_of_fitted_coeffs : list of the whole coeffs after fitting by the polyfit_dispersion2D algorithm
Outputs :
- order_lambdas_Angstroms : list of converted wavelengths in Angstroms, using the 2D fit.
""" 

def convert_arturo_angstroms2D(order,lists_of_fitted_coeffs):
    
    print(" ")
    print(" _________________________ Converting  _________________________ ")
    print(" ")
    
    # Selecting the rigth coeffs in the input list
    order_fitted_coeffs = lists_of_fitted_coeffs[order]

    
    # Generating the arturo scaled wavelengths list to convert
    arturos_list = [i for i in range(0,order_total*order_len)]
    order_lambdas_arturos = arturos_list[order_len*order:order_len*(order+1)]
    
    # Creating the output list
    order_lambdas_Angstroms = []
    
    # Convertion using the 2Dfit's results (the coeffs) :
    #print("Polynomial coefficients",order_fitted_coeffs)
    pol = np.poly1d(order_fitted_coeffs)
    order_lambdas_Angstroms = pol(order_lambdas_arturos)

    return(order_lambdas_Angstroms)
   
   
    
"""
Now we have written the conversion function so we can basically use it the way we used the interpolated_conversion function. We are gonna write the algorithm which will use this function, convert all orders wavelengths, use those wavelengths and compute a matching to see if there is an improvevment or not.
"""


"""
The following function will compute a conversion using the function above, and then compute the matching or a given order. It will plot the results so that we can compare the effiency of the 2D fit versus the itered conversion. 
Inputs : 
- order : the number of the order to compute
Outputs : 
- None
"""

def order_conversion2D(order):
    
    # Computing all the fitted coefficients
    fitted_coeffs = polyfit_dispersion2D()
    
    # Converting the wavelengths for the given order
    order_lambdas = convert_arturo_angstroms2D(order,fitted_coeffs)
    print(order_lambdas)
    lambdas_file = open("temporary_file_for_2Dconversion",'w')
    pickle.dump(order_lambdas,lambdas_file)
    lambdas_file.close()
    
    matching_results = mtch.order_matching("temporary_file_for_2Dconversion",0.1,order,0.1)
    
    return(None)


def all_order_conversion2D(n):
    for i in range(n):
        order_conversion2D(i)
        


"""
The following function computes the polynomial coefficients of the conversion polynom for the given order.
Inputs :
- order : the order
Outputs :
None.
"""
def coeffs(order):
    
    alpha=63.495*np.pi/180.
    gamma=0.6*np.pi/180.
    G=79. # grooves/mm
    F=388.  #mm focal length
    p=13.5e-3 # 12 micron pixel in mm

    m = 57 -order
    
    # On multiplie par 1e7 pour passer des mm aux Angstroms
    a0 = 1e7*(1.0/m)*2*np.cos(gamma)*np.sin(alpha)/G
    a1 =  1e7*(1.0/m)*np.cos(gamma)*np.cos(alpha)*p/G/F
    a2 =  1e7*-(1.0/m)*np.cos(gamma)*np.sin(alpha)*p**2/2/G/(F**2)
    a3 =  1e7*-(1.0/m)*np.cos(gamma)*np.cos(alpha)*p**3/6/G/(F**3)
    a4 =  1e7*(1.0/m)*np.cos(gamma)*np.sin(alpha)*p**4/24/G/(F**4)
    a5 =  1e7*(1.0/m)*np.cos(gamma)*np.cos(alpha)*p**5/120/G/(F**5)
    
    print(a5,a4,a3,a2,a1,a0)
    
    
    coeffs_file = open("results/Interpolation coefficients/Interpolation_coefficients_order_"+str(order)+"_"+str(0.1)+"_"+str(0.1),'r')
    old_coeffs  = pickle.load(coeffs_file)
    coeffs_file.close()
    print(old_coeffs)



"""
Function which cleans the screen by closing all the plots.
"""
def clean_plt():
    for i in range(100):
        plt.close()



"""
Another try of a different fit :
2D fit of the f(m,indice) : ex : f(1,10) = 0.5[f(0,10)+f(2,10)]
"""

def grid2D():
    
    # Creating the grid
    param = (34,order_len)
    grid = np.zeros(param)
    
    # Filling the grid
    for i in range(34):
        
        # Loading the better polynomial fit for each order
        coeffs_file = open("results/Interpolation coefficients/Interpolation_coefficients_order_"+str(i)+"_"+str(0.1)+"_"+str(0.1),'r')
        old_coeffs  = pickle.load(coeffs_file)
        coeffs_file.close()
        #print(old_coeffs)
        
        order_polynom = np.poly1d(old_coeffs)
        indices = [k for k in range(order_len*i , order_len*(i+1))]
        for j in range(order_len):
            grid[i,j] = order_polynom(indices[j])
    
    print(grid)
    # Now we have that big matrix containing all the informations about the conversions of all orders. We can easily use it to interpolate between the orders and find a new "vertical" fit. The result will be a 2D cross fit between the 1D horizontal fit for each order and the 1D vertical fit between the orders.
     
    orders = [o for o in range(34)]
    
    new_grid = np.zeros(param)
    
    for ind in range(order_len):
        
        vertical_values = [ grid[order,ind] for order in orders ] # Creating the list of values to fit for each ind
        
        try :    
            coeffs = np.polyfit(orders,vertical_values,10)
            
        except :
            print("Polynomial fitting failed!")
                        
        # Computation of the fit
        pol = np.poly1d(coeffs)
        
        for o in orders : 
            new_grid[o,ind] = pol(o)  
    
    print(new_grid)
    # Now the new grid has been computed, we can compute a new matching for each order.
    for o in orders :
        
        order_lambdas = [ new_grid[o,ind] for ind in range(order_len) ]
        order_lambdas_file = open("temporary_file_for_vertical_fit",'w')
        pickle.dump(order_lambdas,order_lambdas_file)
        order_lambdas_file.close()
        order_matching_results = mtch.order_matching("temporary_file_for_vertical_fit",0.1,o,0.1)
    
    return(None)    
        
    