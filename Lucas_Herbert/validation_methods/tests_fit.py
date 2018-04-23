#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Python's modules imports :
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit



"""
This function tries to fit the input data with a function like "ax+b".
Inputs :
- lambdas : list of wavelengths
- data : list of associated values
Output : 
- best_a : float, coefficient of the fitted affine function.
- best_b : float, coefficient of the fitted affine function.
"""

def fit_affine(lambdas,data):
    
    X = np.copy(lambdas)
    y = np.copy(data)
    plt.figure(3)
    plt.plot(X,y,'.',color='blue',label='Original data')
    plt.show()
    
    best_a = 0
    best_b = 0
    
    def affine(x,a,b):
        return(a*x+b)
        
    try :
        affine_model = lmfit.Model(affine)
        params = affine_model.make_params(a=0,b=0.02)
        result = affine_model.fit(y,params,x=X)
        
        best_a = result.best_values['a']
        best_b = result.best_values['b']
        
        best_affine = [ affine(x,best_a,best_b) for x in X]
        plt.figure(3)
        plt.plot(X,best_affine,color='red',label='Best affine fit')
        plt.show()
        report = result.fit_report()
        print(report)
    
    except :
        print("Data affine fitting failed")
        pass
    
    return(best_a,best_b)
    
"""
This functions tries to fit the input data with a degree 2 polynome.
Inputs :
- lambdas : list of wavelengths
- data : list of associated values
Output : 
- best_a : float, coefficient of the fitted polynomial function.
- best_b : float, coefficient of the fitted polynomial function.
- best_c : float, coefficient of the fitted polynomial function.
"""

def fit_pol2(lambdas,data):
    
    X = np.copy(lambdas)
    y = np.copy(data)
    
    # plt.figure(3)
    # plt.plot(X,y,'.',color='blue',label='Original data')
    # plt.show()
    
    best_a = 0
    best_b = 0
    best_c = 0
    
    def pol2(x,a,b,c):
        return(a*x*x+b*x+c)
    
    try :
        model = lmfit.Model(pol2)
        params = model.make_params(a=0,b=0,c=0)
        result = model.fit(y,params,x=X)
        
        best_a = result.best_values['a']
        best_b = result.best_values['b']
        best_c = result.best_values['c']
            
        best_pol2 = [ pol2(x,best_a,best_b,best_c) for x in X]
        plt.figure(3)
        plt.plot(X,best_pol2,color='red',label='Best polynomial fit')
        plt.show()
        report = result.fit_report()
        print(report)

    
    except :
        print("Data 2d degree polynomial fitting failed")
        pass

    return(best_a,best_b,best_c)
    
    
    
    
    
    
    


