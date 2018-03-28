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
"""

def fit_affine(lambdas,data):
    
    X = np.copy(lambdas)
    y = np.copy(data)
    plt.figure(3)
    plt.plot(X,y,'.',color='blue',label='Original data')
    plt.show()
    
    def affine(x,a,b):
        return(a*x+b)
        
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
    return(best_a,best_b)