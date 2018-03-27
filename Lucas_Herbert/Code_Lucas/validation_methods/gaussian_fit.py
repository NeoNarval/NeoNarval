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
This function finds the gaussian function which is the closest from our data. It's a fitting algorithm which uses a least squares method. It returns the centre and the width of the fitted gaussian. It takes an an input the necessary data for the fit : the wavelengths, lambdas, and their intensities (or whatever we need to fit), data. It returns the centre and width of the fitted gaussian.
"""

    
def fit_the_spike(lambdas,data):
    
    y = np.copy(data) - np.min(data)
    X = lambdas
    def gaussian(x,cen,amp,wid):
        return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
    
    # printing the initial gaussian (before the optimization of the fit)
    naive_center = float(np.sum(lambdas*y))/np.sum(y)
    naive_width = np.sqrt(abs((np.sum((lambdas-naive_center)**2*y)/np.sum(y))))
    naive_ampl = np.max(y)
    naive_gaussian = [ gaussian(x,naive_center,naive_ampl,naive_width) for x in lambdas ]
    # plt.plot(lambdas,naive_gaussian,color='green')
    # plt.axvline(naive_center,color='green')
    # plt.show()
    lambda_centre = naive_center
    lambda_width = naive_width
    try :
        # we use the lmfit algorithm to improve our fit's precision
        gaussian_model = lmfit.Model(gaussian)
        params = gaussian_model.make_params(cen=naive_center,amp=naive_ampl,wid=naive_width)
        result = gaussian_model.fit(y,params,x=X) 
        # printing the best gaussian fit
        best_gaussian_fit = result.best_fit
        best_cen = result.best_values['cen']
        best_wid = result.best_values['wid']
        best_amp = result.best_values['amp']
        best_params_gaussian = [ gaussian(x,best_cen,best_amp,best_wid) + np.min(data) for x in lambdas ]
        plt.figure(1)
        plt.plot(lambdas, best_params_gaussian, 'b--', color='purple')
        plt.plot(best_cen,best_amp,'.',color='purple')
        plt.show()
        #computed_centre = float(np.sum(X*best_gaussian_fit))/np.sum(best_gaussian_fit) 
        #plt.plot(lambdas, best_gaussian_fit, 'b--' ,color='purple')
        #plt.axvline(computed_centre, color='purple')
        lambda_centre = best_cen
        lambda_width = best_wid
        # we need the report data to improve our understanding of the results
        report = result.fit_report()
        
        
    except : 
        report = "Computation failed for this spike : default data = naive fit"
        pass
        
    return(lambda_centre,lambda_width,report)