import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

"""
This function finds the gaussian function which is the closest from our data. It's a fitting algorithm which uses a least squares method. It returns the centre and the width of the fitted gaussian.
"""


def gaussian_fit(lambdas,data):

    y = np.copy(data)
    X = np.copy(lambdas)
    
    def gaussian(x,cen,amp,wid):
        return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
    
    # We can compute the intitial gaussian, which is the "naïve" way to fit the data with a gaussian
    naive_center = float(np.sum(lambdas*y))/np.sum(y)
    naive_width = float(np.sqrt(abs((np.sum((lambdas[i]-naive_center)**2*y[i] for i in range(len(lambdas)))/np.sum(y)))))
    naive_ampl = np.max(y)
    naive_gaussian = [ gaussian(x,naive_center,naive_ampl,naive_width) for x in lambdas ]
    # plt.plot(lambdas,naive_gaussian,color='green')
    # plt.axvline(naive_center,color='green')
    # plt.show()
    
    # if the following fit fails, we can keep this naïve data as the fitted parameters :
    lambda_centre = naive_center  
    lambda_width = naive_width
    
    
    try :
        # we use the lmfit algorithm to improve our fit's precision
        gaussian_model = lmfit.Model(gaussian)
        params = gaussian_model.make_params(cen=naive_center,amp=naive_ampl,wid=naive_width)
        result = gaussian_model.fit(y,params,x=X) 
        
        # computing the best gaussian fit parameters
        best_gaussian_fit = result.best_fit
        best_cen = result.best_values['cen']
        best_wid = result.best_values['wid']
        best_amp = result.best_values['amp']
        
        # computing and plotting the associated gaussian
        best_params_gaussian = [ gaussian(x,best_cen,best_amp,best_wid) for x in lambdas ]
        plt.figure(1)
        plt.plot(lambdas, best_params_gaussian, color='red')
        plt.axvline(best_cen, color = 'red')
        plt.show()
        
        #computed_centre = float(np.sum(X*best_gaussian_fit))/np.sum(best_gaussian_fit) 
        #plt.plot(lambdas, best_gaussian_fit, 'b--' ,color='purple')
        #plt.axvline(computed_centre, color='purple')
        
        lambda_centre = best_cen
        lambda_width = best_wid
        
        # we can save the report of the fit to improve our understanding of the results
        report = result.fit_report()
        
        
    except : 
        report = ""
        pass
        
    return(lambda_centre,lambda_width)