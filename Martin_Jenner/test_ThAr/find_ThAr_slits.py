
#-*- coding: utf-8 -*-

# Python's global modules imports :
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit
import collections
from scipy.optimize import leastsq
from scipy.optimize import curve_fit


""" 
In this script, we will search for all the spikes of a normalized spectrum and then fit a gaussian to each of the interesting spikes in order to find its precise centre and width.
"""
def rint(x): return int(round(x))

def sum_on_y_axis(left_lane, thickness, ini_index, end_index):
    lenX = end_index - ini_index
    intensity = np.zeros(lenX)
    nb_col = rint(thickness)
   
    for i in range(ini_index, end_index):
        intensity[i-ini_index] = sum(ThAr_data[i][min(rint(left_lane[i]+j), lenDataY-1)] for j in range(nb_col))
    plt.figure()
    #abscisse = np.linspace(ini_index,end_index,lenX)
    plt.plot(intensity)
    plt.show()    
    return intensity

""" 
First function : we need to find the lists of wavelengths and intensities corresponding to each spike of the spectrum in order to go further.
"""

def find_spikes_data(lambdas,intensities):
    
    spikes_data = []
    
    intensities_maxima = []
    lambdas_maxima = []
    maxima_indices = []
    
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] and intensities[i] >= 0  and intensities[i] >= threshold) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
            maxima_indices.append(i)
    # Now we need to find the spike around each maximum
    #print(maxima_indices)
    for j in range(len(lambdas_maxima)):
        
        local_spike_data = []
        
        # left minimum research
        index = maxima_indices[j]
        while (index>=0) and (intensities[index] > intensities[index-1]):
            local_spike_data.append([lambdas[index],intensities[index]])
            index -= 1
        local_spike_data.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        # right minimum research
        index = maxima_indices[j]
        while (index<len(lambdas)-1) and (intensities[index] > intensities[index+1]):
            local_spike_data.append([lambdas[index],intensities[index]])
            index += 1
        local_spike_data.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        local_spike_data.sort() # We sort the list according to the lambdas order
        local_spike_lambdas = [ local_spike_data[i][0] for i in range(len(local_spike_data)) ]
        local_spike_intensities = [ local_spike_data[i][1] for i in range(len(local_spike_data)) ]
        spikes_data.append([local_spike_lambdas,local_spike_intensities])
    
    return(spikes_data)
    

"""
For each spike we have found, we need to localize more precisly its position to know its exact wavelengths. That's why we are gonna fit a gaussian on each spike thanks to another algorithm using lmfit and the least squares method. The centre of the gaussian will be the precise wavelength of the spike, and we can also have an access to informations like its width, all those data helping us to discriminate among the spikes. We will for example be able to filter the spikes which don't have a right fitting or which have a too large width, because they are not physically realistic enough.
"""
"""
This function finds the gaussian function which is the closest from our data. It's a fitting algorithm which uses a least squares method. It returns the centre and the width of the fitted gaussian.
"""    
def fit_the_spike(lambdas,data):
    
    y = np.copy(data)
    maxY = np.max(y)
    minY = np.min(y)
    X = lambdas
    #plt.figure(1)
    #plt.plot(X, y, color = 'red')
    
    def gaussian(x,cen,amp,wid):
        return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
    
    # printing the initial gaussian (before the optimization of the fit)
    naive_center = float(np.sum(X*y))/np.sum(y)
    naive_width = np.sqrt(abs((np.sum((X-naive_center)**2*y)/np.sum(y))))
    naive_ampl = maxY-minY
    naive_gaussian = [ gaussian(x,naive_center,naive_ampl,naive_width) for x in X ]
    # plt.plot(lambdas,naive_gaussian,color='green')
    # plt.axvline(naive_center,color='green')
    plt.show()
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
        best_params_gaussian = [ gaussian(x,best_cen,best_amp,best_wid) for x in lambdas ]
        plt.figure(100)
        #plt.plot(lambdas, best_params_gaussian, 'b', color='purple')
        #plt.axvline(best_cen, color = 'red')
        computed_centre = float(np.sum(X*best_gaussian_fit))/np.sum(best_gaussian_fit) 
        plt.plot(lambdas, best_gaussian_fit, 'b--' ,color='green')
        #plt.axvline(computed_centre, color='green')
        plt.show()
        lambda_centre = best_cen
        lambda_width = best_wid
        # we need the report data to improve our understanding of the results
        report = result.fit_report()
        
        
    except : 
        report = "Computation failed for this spike : default data = naive fit"
        print(report)
        pass
        
    return(lambda_centre,lambda_width)
   

    
"""
main script that open and compute all the lists needed for the algorithm
"""
def find_ThAr_slits(test):
    path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt"
    global ini_index       # first index of the window of ThAr to be processed
    global end_index       # last index of the window of ThAr to be processed
    global lane             # considered lane
    global nbr_lanes        # Number of lanes per order
    global lenDataX         # Dimension of the CCD (length)
    global lenDataY         # Dimension of the CCD (width)
    global ThAr_data        # The ThAr data
    global threshold        # detection threshold 
    
    # Import of data from data file
    dic = collections.OrderedDict()
    with open(path, 'r') as file:
        content = file.readlines()
    content = [line.strip() for line in content]
    for line in content:
        param = line.split(" : ")
        if len(param) == 2:
            nom_param = param[0]
            value_param = param[1]
            dic[nom_param] = value_param
                 
    envelope_data_file  = dic["Lane envelope file"]
    thickness_data_file = dic["Lane thickness file"]
    ini_index          = int(dic["initial index"])
    end_index          = int(dic["final index"])
    order               = int(dic["order"])
    ThAr_file           = dic["ThAr fts file"]
    test_file            =dic["test file"]
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    
    envelope_data  = cPickle.load(open(envelope_data_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    j = 2*order + lane - 1
    left_lane = envelope_data[:, j]                     # Upper envelope of the current lane
    thickness = rint(thickness_data[j])                 # Thickness of the current lane
    
    
    if test == 'ThAr':
        image_file = pyfits.open(ThAr_file)
        ThAr_data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()
    elif test == 'test':
        ThAr_data = cPickle.load(open(test_file, 'r'))
        
    (lenDataX, lenDataY) = ThAr_data.shape
    
    intensities = sum_on_y_axis(left_lane, thickness, ini_index, end_index)
    int_min = np.min(intensities)
    print(int_min)
    norm_int = np.zeros(end_index-ini_index)
    for i in range(end_index-ini_index):
        norm_int[i]=intensities[i]-int_min
    threshold = np.mean(norm_int)  
    print(threshold)
    lambdas = np.linspace(ini_index,end_index,len(intensities), dtype = 'int')
    
    
    
    plt.figure(100)
    plt.plot(lambdas, norm_int, color = 'black')
    plt.show()
    #print(lambdas)
    #print(intensities)
    
    slits = find_spikes_data(lambdas, norm_int)
    tab_pos = []
    for s in slits:
        (x_pos, _) = fit_the_spike(s[0], s[1])
        y_mean = np.mean(s[1])
        tab_pos.append([x_pos, y_mean])
        
    if test == 'ThAr':
        cPickle.dump(tab_pos, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_slits_position_old', 'w'))
    elif test == 'test':
        cPickle.dump(tab_pos, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\test_slits_position_old', 'w'))
    
    # to plot the pos  
    tab_x = []
    tab_y = []  
    for a in tab_pos:
        tab_x.append(a[0])
        tab_y.append(a[1])
    print(tab_x)
    print(tab_y)
    
    

find_ThAr_slits('ThAr')
        
    
    












