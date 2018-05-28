# -*- coding: utf-8 -*-

"""
    This module allows to find the centroids of the ThAr slits. 
    It uses a Data file which gatherx all the paths toward the needed data.
    
    Those are :
        - the fts file of a ThAr (already pre-processed)
        - the envelope and the thickness of the lanes
        - the number of lanes per order
"""
import lmfit
import pyfits
import cPickle
import peakutils
import numpy as np
from scipy.optimize import leastsq
from multiprocessing import Pool
import collections
import matplotlib.pyplot as plt


def rint(x): return int(round(x))

def sum_on_y_axis(left_lane, thickness, ini_index, end_index):
    lenX = end_index - ini_index
    intensity = np.zeros(lenX)
    nb_col = rint(thickness)
    for i in range(ini_index, end_index):
        intensity[i-ini_index] = sum(ThAr_data[i][min(rint(left_lane[i]+j), lenDataY-1)] for j in range(nb_col))
    plt.figure(104)
    abscisse = np.linspace(ini_index,end_index,lenX,dtype ='int')
    plt.plot(abscisse, intensity)
    plt.show()    
    return intensity
    
    
def cut_border_edge(left_lane, thickness):
    """
        Determine the first and the last index of the lane.
        It is important for the last lanes (orders 60 & 61), when the CCD's edge cut the lanes.

        : left_lane : the upper enveloppe of the lane
        : thickness : the thickness of the considered lane

        : return : the first and the last index of the lane    
    """
    init = 0
    end = lenDataX - 1
    test = False        # False while data not encountered along the CCD

    for i in range(lenDataX):

        if (not left_lane[i]) or (left_lane[i]+thickness >= lenDataY):
            # Current index is at the left of the data
            if not test and i >= init:
                init = min(i+1, end)
            # Current index is at te right of the data
            elif test and i >= init:
                end = max(init, i-1)
                break# There is no need to continue as data index has already been exceeded

        else: # data encountered
             test = True

    return (init, end) 
    
def approx_center(left_lane, thickness, ini_index, end_index):
    
    intensity = sum_on_y_axis(left_lane,thickness, ini_index, end_index)
    
    indexes = peakutils.indexes(intensity, thres = threshold, min_dist = 1)
    
    return indexes
    
def find_gaussian_slit(data):
    """
        Returns the centre of a gaussian fit applied to the given data.

        : data : a 1D array where there is a unique gaussian

        : return : the centre of the gaussian
    """
    data_copy = np.copy(data) - np.min(data)
    X = np.arange(len(data_copy))
    
    # The estimated centroid position is derived with the moments method
    centre = float(np.sum(X*data_copy))/np.sum(data_copy)
    height = np.max(data_copy)
    width = np.sqrt(np.abs(np.sum((X-centre)**2*data_copy)/np.sum(data_copy)))

    def gaussian(x, cen, amp, wid):
        return amp * np.exp(-(x-cen)**2 / wid)

    try:
        # Using a fit with the estimated moments increase the precision
        gmod = lmfit.Model(gaussian)
        params = gmod.make_params(cen=centre, amp=height, wid=width)
        result = gmod.fit(data_copy, params, x=X)
        [cen, wid, amp] = list(result.best_values.values())

        if (0 < cen < len(data_copy)):
            centre = cen
    except:
        pass
    return centre
    
def ThAr_local_center(left_lane, thickness, x):
    

    max_xsearch = max(0, x-search_window)
    #print((rint(left_lane[x]), rint(left_lane[x] + thickness)))   
    #print((max_xsearch, min(lenDataX, x+search_window+1)))
    #print(x)
    #for k in range(max_xsearch, min(lenDataX, x+search_window+1)):
      #  print(ThAr_data[k, rint(left_lane[k]): rint(left_lane[k] + thickness)])
    
    # The data in which we look for the slit
    srch_data = []
    for k in range(max_xsearch, min(lenDataX, x+search_window+1)):
        srch_data.append(sum(ThAr_data[k, rint(left_lane[k]): rint(left_lane[k] + thickness)]))
    
    
# A Gaussian is then fitted to the flattened data to find the centre
    gaussian_slit = find_gaussian_slit(srch_data)
    if (0 <= gaussian_slit < len(srch_data)):
        x_centr = min(gaussian_slit + max_xsearch, lenDataX-1)      
        # We then proceed to find the centre on the Y-axis
        x_rint = rint(x_centr)
        x_left = rint(left_lane[x_rint])
        x_right = rint(left_lane[x_rint] + thickness)
        #print(x_rint, x_left, x_right)
        # The centre column of the split is selected
        # dtype is given to have a signed array
        y_search = np.array(ThAr_data[x_rint, x_left:1+x_right], dtype='int')
        #print(y_search)
        # The second local maxima is selected as the estimated centre
        # Indeed, there should be 3 maxima (1 for each slice of the fibre)
        min_norm = float(np.mean(y_search)) / (1.5*np.max(y_search)+0.001)
        indexes = peakutils.indexes(np.array(y_search), thres=min_norm)

        y = indexes[min(1, abs(len(indexes) - 1))] if len(indexes) else (thickness)/2
    return (x_centr, y)
                 
            
def ThAr_global_center(left_lane, thickness, ini_index, end_index):
    """
        Returns the precise position (on X axis) and the estimated one (Y axis)
        of all the ThAr split of a given lane.

        : left_lane : the upper envelope of the lane
        : thickness : the thickness of the lane
        : ini_index : the first index to be considered
        : end_index : the last index to be considered

        : return : the list of the position (X axis) and the estimated relative position on Y axis
    """
    # To find the precise position of the centroids, we first estimate it
    ThAr_center = approx_center(left_lane, thickness, ini_index, end_index)
        
    
    #print(ThAr_center)
    x_tab, y_tab = [], []
    plt.figure(105)
    for x in ThAr_center:
        (x, y) = ThAr_local_center(left_lane, thickness, x + ini_index)
        plt.axvline(x, color='red')
        x_tab.append(x)
        y_tab.append(y)
    plt.show()
    # As the distance between the centroids and the lane (left or right) must
    # not change, the mean value of this distance is chosen as the true value.
    y_mean = np.mean(y_tab)

    return x_tab, y_mean
                
def launch(ini, end, order):
    """
        This function launch the appropriate functions to find the positions of the slits
        for one lane. 

        : j : the cosnidered lane

        : return : an array with all the positions (x,y) for the lane
    """
    j = nbr_lanes*order + lane-1
    left_lane = envelope_data[:, j]                     # Upper envelope of the current lane
    thickness = rint(thickness_data[j])                 # Thickness of the current lane

    #(ini, end) = cut_border_edge(left_lane, thickness)  # First and last index of the lane

    # First research of the slits, with potential errors on X (x_pos is an array, y is a real)
    # Y coordinates are relative to the envelope (not absolute coordinates)
    x_pos, y = ThAr_global_center(left_lane, thickness, ini, end)

    # Y coordinates are set as absolute
    tab_pos = [[x, y + left_lane[rint(x)]] for x in x_pos if (ini <= x <= end)]
    #print(tab_pos)
    return tab_pos                 
                 
def find_ThAr_slits(test):
    
    path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt"
    global ini_index       # first index of the window of ThAr to be processed
    global end_index       # last index of the window of ThAr to be processed
    global envelope_data    # Upper envelope of each lane
    global thickness_data   # Thickness of the lanes
    global lane             # considered lane
    global nbr_lanes        # Number of lanes per order
    global lenDataX         # Dimension of the CCD (length)
    global lenDataY         # Dimension of the CCD (width)
    global ThAr_data        # The ThAr data
    global threshold        # detection threshold for peakutil
    global search_window
    search_window = 8
    threshold = 0.05
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
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    order               = int(dic["order"])
    ThAr_file           = dic["ThAr fts file"]
    FP_file             = dic["FP fts file"]
    test_file           = dic["test file"]
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])

    envelope_data  = cPickle.load(open(envelope_data_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    if test == 0:
        image_file = pyfits.open(ThAr_file)
        ThAr_data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()
    elif test == 1:
        ThAr_data = cPickle.load(open(test_file, 'r'))
    elif test == 2:
        image_file = pyfits.open(FP_file)
        ThAr_data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()

    (lenDataX, lenDataY) = ThAr_data.shape
    
    tab_pos = launch(ini_index, end_index, order)
    print(tab_pos)
    print(len(tab_pos))
    if test ==0:
        cPickle.dump(tab_pos, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_slits_position_old', 'w'))
    elif test == 1:
        cPickle.dump(tab_pos, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\test_slits_position_old', 'w'))
    elif test == 2:
        cPickle.dump(tab_pos, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\FP_slits_position_old', 'w'))
    #print(tab_pos)
    # to plot the pos             
    tab_x = []
    tab_y = []  
    for a in tab_pos:
        tab_x.append(a[0])
        tab_y.append(a[1])
    #print(tab_x)
    #print(tab_y)
    #plt.figure()
    #plt.plot(tab_x,tab_y,'r.')
    #plt.show()
    
    
# test = 0 : ThAr
# test = 1 : CCD random issue de CCD _creator
# test = 2 : FP
find_ThAr_slits(1)
        
        
    
                 
                 
                 
                 
                 
                 
                 
                 