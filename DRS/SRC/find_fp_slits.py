# -*- coding: utf-8 -*-

"""
    This module allows to find the centroids of the Fabry-Perot slits. 
    It uses a Data file which gatherx all the paths toward the needed data.
    
    Those are :
        - the fts file of a Fabry-Perot (already pre-processed)
        - the envelope and the thickness of the lanes
        - the number of lanes per order
"""

import utils
import lmfit
import pyfits
import cPickle
import peakutils
import numpy as np
from scipy.optimize import leastsq
from multiprocessing import Pool


# Method used to get the closest integer of a real number
def rint(x): return int(round(x))


def approx_fp_centres(left_lane, thickness, ini_index, end_index):
    """
        Returns the estimated position (on X axis) of the Fabry-Perot slits of
        a given lane on a given order. The slits are situated with the maxima of intensity in the data.
        First, those maximas are roughly situated to find an estimation of the gap between them 
        and then the maximas are situated with a precision of 1 pixel.

        : left_lane : the upper envelope of the lane
        : thickness : the thickness of the lane
        : ini_index : the first index to be considered
        : end_index : the last index to be considered

        : return : a list of the indices of the slits (roughly situated)
    """
    lenX = end_index - ini_index
    rowIntensity = np.zeros(lenX)
    nb_col = rint(thickness)

    # For each row, the values of the corresponding column are added
    for i in range(ini_index, end_index):
        rowIntensity[i-ini_index] = sum(fp_data[i][min(rint(left_lane[i]+j), lenDataY-1)] for j in range(nb_col))

    # For a better detection of the slits, the background is set to 0 using windows of 15 pixels wide
    local_min_win = 15
    for i in range(0, lenX, local_min_win):
        local_min = np.min(rowIntensity[i:min(i+local_min_win, lenX)])
        rowIntensity[i:min(i+local_min_win, lenX-1)] -= local_min

    rowIntensity2 = rowIntensity[50:-50]

    # The local maxima are the estimated position of the split's centroids
    threshold = 0.08
    indexes = peakutils.indexes(rowIntensity2, thres=threshold, min_dist=10)

    # In case not enough peaks were detected, the threshold is reduced
    if (len(indexes) < 50):
        threshold = 0
        indexes = peakutils.indexes(rowIntensity2, thres=threshold, min_dist=10)
    
    # The first estimation give an approximate gap between each slit
    # This mean gap is used to detect more precisely the maxima
    gap = rint(np.mean(indexes[1:]-indexes[0:-1])/2) # gap/2 gives the radius of the search window
    indexes = peakutils.indexes(rowIntensity2, thres=threshold, min_dist=gap)

    gap = rint(np.mean(indexes[1:]-indexes[0:-1]))
    init_gap = gap

    # We are looking for the first local minima in a window where we are sure there is one slit (0.7*gap)
    # This way, the first slit may be avoided if it is not clear enough
    first_min = np.argmin(rowIntensity[0:rint(0.7 * gap)])
    # The first slit is then easily found in the window
    first_index = np.argmax(rowIntensity[first_min:first_min+gap]) + first_min

    # ini_window is the variable first index of the search window
    ini_window = max(0, first_index - rint(gap/2))
    indexes = []

    while (ini_window < lenX - gap):
        # The maxima are precisely situated in each search window
        index = np.argmax(
            rowIntensity[ini_window:min(lenX-1, ini_window + gap)])
        indexes.append(ini_window + index + ini_index)

        # As the gap between two slits slowly changes, it is updated regularly
        if (len(indexes) > 1):
            gap = rint(0.5*(max(indexes[-1] - indexes[-2], gap) + init_gap))

        ini_window = indexes[-1] - ini_index + rint(gap * 0.4)

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


def fp_centre_local(left_lane, thickness, x):
    """
        Returns the precise position of one slit of the F-P.

        : left_lane : the upper envelope of the lane
        : thickness : the thickness of the lane
        : x : the approxiative position of one slit (X axis)

        : return : the position of the slit (precise on X, approximate on Y)
    """
    search_window = 8

    max_xsearch = max(0, x-search_window)

    # The data in which we look for the slit
    srch_data = [sum(fp_data[k, rint(left_lane[k]): rint(left_lane[k] + thickness)])
                 for k in range(max_xsearch, min(lenDataX, x+search_window+1))]

    # A Gaussian is then fitted to the flattened data to find the centre
    gaussian_slit = find_gaussian_slit(srch_data)
    if (0 <= gaussian_slit < len(srch_data)):
        x_centr = min(gaussian_slit + max_xsearch, lenDataX-1)

        # We then proceed to find the centre on the Y-axis
        x_rint = rint(x_centr)
        x_left = rint(left_lane[x_rint])
        x_right = rint(left_lane[x_rint] + thickness)

        # The centre column of the split is selected
        # dtype is given to have a signed array
        y_search = np.array(fp_data[x_rint, x_left:1+x_right], dtype='int')

        # The second local maxima is selected as the estimated centre
        # Indeed, there should be 3 maxima (1 for each slice of the fibre)
        min_norm = float(np.mean(y_search)) / (1.5*np.max(y_search))
        indexes = peakutils.indexes(np.array(y_search), thres=min_norm)

        y = indexes[min(1, abs(len(indexes) - 1))] if len(indexes) else (thickness)/2

    return (x_centr, y)


def fp_centre_global(left_lane, thickness, ini_index, end_index):
    """
        Returns the precise position (on X axis) and the estimated one (Y axis)
        of all the Fabry-Perot split of a given lane.

        : left_lane : the upper envelope of the lane
        : thickness : the thickness of the lane
        : ini_index : the first index to be considered
        : end_index : the last index to be considered

        : return : the list of the position (X axis) and the estimated relative position on Y axis
    """
    # To find the precise position of the centroids, we first estimate it
    fp_centre = approx_fp_centres(left_lane, thickness, ini_index, end_index)

    x_tab, y_tab = [], []

    for x in fp_centre:
        (x, y) = fp_centre_local(left_lane, thickness, x)
        x_tab.append(x)
        y_tab.append(y)

    # As the distance between the centroids and the lane (left or right) must
    # not change, the mean value of this distance is chosen as the true value.
    y_mean = np.mean(y_tab)

    return x_tab, y_mean


def fit_deviation(X_pos):
    """
        Fit the position with a given polynomial to improve the results.

        : X_pos : 1D array of the position on X axis for one lane

        : return : fitted 1D array of the position on X axis
    """
    lenD = len(X_pos) - 1
    X = np.arange(lenD).astype(float)

    # Distance between two successive slits
    deviation = [X_pos[i+1] - X_pos[i] for i in range(lenD)]

    # Coefficient of the polynomial used to fit
    a, b, c, d, e = 25, 0.05, 0, -0.005, 0
    
    # Function to apply the polynomial with given coefficients
    calc_func       = lambda x:  x[0] + x[1]*X + x[2]*X**2 + x[3]*X**3 / (1+x[4]*X)
    # Function to be optimized
    optimize_func   = lambda x:  calc_func(x) - deviation

    # First estimation of the coefficients
    moments = leastsq(optimize_func, (a, b, c, d, e))[0]
    data_fit1 = calc_func(moments)

    # Gap between the real data and the fitted one
    fit_gap = np.array(data_fit1) - np.array(deviation)
    std = np.std(fit_gap)

    # The top 80% of the fitted data must not have a stp > 0.05
    # If not, we fit a second time the data
    std_test = np.std(np.sort(fit_gap)[:int(0.8*len(fit_gap))])

    if std_test > 0.05:
        # Proportion of the fitted data we trust (center of the lane)
        # It's bigger if the fitted-data was not so bad
        coeff = 0.7 if std_test < 2 else 0.5
        len_ok = int(lenD * coeff)

        # We look for for the best part of the lane with the given size
        init, end = 0, len_ok
        std_min = std
        init_min = int((lenD - len_ok) * coeff)

        while (end < lenD):
            std_temp = np.std(fit_gap[init:end])
            if std_temp < std_min:
                init_min = init
                std_min = std_temp
            init += 1
            end += 1

        init = init_min
        end = init + int(lenD * coeff)

        # We then re-fit only the trusted part of the lane
        X = np.arange(len_ok).astype(float) + init
        deviation = deviation[init:init+len_ok]
        moments = leastsq(optimize_func, (a, b, c, d, e))[0]

        # The new fit is extended to the whole lane
        X = np.arange(lenD).astype(float)
        fitted_data = calc_func(moments)

        # The positions that were not ok are then changed to the fitted ones
        for i in range(init, -1, -1):           # Before the trusted part
            X_pos[i] = X_pos[i+1] - fitted_data[i]

        for i in range(init+len_ok, lenD+1):    # After the trusted part
            X_pos[i] = X_pos[i-1] + fitted_data[i-1]

    return X_pos


def cut_border_edge(left_lane, thickness):
    """
        Determine the first and the last index of the lane.
        It is important for the last lanes (orders 60 & 61), when the CCD's edge cut the lanes.

        : left_lane : the lupper enveloppe of the lane
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


def launch(j):
    """
        This function launch the appropriate functions to find the positions of the slits
        for one lane. 

        : j : the cosnidered lane

        : return : an array with all the positions (x,y) for the lane
    """

    print("Cherche FP dans Ordre n°{} - Voie n°{}".format(j/nbr_lanes, j % nbr_lanes+1))

    left_lane = envelope_data[:, j]                     # Upper envelope of the current lane
    thickness = rint(thickness_data[j])                 # Thickness of the current lane

    (ini, end) = cut_border_edge(left_lane, thickness)  # First and last index of the lane

    # First research of the slits, with potential errors on X (x_pos is an array, y is a real)
    # Y coordinates are relative to the envelope (not absolute coordinates)
    x_pos, y = fp_centre_global(left_lane, thickness, ini, end)

    # Fit of the position on X axis to increase precision
    x_tab = fit_deviation(x_pos)

    # Y coordinates are set as absolute
    tab_pos = [[x, y + left_lane[rint(x)]] for x in x_tab if (ini <= x <= end)]

    return tab_pos


def find_fp_slits():
    """
        Import all the data needed to find the F-P slits and eventually save the positions with correct format.

        : return : save the array with all the positions
    """
    global envelope_data    # Upper envelope of each lane
    global thickness_data   # Thickness of the lanes
    global nbr_lanes        # Number of lanes per order
    global lenDataX         # Dimension of the CCD (length)
    global lenDataY         # Dimension of the CCD (width)
    global fp_data          # The F-P data

    # Import of data from data file
    cfgFile = utils.CfgFile("../DATA/find_fp_slits_DATA.txt")

    envelope_data_file  = cfgFile.get_param("Lane envelope file")
    thickness_data_file =  cfgFile.get_param("Lane thickness file")

    fp_file     = cfgFile.get_param("Fabry Perot fts file")
    nbr_lanes   = int(cfgFile.get_param("Lanes per order"))

    envelope_data  = cPickle.load(open(envelope_data_file, 'rb'))
    thickness_data = cPickle.load(open(thickness_data_file, 'rb'))

    image_file = pyfits.open(fp_file)
    fp_data = image_file[0].data.astype(np.float32) # Data of the Fabry-Perot fts file
    image_file.close()

    (lenDataX, lenDataY) = fp_data.shape

    init_order, final_order = 0, 39    # We assume there are 40 orders

    # Multiprocessing : each lane is independant from the others
    p = Pool(10)
    pos_global = p.map(launch, range(nbr_lanes*init_order, nbr_lanes*final_order + nbr_lanes, 1))
    p.terminate()
    p.join()


    cPickle.dump(np.array(pos_global), open('../TEMP/fp_slit_position.p', 'wb'))
    print('FP lines info written in {0}'.format('TEMP/fp_slit_position.p'))

find_fp_slits()
