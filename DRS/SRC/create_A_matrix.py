# -*- coding: utf-8 -*-

"""
    This module creates the A matrix of the current lane. It uses a Data file
    which gather all the paths toward the needed data.
    
    Those are :
        - the fts file of a Fabry-Perot (already pre-processed)
        - the positions of the F-P slits
        - the envelope and the thickness of the lanes
        - the number of lanes per order
        - the considered lane (order + lane)
        - the version of the matrix (pre or post cosmic-process)
    
    The matrix is sparse (CSR format), and is saved with a header written
    on its upper-right hand corner.
"""

import pyfits
import cPickle
import utils
import datetime
import numpy as np
from scipy import sparse
from multiprocessing import Pool

import matplotlib.pyplot as plt


# Method used to get the closest integer of a real number
def rint(x): return int(round(x))


def fill_A_fp(ref_data):
    """
        Creates a matrix only filled at Fabry-Perot wavelength

        : ref_data : The Fabry-Perot data

        : return : Sparse A matrix
    """
    # Reminder : in the A matrix, the columns are the wavelength and the rows are the different pixels of the CCD
    row, col, data = [], [], []
    search_window = 0 # Width of the window containing the slit 

    rng_thickness = range(thickness)

    # # Artificial PSF
    # r=search_window
    # PSF=np.zeros((2*r+1,thickness))
    # X=np.arange(-10*r-5,(10*r)+5)
    # x0=(np.arange(0,thickness)-thickness/2.)*np.sin(11*np.pi/180.)


    for [x_pos, _] in tab_centre:   # We use each F-P slit of the lane, one after another

        if (ini_index <= x_pos <= end_index):   # Check that the slit really is on the CCD
            # Minimum and maximum indices of the search window
            min_x = rint(max(ini_index, x_pos-search_window))
            max_x = rint(min(end_index + 1, x_pos+search_window+1))

            # # Creation of the artificial leaning slit
            # residu = x_pos*pix - rint(x_pos*pix)
            # for i in rng_thickness:
            #     perfh=np.exp(-(X-x0[i]-10*residu)**2/(2*(2.**2)))
            #     PSF[:,i]=np.sum(perfh.reshape((2*r+1,10)),1)
            # PSF = PSF * 1000

            wavelength = rint(x_pos * pix)      # Considered wavelength in arturos
            delta_x_pos = x_pos - rint(x_pos)   # Error due to round number

            for x in range(min_x, max_x):

                for y in rng_thickness: 

                    row.append(rint((x + delta_x_pos) * thickness + y))
                    col.append(wavelength)

                    # Real data
                    data.append(float(ref_data[rint(x + delta_x_pos), rint(y + left_lane(x))]))
                    
                    # # Artificial data
                    # data.append(float(PSF[rint(x - x_pos + search_window), rint(y)]))

    # Creation of the sparse A matrix filled at F-P wavelengths
    A_shape = (lenX*thickness, rint(lenX*pix))
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float')

    return A


def interpolate(i):
    """
        Interpolate the A matrix between two F-P slits.

        : i : the index of the considered left slit in the slit-position tab

        : return :  the information for each interpolated pixel (row, column and intensity) (list format)
                    and the first and last wavelength (in pixel, not arturos)
    """
    row, col, data = [], [], []

    x0, x1 = tab_centre[i][0], tab_centre[i+1][0]   # We consider two successive slits of the F-P.

    if (x0 < end_index and x1 < end_index+1):
        arx0, arx1 = rint(x0 * pix), rint(x1 * pix)                 # Rounded wavelength of the slits in arturos
        dy = -float(left_lane(x1) - left_lane(x0)) / float(x1-x0)   # Gap on the Y axis between the slits because of the curved lane

        # Those values are used after. They are created here to reduce the complexity.
        search_window = 0
        f_diff = float(x1 - x0)
        rng_thickness = range(thickness)
        frac_dy = (1. / (1 + 2 * (dy ** 2)))

        # We successively consider each wavelength (hence each column) between the two slits.
        # for j in range(arx0 + 1, arx1):     # Case 1 : real value of the F-P are kept (discontinuity)
        for j in range(arx0, arx1):         # Case 2 : the F-P wavelenfgth are interpolated (continuous)

            # Minimum and maximum indices of the search window
            min_x = rint(max(ini_index * pix, j - search_window))
            max_x = rint(min((end_index + 1) * pix, j + search_window + 1))

            # Coefficient for the interpolation (equals to 0 at lambda0 and 1 at lambda1)
            coeff = art * float(j-arx0) / f_diff
            test_tab = []

            for x in range(min_x, max_x):   # Cross of the search window on the X axis
                # row_v = rint(x*art*thickness) 

                for y in rng_thickness:     # Cross of the search window on the Y axis
                    # row_value = row_v + y   # Considered pixel
                    row_value = rint(x*art*thickness + y - left_lane(x) + rint(left_lane(i)))   # Considered pixel

                    if (row_value not in test_tab):
                        pos_0 = rint( ((x-j) * art + x0) * thickness ) + y # Pixel equivalent positon in the F-P slit column (left one)
                        pos_1 = rint( ((x-j) * art + x1) * thickness ) + y # Pixel equivalent positon in the F-P slit column (right one)

                        p0 = A_fp[pos_0, arx0]  # Pixel participation intensity (left F-P column)
                        p1 = A_fp[pos_1, arx1]  # Pixel participation intensity (rigth F-P column)

                        # Quadratic interpolation
                        p1_sup1 = A_fp[pos_1 + 1, arx1] # Intensity participation of the upper pixel of the rigth F-P column
                        p1_inf1 = A_fp[pos_1 - 1, arx1] # Intensity participation of the lower pixel of the rigth F-P column

                        p1_sup2 = A_fp[pos_1 + 2, arx1]
                        p1_inf2 = A_fp[pos_1 - 2, arx1]

                        p1_tilde_i     = frac_dy * (p1      + dy * (p1_inf1 - p1_sup1 ))    # Expected intensity of the pixel
                        p1_tilde_i_inf = frac_dy * (p1_inf1 + dy * (p1_inf2 - p1))          # Expected intensity of the lower pixel
                        p1_tilde_i_sup = frac_dy * (p1_sup1 + dy * (p1 - p1_sup2))          # Expected intensity of the upper pixel

                        a = -dy * (p1_tilde_i_sup - p1_tilde_i_inf) / f_diff
                        p = float(p0 + coeff * float(abs(a + p1_tilde_i) - p0)) # Interpolated intensity

                        # Column, row, and intensity of the current pixel added
                        col.append(j)
                        row.append(row_value)
                        data.append(p)

                        test_tab.append(row_value)

    return [col, row, data, x0, x1]


def fill_A_matrix(ref_data, date):
    """
        Creates the full A matrix

        : ref_data : The Fabry-Perot data
        : date : The date of the F-P file

        : return : The full interpolated A matrix
    """
    global art      # Number of pixel per arturo
    global pix      # Number of arturos per pixel
    global A_fp     # A matrice only filled for F-P wavelengths

    pix = 1.       # 1.5 art/pix <=> 1.73km/s / pix 
    art = 1./pix

    A_fp = fill_A_fp(ref_data).tocsr()  # Creation of the A matrix only filled with F-P

    # Interpolation between the slits with multiprocessing
    p = Pool(10)
    DATA = p.map(interpolate, range(0, len(tab_centre) - 1))
    p.terminate()
    p.join()

    # Retrieval of the data (concatenate all the lists)
    zip_data = zip(*DATA)
    col = sum(list(zip_data[0]), [])
    row = sum(list(zip_data[1]), [])
    data = sum(list(zip_data[2]), [])
    x_min = list(zip_data[3])
    x_max = list(zip_data[4])

    x_min = min(end_index, np.min(x_min))   # First wavelength of the whole interpolaton (in pixel, not arturos)
    x_max = max(ini_index, np.max(x_max))   # Last wavelength of the whole interpolaton (in pixel, not arturos)

    rp_lenX = rint(lenX * pix) # Round number of arturos in one lane

    # Adding the header in the upper-rigth hand corner of the matrix (first line)
    # The header has : the date of the F-P, the order and the lane, the number of pixel per aturo,
    # and the first and the last inteprolated wavelengths (in pixels, not arturos).
    row += [0] * 6
    col += range(rp_lenX - 6, rp_lenX)
    data += [date.strftime("%Y%m%d%H%M"), order + 21, lane, art, x_min, x_max]

    # Creation of the full interpolated matrix with the header, in CRS format
    A_shape = (rint(lenX*thickness), rp_lenX)
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float')
    A_full = A.tocsr()

    return A_full


def cut_border_edge(lenDataY):
    """
        Determine the first and the last index of the lane.
        It is important for the last lanes (orders 60 & 61), when the CCD's edge cut the lanes.

        : lenDataY : the height of the CCD

        : return : creates two globals variables of the first and last index    
    """
    global ini_index
    global end_index 

    ini_index = 0
    end_index = lenX - 1
    data_started = False    # False while data not encountered along the CCD

    for i in range(lenX):

        if (not left_lane(i)) or (left_lane(i) + thickness >= lenDataY):
            # Current index is at the left of data
            if not data_started and i >= ini_index:
                ini_index = min(i+1, end_index)

            # Current index is at te right of data
            elif data_started and i >= ini_index:
                end_index = max(ini_index, i-1)
                break # There is no need to continue as data index has already been exceeded

        else: # data encountered
            data_started = True


def create_A_matrix():
    """
        Import all the data needed to create the A matrix and eventually save the matrix with correct syntax.
        
        : return : save A matix in CSR format
    """
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global tab_centre   # positions of the Fabry-Perot slits
    global order        # Considered order
    global lane         # Considered lane within the order
    global lenX         # Dimension of the CCD (length)

    # Import of data from data file
    cfgFile = utils.CfgFile("../DATA/Amatrix_DATA.txt") # Path to the data file

    thickness_data_file = cfgFile.get_param("Lane thickness file")
    polys_data_file     = cfgFile.get_param("Polys envelope file")
    fp_position_file    = cfgFile.get_param("FP slits file")
    fp_file             = cfgFile.get_param("Fabry Perot fts file") 

    order           = int(float(cfgFile.get_param("Order n°"))) # Orders are counted from 0 in data file (0 <-> 21)
    lane            = int(float(cfgFile.get_param("Lane n°")))  # Lanes are counted from 1
    nbr_lanes       = int(cfgFile.get_param("Lanes per order"))
    evolution_level = cfgFile.get_param("Evolution level")

    polys_data     = cPickle.load(open(polys_data_file, 'rb'))
    pos_global     = cPickle.load(open(fp_position_file, 'rb'))
    thickness_data = cPickle.load(open(thickness_data_file, 'rb'))

    image_file = pyfits.open(fp_file)
    fp_data = image_file[0].data.astype(np.float32)     # Data of the Fabry-Perot fts file
    image_file.close()

    thickness_float = thickness_data[nbr_lanes * order + lane - 1]
    thickness = rint(thickness_float)   # The thickness of the considered lane

    left_lane = polys_data[order]       # A polynomial giving the upper enveloppe of the lane
    left_lane[0] += thickness_float * (lane - 1)

    tab_centre = pos_global[nbr_lanes * order + lane - 1]   # The coordinates of the centres of the Fabry-Perot slits

    (lenX, lenY) = fp_data.shape
    cut_border_edge(lenY)               # Locating the edge of the lane

    date_txt = fp_file[-23:-10]         # Extraction of the date/time frome the name of the F-P file
    date = datetime.datetime.strptime(date_txt, "%Y%m%d_%H%M") 

    A_matrix = fill_A_matrix(fp_data, date)     # Creation of the full A matrix
    A_matrix_name = '../TEMP/Amatrix_' + 'v' + evolution_level + '_' + 'OR' + str(order) + '_LA' + str(lane) + '.npz'
    sparse.save_npz(A_matrix_name, A_matrix)    # A matrix saved with sparse format

    cfgFile.modif_param('A matrix file', A_matrix_name)     # The path to the matrix is updated in the data file


create_A_matrix()
