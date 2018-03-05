# -*- coding: utf-8 -*-

"""
    This module creates a spectrum of the data given using the A matrix. 
    It uses a Data file which gathers all the paths toward the needed data.
    
    Those are :
        - the A matrix of the lane
        - the fts file of which we want the spectrum
        - the envelope and the thickness of the lanes
        - the number of lanes per order
        - the considered lane (this is directly read in the matrix header)
    
    The generated spectrum is a list saved with the .p format of cPickle.
"""

import utils
import pyfits
import cPickle
import numpy as np

import matplotlib.pyplot as plt


# Method used to get the closest integer of a real number
def rint(x): return int(round(x))


def generate_spectrum(test_data, A_matrix_name):
    """
        Generate the spectrum thanks to the method of Normal Equations.A_matrix

        : test_data : the data of which he spectrum has to be generated
        : A_matrix_name : the path to the A matrix to retrieve its data

        : return : the generated spectrum as a 1D array
    """
    # Retrieval of the data from the matrix
    (A_matrix, art, x_min, x_max) = utils.read_A_matrix(A_matrix_name)

    pix = 1. / art

    # First and last wavelength in arturos
    x_inf = int(rint(x_min) * pix) + 1
    x_sup = int(rint(x_max) * pix) - 2
    # As the interpolation was not done up to the edge, we have to cut the matrix
    # Indeed, the C matrix (see below) has to be symetrical and positive
    A_matrix = A_matrix[:, x_inf:x_sup]

    # The CCD is set as a vector (from 2D to a 1D matrix)
    CCD = []
    for i in range(lenX):
        for j in range(thickness):
            CCD.append(test_data[i, rint(left_lane(i)+j)])
    CCD = np.matrix(CCD).T

    # We use the method of Normal Equations to retrieve the spectrum of the lane
    C = np.dot(A_matrix.T, A_matrix)            # C = At * A
    C = C + np.identity(C.shape[0]) * 0.000001  # We get rid of the potential null values on the diagonal

    G = np.linalg.cholesky(C)   # Cholesky factorization
    d = (A_matrix.T).dot(CCD)
    y = np.linalg.solve(G, d)
    X = np.linalg.solve(G.T, y)

    spectrum = ([-1] * x_inf + list(np.array(X.T.tolist()[0])) + [-1] * (rint(lenX * pix) - x_sup)) 

    # The spectrum is returned as an 1D array.
    return spectrum


def create_spectrum():
    """
        Import all the data needed to create the spectrum and eventually save it with correct syntax.
        
        : return : save spectrum in .p format
    """
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global lenX         # Dimension of the CCD (length)

    # Import of data from data file
    cfgFile = utils.CfgFile("../DATA/Amatrix_DATA.txt")

    thickness_data_file = cfgFile.get_param("Lane thickness file")
    polys_data_file     = cfgFile.get_param("Polys envelope file")
    test_file           = cfgFile.get_param("Test fts file")
    A_matrix            = cfgFile.get_param("A matrix file")
    nbr_lanes           = int(cfgFile.get_param("Lanes per order"))

    polys_data      = cPickle.load(open(polys_data_file, 'rb'))
    thickness_data  = cPickle.load(open(thickness_data_file, 'rb'))

    image_file = pyfits.open(test_file)
    test_data  = image_file[0].data.astype(np.float32)
    image_file.close()

    # The considered lane and order are directly read in the A matrix header, to avoid any possible problem with the matrix name
    [order, lane] = np.array(utils.header_A_matrix(A_matrix).astype(int)) + [-21, 0]

    thickness_float = thickness_data[nbr_lanes * order + lane - 1]
    thickness = rint(thickness_float)   # The thickness of the considered lane

    left_lane = polys_data[order]       # A polynomial giving the upper enveloppe of the lane
    left_lane[0] += thickness_float * (lane - 1)

    (lenX, _) = test_data.shape

    # creation of the spectrum of the lane
    spectrum = generate_spectrum(test_data, A_matrix)
    # Saving the spectrum
    Pickle_name = '../TEMP/SP_' + test_file.split('/')[-1][:-4] + '_OR' + str(order) + '_LA' + str(lane) + '.p'
    cPickle.dump(spectrum, open(Pickle_name, 'wb'))

create_spectrum()