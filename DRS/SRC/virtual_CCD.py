# -*- coding: utf-8 -*-

import cPickle
import utils
import numpy as np

"""
    Create a virtual image of a lane using the spectrum and the A matrix (AX = Y)
"""

import matplotlib.pyplot as plt

def rint(x):
    return int(round(x))


def virtual_CCD():
    cfgFile = utils.CfgFile("../DATA/Amatrix_DATA.txt")

    thickness_data_file = cfgFile.get_param("Lane thickness file")
    lanes_per_order = int(cfgFile.get_param("Lanes per order"))
    spectrum_file = cfgFile.get_param("Spectrum")
    A_matrix_file = cfgFile.get_param("A matrix file")

    print(A_matrix_file)

    thickness_data = cPickle.load(open(thickness_data_file, 'rb'))
    spectrum = cPickle.load(open(spectrum_file, 'rb'))

    (A_matrix, art, x_min, x_max) = utils.read_A_matrix(A_matrix_file)
    [order, lane] = np.array(utils.header_A_matrix(A_matrix_file).astype(int)) + [-21, 0]
    pix = 1. / art

    # First and last wavelength in arturos
    x_inf = int(rint(x_min) * pix) + 1
    x_sup = int(rint(x_max) * pix) - 2
    thickness = rint(thickness_data[lanes_per_order * order + lane - 1])
    A_matrix = A_matrix[:, x_inf:x_sup]

    spectrum = spectrum[x_inf:x_sup]
   
    img = A_matrix.dot(spectrum)
    lenImg = len(img) / thickness
    img_tilde = np.matrix(np.zeros((lenImg, thickness)))
    for i in range(lenImg):
        img_tilde[i, :] = img[i * thickness: (i + 1) * thickness].T

    return img_tilde
