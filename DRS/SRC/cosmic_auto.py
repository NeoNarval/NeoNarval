#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import pyfits
import numpy as np
from math import *
import cPickle as pickle
import utils
import virtual_CCD
import sys
import os.path

# Argument 1 : fichier output
# Argument 2 : le spectre est-il normalise


def erase_cosmic(fichier, file_envelope, order):
    """
                             This method get us rid off cosmic impacts one the CCD
                                    fichier : the fits file we use
                                    file_Y: the pickle file which contains the CCD matrix create by the matrix A
                                    file_envelope : the pickle file which contains the border of the orders
                            It creates a new fits files with the correct name.
    """
    # We get the data in the fits file
    hdulist = pyfits.open(fichier)
    image = hdulist[0].data.astype(np.float32)
    (lenX, lenY) = image.shape
    hdulist.close()
    # We get the envelope of the orders
    pkl_file = open(file_envelope, "rb")
    envelope = pickle.load(pkl_file)
    # We prepare the new image to be create
    new_image = image
    diff = []
    cfg = utils.CfgFile("../DATA/Amatrix_DATA.txt")
    lvl = cfg.get_param("Evolution level")
    fichier_fts = os.path.basename(fichier)
    for order in xrange(40):
        cfg.modif_param("Order n°", str(order))
        for lane in range(1, 3):
            # Configuration pour virtual_CCD
            cfg.modif_param("Lane n°", str(lane))
            cfg.modif_param("A matrix file",
                            "../TEMP/Amatrix_v{}_OR{}_LA{}.npz".format(lvl, order, lane))
            try:
                normalized = int(sys.argv[2])
            except:
                normalized = 0
            if normalized:
                spectrum = "../TEMP/SN_" + fichier_fts[:-4] + "_OR{}_LA{}.p".format(order, lane)
            else:
                spectrum = "../TEMP/SP_" + fichier_fts[:-4] + "_OR{}_LA{}.p".format(order, lane)

            cfg.modif_param("Spectrum",
                            spectrum)
            Y = virtual_CCD.virtual_CCD()
            j = order
            max_l = []
            max_Y = []
            # For each order, assuming the cosmic that bother us is much more intense than the true data, we search the maximum in the order for each line
            for i in range(lenX):
                if (envelope[i, j] == 0 or envelope[i, j + 1] == 0):
                    max_l.append(0)
                    max_Y.append(0)
                else:
                    max_l.append(
                        np.max(image[i, int(envelope[i, j]):int(envelope[i, j + 1])]))
                    max_Y.append(np.max(Y[i, :]))

            local_min_win = 15
            lenX = len(max_Y)
            for i in range(0, lenX, local_min_win):
                local_minY = np.min(max_Y[i:min(i+local_min_win, lenX)])
                max_Y[i:min(i+local_min_win, lenX-1)] -= local_minY
                local_minl = np.min(max_l[i:min(i+local_min_win, lenX)])
                max_l[i:min(i+local_min_win, lenX-1)] -= local_minl

            # max_Y = pickle.load(open(spectrum, 'rb'))
            # Then we create a second array which is sorted in order to get the 10th maximum ( not the first it may be a cosmic)
            max_l_order = np.copy(max_l)  # [max_l[k] for k in range(len(max_l))]
            max_l_order.sort()
            max_Y_order = np.copy(max_Y)  # [max_Y[k] for k in range(len(max_Y))]
            max_Y_order.sort()
            # Then, beacause the true data and the created one haven't the same intensity, we normalise them
            for i in range(lenX):
                max_l[i] = max_l[i] / max_l_order[len(max_l_order) - 10]
                max_Y[i] = max_Y[i] / max_Y_order[len(max_Y_order) - 10]
                diff.append(abs(max_l[i] - max_Y[i]) if max_Y[i] else 0)				# And we substract them ...
            # ... In order to calculate the standard deviation

            dev = np.std(diff)
            
            for i in range(lenX):
                # If, in the difference of the two images, one is greater than 5times the standard deviation, we assume it's a cosmic
                if diff[i] > 10 * dev:
                    # And we get rid of it
                    cos_loc = int(envelope[i, j]) + np.argmax(image[i,
                                                                    int(envelope[i, j]):int(envelope[i, j + 1])])
                    new_image[i, cos_loc] = (new_image[max(
                        i - 1, 0), cos_loc] + new_image[min(i + 1, lenX - 1), cos_loc]) / 2
                    print "Cosmic retire i = {} cos_loc =  {} diff = {}".format(i, cos_loc, diff[i])
            # And we store the new image
            # la il faut finir pour creer les fichiers avec les bonnes extensions tout ca tout ca....

    new_file = sys.argv[1]
    hdu = pyfits.PrimaryHDU(new_image)
    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(new_file, clobber=True)
    hdulist.close()

       
cfg = utils.CfgFile("../DATA/Amatrix_DATA.txt")

lanes_per_order = int(cfg.get_param("Lanes per order"))
envelope = cfg.get_param("Lane envelope file")

order = int(cfg.get_param("Order n°"))
lane = int(cfg.get_param("Lane n°"))
fts_file = cfg.get_param("Test fts file")

erase_cosmic(fts_file, envelope, order * lanes_per_order + lane - 1)
