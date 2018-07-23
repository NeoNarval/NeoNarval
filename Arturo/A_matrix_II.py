# -*- coding: utf-8 -*-

"""

"""

import os
import gc
import utils
import pyfits
import tables
import cPickle
import numpy as np
from scipy import sparse, linalg
from scipy.sparse import linalg as slinalg
import matplotlib.pyplot as plt

def rint(x): return int(round(x))


def sign(x): return int(np.sign(x))


def fill_A_fp(left_lane, thickness, ref_data, tab_centre, art, ini_index, end_index):
    pix = 1. / art

    (lenX, _) = ref_data.shape
    nb_col = rint(thickness)

    row, col, data = [], [], []
    search_window = 5 * pix

    for [x_pos, _] in tab_centre[0:10]:
        if (ini_index < x_pos < end_index):
            min_x = rint(max(ini_index, x_pos-search_window))
            max_x = rint(min(end_index, x_pos+search_window+1))

            for x in range(min_x, max_x):
                left = left_lane[x]
                for y in range(rint(left), rint(left+thickness)):
                    row.append((x-ini_index)*nb_col + y-rint(left))
                    col.append(rint((x_pos-ini_index) * pix))

                    data.append(float(ref_data[x, y]))

    A_shape = (lenX*nb_col, rint(lenX*pix))
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float')

    return A


def fill_A_matrix(left_lane, thickness, tab_centre, ref_data, ini_index, end_index):
    art = 1
    pix = 1. / art
    rpix = rint(pix)
    (lenX, _) = ref_data.shape
    nb_col = rint(thickness)

    A_fp = fill_A_fp(left_lane, thickness, ref_data, tab_centre, art, ini_index, end_index)
    print(A_fp.todense().shape)	
    print(np.shape(np.squeeze(A_fp.getrow(0).todense()))) 
    A=np.zeros((2000,4612))
    for i in np.arange(0,2000):
	A[i,:]=  np.squeeze(A_fp.getrow(i).todense())

    plt.imshow(A,aspect='auto')
    plt.show()
    
    u, s, v = slinalg.svds(A_fp, k=len(tab_centre) + 10)

    del(A_fp)
    gc.collect()

    C1 = np.matrix(v.transpose()) * np.matrix(linalg.pinv(np.diag(s)))
    U = np.matrix(u.transpose())

    A_fp_inv = np.dot(C1, U)

    row, col, data = [], [], []
    search_window = 5 * pix

    for i in range(0, len(tab_centre) - 1):
        x0, x1 = tab_centre[i][0], tab_centre[i+1][0]

        if (ini_index < x0 < end_index and ini_index < x1 < end_index):
            ax0, ax1 = pix * x0, pix * x1
            arx0, arx1 = rint(ax0), rint(ax1)

            for j in range(arx0+1, arx1):
                g = j
                j = rint(j*art)
                min_x = rint(max(ini_index, g - search_window))
                max_x = rint(min(end_index * rpix, g + search_window+1))

                for x in range(min_x, max_x):
                    for y in range(thickness):

                        row.append(g - rint(pix * ini_index))
                        col.append(rint((x-ini_index)*art)*nb_col + y)

                        pos_0 = rint((x-g)*art + x0 - (1+art) * ini_index)*nb_col + y
                        pos_1 = rint((x-g)*art + x1 - (1+art) * ini_index)*nb_col + y

                        p0 = A_fp_inv[arx0, min(rint(lenX*thickness)-1, rint(pos_0))]
                        p1 = A_fp_inv[arx1, min(rint(lenX*thickness)-1, rint(pos_1))]

                        coeff = float(g-arx0) / (rint(x1-x0) * rpix)
                        p = (p0 + coeff * float((p1-p0))) * art

                        data.append(p)

    A_shape = (rint(lenX*pix), rint(lenX*thickness))
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float').tolil()

    for [x_pos, _] in tab_centre:
        if (ini_index < x_pos < end_index):
            min_x = rint(max(0, x_pos-search_window))
            max_x = rint(min(lenX, x_pos+search_window+1))

            for x in range(min_x, max_x):
                for y in range(thickness):
                    A[rint(x_pos * pix), x*nb_col +
                      y] = A_fp_inv[rint(x_pos * pix), x*nb_col + y]

    A_full = A.tocsr()

    del(A)
    del(A_fp_inv)
    gc.collect()

    return A_full


def cut_border_edge(left_lane):
    lenX = len(left_lane)
    init = 0
    end = lenX - 1
    test = False

    for i in range(lenX):
        if not left_lane[i]:
            if not test and i >= init:
                init = min(i+1, end)
            elif test and i >= init:
                end = max(init, i-1)
                break
        else:
            test = True

    return (init, end)


def find_path(path_file, key):
    """
        A method to find the path of the directory we want in a txt file
            path_file   : the path of the file which contains all the paths for all the directory
            key         : the name of the directory we want the path
        return : 
            value       : the value corresponding to the key
    """

    file = open(path_file, "r")
    text = file.read().split("\n")
    for j in range(len(text)):
        # We assume that in the file : "Key : value"
        line = text[j].split(" : ")
        if line[0] == key:
            return line[1]

    return ("not found")



    """
        Returns an array containing the FP centroids of all considered orders.
        it also gives a colourized of the image_data to have a quick overview
        of the result.
    """
    # path_file = "Amatrix_DATA.txt"

cfgFile = utils.CfgFile("./Amatrix_DATA.txt") #ALA: use my own config file

thickness_data_file = cfgFile.get_param("Lane thickness file")
envelope_data_file = cfgFile.get_param("Lane envelope file")
fp_position_file = cfgFile.get_param("FP slits position file")
fp_file = cfgFile.get_param("Fabry Perot fts file")
order = int(float(cfgFile.get_param("Order n°")))
lane = int(float(cfgFile.get_param("Lane n°")))

    # thickness_data_file = find_path(path_file, "Lane thickness file")
    # envelope_data_file = find_path(path_file, "Lane envelope file")
    # fp_position_file = find_path(path_file, "FP slits position file")
    # fp_file = find_path(path_file, "Fabry Perot fts file")
    # order = int(float(find_path(path_file, "Order n°")))
    # lane = int(float(find_path(path_file, "Lane n°")))

pos_global = cPickle.load(open(fp_position_file, 'rb'))
envelope_data = cPickle.load(open(envelope_data_file, 'rb'))
thickness_data = cPickle.load(open(thickness_data_file, 'rb'))

image_file = pyfits.open(fp_file)
fp_data = image_file[0].data.astype(np.float32)
image_file.close()

left = envelope_data[:,2*order + lane - 1]
tab_pos = pos_global[2*order + lane - 1]
thickness = rint(thickness_data[2*order + lane - 1])

#print(np.shape(tab_pos))	
#plt.plot(tab_pos[0:30],'o')
#plt.show()	
    		

(ini, end) = (2000,3000)#cut_border_edge(left)
img=np.zeros((end-ini+1,thickness))
for i in range(ini,end):
	img[i-ini,:]=fp_data[i,rint(left[i]):rint(left[i]+thickness)]


A_matrix=np.zeros(((end-ini+1)*thickness,end-ini+1))
for [xpos,_] in tab_pos:
	if (ini<xpos<end):
		B=np.copy(img)*0.
		wvl=rint(xpos)-ini
		B[wvl-5:wvl+5,:]=img[wvl-5:wvl+5,:]
		A_matrix[:,wvl]=np.hstack(B)

u, s, v = slinalg.svds(A_matrix, k=len(tab_pos) + 10)

C1 = np.matrix(v.transpose()) * np.matrix(linalg.pinv(np.diag(s)))
U = np.matrix(u.transpose())

A_fp_inv = np.dot(C1, U)


for i,[xpos,_] in enumerate(tab_pos):
	if (ini<xpos<end):
		wvl0=rint(xpos)-ini
		wvl1=min(rint(tab_pos[i+1][0])-ini,len(A_fp_inv)-1)	
		A0=np.squeeze(np.asarray(A_fp_inv[wvl0,:]))
		A1=np.squeeze(np.asarray(A_fp_inv[wvl1,:]))
		pend=(np.roll(A1,-thickness*(wvl1-wvl0))-A0)/(wvl1-wvl0)
		for wvl in np.arange(wvl0+1,wvl1):
			A_fp_inv[wvl,:]=np.roll(A0+pend*(wvl-wvl0),thickness*(wvl-wvl0))



spectra=np.squeeze(np.asarray(A_fp_inv.dot(np.hstack(img))))
stop
#A_matrix = fill_A_matrix(left, thickness, tab_pos, fp_data, ini, end)

#name_A_matrix = './Amatrix_' + 'OR' + str(order) + '_LA' + str(lane) + '.npz'
#sparse.save_npz(name_A_matrix, A_matrix)

#cfgFile.modif_param('A matrix', name_A_matrix)


