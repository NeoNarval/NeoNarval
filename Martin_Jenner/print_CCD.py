import numpy as np
import pyfits
import collections
import numpy.core.multiarray as multiarray
import cPickle
from multiprocessing import Pool
import matplotlib.pyplot as plt
"""
a way to visualize a ccd
"""
def rint(x): return int(round(x))

def sum_on_y_axis(left_lane, thickness, ini_index, end_index):
    lenX = end_index - ini_index
    intensity = np.zeros(lenX)
    nb_col = rint(thickness)
    
    for i in range(ini_index, end_index):
        intensity[i-ini_index] = sum(CCD_data[i][min(rint(left_lane[i]+j), lenDataY-1)] for j in range(nb_col))
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
    
def launch(j):
    """
    This function launch the appropriate functions to print an image of the ccd
    for one lane. 
    : j : the considered lane
    
    : return : an array with all the intensity for each column in a lane        
    """
    left_lane = envelope_data[:, j]       # Upper envelope of the current lane
    thickness = rint(thickness_data[j])   # Thickness of the current lane
    
    (ini, end) = cut_border_edge(left_lane,thickness)
    
    return sum_on_y_axis(left_lane, thickness, ini, end)
    
    
def print_CCD(path,j):
    
    
    global envelope_data    # Upper envelope of each lane
    global thickness_data   # Thickness of the lanes
    global nbr_lanes        # Number of lanes per order
    global lenDataX         # Dimension of the CCD (length)
    global lenDataY         # Dimension of the CCD (width)
    global CCD_data         # The CCD data

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

    CCD_file     = dic["CCD fts file"]
    nbr_lanes   = int(dic["Lanes per order"])

    envelope_data  = cPickle.load(open(envelope_data_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))

    image_file = pyfits.open(CCD_file)
    CCD_data = image_file[0].data.astype(np.float32) # Data of the CCD fts file
    image_file.close()

    (lenDataX, lenDataY) = CCD_data.shape

    init_order, final_order = 0, 39    # We assume there are 40 orders

    # Multiprocessing : each lane is independant from the others

    # if __name__ == '__main__':
    #     p = Pool(10)
    #     intensity_global = p.map(launch, range(nbr_lanes*init_order, nbr_lanes*final_order + nbr_lanes, 1))
    #     p.terminate()
    #     p.join()
    
    # for j in range(nbr_lanes*init_order, nbr_lanes*final_order + nbr_lanes, 1):
    #     plt.plot(launch(j))
    # plt.show()
    plt.plot(launch(j))
    plt.show()
    
j = input("numéro de l'ordre")
path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\CCD_data_sheet.txt"
print_CCD(path, j)