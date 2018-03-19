import pyfits
import cPickle
import numpy as np
from scipy import sparse
import collections
import matplotlib.pyplot as plt

def rint(x): return int(round(x))


def fill_B_ThAr(ref_data):
    
    row, col, data = [], [], []
    search_window = 0
    rng_thickness = range(thickness)
    
    for [x_pos, _] in tab_centre:
        if (ini_index<= x_pos <= end_index):
            min_x = rint(max(ini_index, x_pos-search_window))
            max_x = rint(min(end_index + 1, x_pos+search_window +1))
            
            wavelength = rint(x_pos * pix) # wavelength in arturos
            delta_x_pos = x_pos - rint(x_pos) #error due to round number
            
            for x in range(min_x, max_x):
                for y in rng_thickness:
                    row.append(rint((x+delta_x_pos)*thickness + y))
                    col.append(wavelength)
                    data.append(float(ref_data[rint(x + delta_x_pos), rint(y+left_lane(x))]))
                    
    B_shape = (lenX*thickness, rint((end_index-ini_index)*pix))
    B = sparse.coo_matrix((data, (row, col)), shape=B_shape, dtype='float')

    return B


def interpolate(i):
    """
        Interpolate the A matrix between two ThAr slits.

        : i : the index of the considered left slit in the slit-position tab

        : return :  the information for each interpolated pixel (row, column and intensity) (list format)
                    and the first and last wavelength (in pixel, not arturos)
    """
    row, col, data = [], [], []

    x0, x1 = tab_centre[i][0], tab_centre[i+1][0]   # We consider two successive slits of the ThAr.

    if (x0 < end_index and x1 < end_index+1):
        arx0, arx1 = rint(x0 * pix), rint(x1 * pix)                 # Rounded wavelength of the slits in arturos
        dy = -float(left_lane(x1) - left_lane(x0)) / float(x1-x0)   # Gap on the Y axis between the slits because of the curved lane

        # Those values are used after. They are created here to reduce the complexity.
        search_window = 0
        f_diff = float(x1 - x0)
        rng_thickness = range(thickness)
        frac_dy = (1. / (1 + 2 * (dy ** 2)))

        # We successively consider each wavelength (hence each column) between the two slits.
        
        for j in range(arx0, arx1):  #the ThAr wavelength are interpolated (continuous)

            # Minimum and maximum indices of the search window
            min_x = rint(max(ini_index * pix, j - search_window))
            max_x = rint(min((end_index + 1) * pix, j + search_window + 1))

            # Coefficient for the interpolation (equals to 0 at lambda0 and 1 at lambda1)
            coeff = art * float(j-arx0) / f_diff
            test_tab = []

            for x in range(min_x, max_x):   # Cross of the search window on the X axis

                for y in rng_thickness:     # Cross of the search window on the Y axis
                    row_value = rint(x*art*thickness + y - left_lane(x) + rint(left_lane(i)))   # Considered pixel

                    if (row_value not in test_tab):
                        pos_0 = rint( ((x-j) * art + x0) * thickness ) + y # Pixel equivalent positon in the ThAr slit column (left one)
                        pos_1 = rint( ((x-j) * art + x1) * thickness ) + y # Pixel equivalent positon in the ThAr slit column (right one)

                        p0 = B_ThAr[pos_0, arx0]  # Pixel participation intensity (left ThAr column)
                        p1 = B_ThAr[pos_1, arx1]  # Pixel participation intensity (rigth ThAr column)

                        # Quadratic interpolation
                        p1_sup1 = B_ThAr[pos_1 + 1, arx1] # Intensity participation of the upper pixel of the rigth ThAr column
                        p1_inf1 = B_ThAr[pos_1 - 1, arx1] # Intensity participation of the lower pixel of the rigth ThAr column

                        p1_sup2 = B_ThAr[pos_1 + 2, arx1]
                        p1_inf2 = B_ThAr[pos_1 - 2, arx1]

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

    return [col, row, data]
    
def fill_B_matrix(ref_data):
    """
        Creates the full B matrix

        : ref_data : The ThAr data

        : return : The full interpolated B matrix
    """
    
   
    # Interpolation between the slits with multiprocessing
    # p = Pool(10)
    # DATA = p.map(interpolate, range(0, len(tab_centre) - 1))
    # p.terminate()
    # p.join()
    DATA = []
    for i in range(0, len(tab_centre)-1):
        DATA.append(interpolate(i))
    zip_data = zip(*DATA)
    col = sum(list(zip_data[0]), [])
    row = sum(list(zip_data[1]), [])
    data = sum(list(zip_data[2]), [])    
    
    rp_lenX = rint((end_index-ini_index) * pix) # Round number of arturos in one lane
    # Creation of the full interpolated matrix with the header, in CRS format
    B_shape = (rint((end_index-ini_index)*thickness), rp_lenX)
    B = sparse.coo_matrix((data, (row, col)), shape=B_shape, dtype='float')
    B_full = B.tocsr()

    return B_full
    
def create_B_matrix():
    global nbr_lanes
    global ini_index    # first index of the window of ThAr to be processed
    global end_index    # last index of the window of ThAr to be processed
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global tab_centre   # positions of the ThAr slits
    global order        # Considered order
    global lane         # Considered lane within the order
    global lenX         # Dimension of the CCD (length)
    global art          # Nu mber of pixel per arturo
    global pix          # Number of arturos per pixel : 1.5 art/pix <=> 1.73km/s / pix 
    global B_ThAr       # B matrice only filled for ThAr wavelengths
    path = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\Bmatrix_data_sheet.txt'
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
                 
    thickness_data_file = dic["Lane thickness file"]
    polys_data_file     = dic["Polys envelope file"]
    ThAr_position_file  = dic["ThAr slits file"]
    ThAr_file           = dic["ThAr fts file"]
    order               = int(dic["order"])
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    pix                 = float(dic["pix/arturo"])
    art                 = 1./pix
    polys_data     = cPickle.load(open(polys_data_file, 'r'))
    tab_centre     = cPickle.load(open(ThAr_position_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    
    thickness_float = thickness_data[nbr_lanes * order + lane - 1]
    thickness = rint(thickness_float)   # The thickness of the considered lane
    print(thickness)
    
    left_lane = polys_data[order]       # A polynomial giving the upper enveloppe of the lane
    left_lane[0] += thickness_float * (lane - 1)

    image_file = pyfits.open(ThAr_file)
    ThAr_data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
    image_file.close()
    [lenX, _]           = ThAr_data.shape
    
        
    B_ThAr = fill_B_ThAr(ThAr_data).tocsr()  # Creation of the B matrix only filled with ThAr
    B_matrix = fill_B_matrix(ThAr_data)
    B_matrix_name = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\Bmatrix\Bmat_ord'+str(order)+r'_lane'+str(lane)+r'.npz'
    sparse.save_npz(B_matrix_name, B_matrix)    # B matrix saved with sparse format
    dic["B matrix"] = B_matrix_name
    with open(path, 'w') as file:
        for value in dic.items():
            line = value[0] + " : " + value[1] + "\n"
            file.write(line)
    
create_B_matrix()
    
    
    
    
    
    
    
    
    
    