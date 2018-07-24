import pyfits
import cPickle
import numpy as np
from scipy import sparse
import collections
import matplotlib.pyplot as plt
import math as m

def rint(x): return int(round(x))


def fill_B_ThAr(ref_data):
    
    row, col, data = [], [], []
    #print(thickness)
    rng_thickness = range(thickness)
    #print(first_pix)
    #print(len(tab_centre))
    #print(len(tab_centre),len(tab_centre[0]))
    to_del = []
    for i in range(1,len(tab_centre)-1):
        x_pos = tab_centre[i]
        x_pos_prev = tab_centre[i-1]
        x_pos_next = tab_centre[i+1]
        if (abs(x_pos[0]-x_pos_prev[0])<2*search_window) :
            if not to_del.__contains__(x_pos):
                to_del.append(x_pos)
            if not to_del.__contains__(x_pos_prev):
                to_del.append(x_pos_prev)
        elif (abs(x_pos[0]-x_pos_next[0])<2*search_window):
            if not to_del.__contains__(x_pos):
                to_del.append(x_pos)
            if not to_del.__contains__(x_pos_next):
                to_del.append(x_pos_next)
    for x in to_del:
        tab_centre.remove(x)
    
    min_data = 0
    for [x_pos, y_pos] in tab_centre:
        if x_pos == -1:
            min_data = y_pos/thickness
            #print('min_data= ', min_data)
        if (ini_index<= x_pos <= end_index):
            
            min_x = rint(max(ini_index, x_pos-search_window))
            max_x = rint(min(end_index, x_pos+search_window)+1)
            
            wavelength = rint(x_pos * pix)# wavelength in arturos
            delta_x_pos = x_pos - rint(x_pos) #error due to round number
            #print(delta_x_pos)
            #delta_x_pos = 0
            
            #print(min_x, max_x, x_pos)
            
            for x in range(min_x, max_x):
                X = rint(x+delta_x_pos)
                curv_pos = int(np.floor(left_lane(x)))
                # if ((ref_data[X,curv_pos]-min_data)/(np.mean(ref_data[rint(x_pos),curv_pos:curv_pos+thickness])-min_data))<=0.0001:
                #     print("!", curv_pos)
                #     curv_pos+=1
                for y in rng_thickness:
                    row.append(rint((x-ini_index)*thickness + y))
                    col.append(rint(wavelength-ini_index*pix))

                    Y = y+curv_pos
                    d = float(((ref_data[X,Y]-Bzero)/Bscale)/(np.sum((ref_data[rint(x_pos-search_window):rint(x_pos+search_window),curv_pos:curv_pos+thickness]-Bzero)/Bscale))) 
                    # d = float(((ref_data[X,Y]-np.min(ref_data[rint(x_pos-search_window):rint(x_pos+search_window),curv_pos:curv_pos+thickness])-Bzero)/Bscale)/(np.max((ref_data[rint(x_pos-search_window):rint(x_pos+search_window),curv_pos:curv_pos+thickness]-Bzero)/Bscale)))
                    # if d<=0.0001:
                    #     d = 0
                    data.append(d)
                    print(x,y,curv_pos,rint((x-ini_index)*thickness + y),rint(wavelength-ini_index*pix),d)
                    #data.append(float(ref_data[rint(x + delta_x_pos), rint(y+left_lane(x))]))
                    #print(float(ref_data[rint(x + delta_x_pos), rint(y+left_lane(x))]), x+delta_x_pos, rint(x+delta_x_pos), y+left_lane(x),rint(y+left_lane(x)))
                
        
    #print(max(row), max(col))
    # row = row[search_window:-search_window]
    # col = col[search_window:-search_window]
    # data = data[search_window:-search_window]
    #print(row)
    
    B_shape = (rint((end_index-ini_index+1)*thickness), rint((end_index-ini_index)*pix))
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

    x0, x1 = tab_centre[i][0], tab_centre[i+1][0] # We consider two successive slits of the ThAr.
    print(i, x0, x1)
    x0 -= ini_index
    x1 -= ini_index
    
    
    if (x0 < end_index-ini_index and x1 < end_index+1-ini_index and x0 != x1):
        arx0, arx1 = rint(x0 * pix), rint(x1 * pix)                 # Rounded wavelength of the slits in arturos
        dy = -float(left_lane(x1) - left_lane(x0)) / float(x1-x0)   # Gap on the Y axis between the slits because of the curved lane
        
        # Those values are used after. They are created here to reduce the complexity. 
        f_diff = float(x1 - x0)
        rng_thickness = range(thickness)
        frac_dy = (1. / (1 + 2 * (dy ** 2)))
        
        # settings to extrapolate at the edge of the matrix
        if i==1:
            interp_window0 = 0
        else:
            interp_window0 = arx0
        if i==len(tab_centre)-2:
            interp_window1 = int(np.floor((end_index- ini_index) * pix))
        else:
            interp_window1 = arx1
           
            
        # We successively consider each wavelength (hence each column) between the two slits.
        # for j in range(arx0 + 1, arx1):     # Case 1 : real value of the F-P are kept (discontinuity)
        for j in range(interp_window0, interp_window1):         # Case 2 : the F-P wavelenfgth are interpolated (continuous)
            
            
            # Minimum and maximum indices of the search window

            min_x = rint(max(0, j - search_window))
            max_x = rint(min((end_index- ini_index) * pix+1, j + search_window+1))

            # Coefficient for the interpolation (equals to 0 at lambda0 and 1 at lambda1)
            coeff = art * float(j-arx0) / f_diff
            test_tab = []

            for x in range(min_x, max_x):   # Cross of the search window on the X axis
                # row_v = rint(x*art*thickness) 
                
                for y in rng_thickness:     # Cross of the search window on the Y axis
                    # row_value = row_v + y   # Considered pixel
                    row_value = x*art*thickness + y   # Considered pixel
                    index_row_value = int(np.ceil(row_value))
                    delta_row_value = row_value - index_row_value
                     #- left_lane(x) + left_lane(i)
                    if (row_value not in test_tab and row_value<(end_index-ini_index)*thickness):
                        pos_0 = rint( ((x-j) * art + x0) * thickness ) + y # Pixel equivalent positon in the F-P slit column (left one)
                        pos_1 = rint( ((x-j) * art + x1) * thickness ) + y # pixel equivalent positon in the f-p slit column (right one)
                        #print(pos_0,pos_1,x,y,min_x,max_x,j)
                        if pos_0<0:
                            p0 = 0.0
                        else:
                            p0 = B_ThAr[pos_0, arx0]  # Pixel participation intensity (left F-P column)
                        if pos_1>=(end_index-ini_index)*pix*thickness:
                            p1 = 0.0
                        else:
                            p1 = B_ThAr[pos_1, arx1]  # Pixel participation intensity (rigth F-P column)

                        #---------------Quadratic interpolation----------------
                        #p1_sup1 = B_ThAr[pos_1 + 1, arx1] # Intensity participation of the upper pixel of the rigth F-P column
                        #p1_inf1 = B_ThAr[pos_1 - 1, arx1] # Intensity participation of the lower pixel of the rigth F-P column

                        #p1_sup2 = B_ThAr[pos_1 + 2, arx1]
                        #p1_inf2 = B_ThAr[pos_1 - 2, arx1]

                        #p1_tilde_i     = frac_dy * (p1      + dy * (p1_inf1 - p1_sup1 ))    # Expected intensity of the pixel
                        #p1_tilde_i_inf = frac_dy * (p1_inf1 + dy * (p1_inf2 - p1))          # Expected intensity of the lower pixel
                        #p1_tilde_i_sup = frac_dy * (p1_sup1 + dy * (p1 - p1_sup2))          # Expected intensity of the upper pixel

                        #a = -dy * (p1_tilde_i_sup - p1_tilde_i_inf) / f_diff
                        #p = float(p0 + coeff * float(abs(a + p1_tilde_i) - p0)) # quadratic Interpolated intensity
                        #------------------------------------------------------
                        
                        p = p0 + (p1-p0)*coeff                                  # linear Interpolated intensity
                        # Column, row, and intensity of the current pixel added
                       
                        last_index = len(row)-1
                        
                        if (row != []) and (row[last_index] == index_row_value and col[last_index] == j):
                            data[last_index]+=p*delta_row_value
                            if index_row_value<rint(end_index-ini_index)*thickness-1:
                                col.append(j)
                                row.append(index_row_value+1)
                                data.append(p*(1-delta_row_value))
                        else :
                            col.append(j)
                            row.append(int(index_row_value))
                            data.append(p*delta_row_value)
                            if index_row_value<rint(end_index-ini_index)*thickness-1:
                                col.append(j)
                                row.append(index_row_value+1)
                                data.append(p*(1-delta_row_value))

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
    
    #print(len(tab_centre), tab_centre)
    for i in range(1, len(tab_centre)-1):
        DATA.append(interpolate(i))
    
    zip_data = zip(*DATA)

    
    col = sum(list(zip_data[0]), [])
    row = sum(list(zip_data[1]), [])
    data = sum(list(zip_data[2]), [])  
    #print(len(row),len(col),len(data))
    rp_lenX = rint((end_index-ini_index) * pix) # Round number of arturos in one lane
    
    #print(max(row),max(col))
    #print(rint((end_index-ini_index)*thickness),rp_lenX)
    # Creation of the full interpolated matrix with the header, in CRS format
    B_shape = (rint((end_index-ini_index+1)*thickness), rp_lenX)
    B = sparse.coo_matrix((data, (row, col)), shape=B_shape, dtype='float')
    B_full = B.tocsr()

    return B_full
    
def create_B_matrix(test):
    global nbr_lanes
    global ini_index    # first index of the window of ThAr to be processed
    global end_index    # last index of the window of ThAr to be processed
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global tab_centre   # positions of the ThAr slits
    global order        # Considered order
    global lane         # Considered lane within the order
    global art          # Nu mber of pixel per arturo
    global pix          # Number of arturos per pixel : 1.5 art/pix <=> 1.73km/s / pix 
    global B_ThAr       # B matrice only filled for ThAr wavelengths
    global first_pix    # first pixel of the slab of order considered
    global Bzero        # BZERO of the fits considered
    global Bscale       # BSCALE of the fits considered
    global search_window
    Bzero = 0
    Bscale = 1
    search_window = 11
    path = r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt'
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
         
    #--------------postion of the lanes---------------
    
    # thickness_data_file = dic["n lanes thic"]
    # polys_data_file     = dic["n polys env"]
    # thickness_data_file = dic["Lane thickness file"]
    # polys_data_file     = dic["Polys envelope file"]
    thickness_data_file = dic["a lanes thic"]
    polys_data_file     = dic["a polys env"]
    
    #------------------------------------------------
    
    position_file       = dic["position file"]
    if test =='test':
        position_file   = dic["test position file"]
        
    ThAr_file           = dic["ThAr fts file"]
    test_file           = dic["test file"]
    FP_file             = dic["FP fts file"]
    simu_file           = dic["simu fts file"]
    flat_file           = dic["flat file"]
    order               = int(dic["order"])
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    pix                 = float(dic["pix/arturo"])
    art                 = 1./pix

    polys_data     = cPickle.load(open(polys_data_file, 'r'))
    tab_centre     = cPickle.load(open(position_file, 'r'))
    #print(tab_centre)
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    #print(thickness_data)
    thickness_float = thickness_data[nbr_lanes * order + lane - 1]
    thickness       = rint(thickness_float)   # The thickness of the considered lane
    #print(thickness)
    global data
    left_lane = polys_data[order]       # A polynomial giving the upper enveloppe of the lane
    left_lane[0] += thickness_float * (lane - 1)
    if test=='ThAr':     
        image_file = pyfits.open(ThAr_file)
        image_header = image_file[0].header
        # Bzero = image_header['BZERO']
        # Bscale = image_header['BSCALE']
        data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()
    elif test=='test':
        data = cPickle.load(open(test_file, 'r'))
    elif test=='FP':
        image_file = pyfits.open(FP_file)
        image_header = image_file[0].header
        # Bzero = image_header['BZERO']
        # Bscale = image_header['BSCALE']
        data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()
    elif test=='simu':
        image_file = pyfits.open(simu_file)
        image_header = image_file[0].header
        # Bzero = image_header['BZERO']
        # Bscale = image_header['BSCALE']
        data = image_file[0].data.astype(np.float32) # Data of the ThAr fts file
        image_file.close()
    # plt.figure(75)
    # plt.matshow(np.array([data[i,int(np.floor(left_lane(i))):int(np.floor(left_lane(i)))+thickness] for i in range(4612)]),aspect = 'auto')
    # plt.show()
    B_ThAr = fill_B_ThAr(data).tocsr()  # Creation of the B matrix only filled with ThAr
    
    # plt.matshow(B_ThAr.toarray(), aspect = 'auto')
    # plt.show()
    
    B_matrix = fill_B_matrix(data)
    
    plt.matshow(B_matrix.toarray(), aspect ='auto')
    plt.show()
    #------ print interpolation ------
    # B= B_matrix.toarray()
    # sum = np.sum(B,axis=0)
    # plt.figure()
    # plt.plot(sum)
    # plt.show()
    #-------------------------------------
    #------print one line--------
    [n,m] = np.shape(B_matrix)
    plt.figure()
    plt.plot(B_matrix.toarray()[:, m//2])
    plt.show()
    plt.figure()
    plt.plot(B_matrix.toarray()[n//2,:])
    plt.show()
    #-------------------------
    #print(Bzero, Bscale)
   
    B_origin = test
    
    B_matrix_name = r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\Bmatrix\Bmat_'+B_origin+'_ord'+str(order)+r'_lane'+str(lane)+r'.npz'
    print(B_matrix_name)
    sparse.save_npz(B_matrix_name, B_matrix)    # B matrix saved with sparse format
    dic["B matrix"] = B_matrix_name
    dic["B matrix origin"] = B_origin
    with open(path, 'w') as file:
        for value in dic.items():
            line = value[0] + " : " + value[1] + "\n"
            file.write(line)

# test = 'ThAr' : ThAr
# test = 'test' : CCD random issue de CCD _creator
# test = 'FP' : FP
create_B_matrix('simu')
    
    
    
    
    
    
    
    
    
    