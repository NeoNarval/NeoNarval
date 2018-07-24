import pyfits
import cPickle
import numpy as np
import collections
from scipy import sparse
import matplotlib.pyplot as plt

plt.close(101)
plt.close(102)
plt.close(103)
def rint(x): return int(round(x))

"""
eta_stop : critere d'arret depandant de si flat ou pas
"""
def landweber(CCD, eta_stop,  B_matrix):
    
    U,max_sig,V = sparse.linalg.svds(B_matrix,1)
    
    iter = 0
    w = 1./max_sig[0] # < 2/sigma_max
    print("relaxation = ",w)
    
    n = np.shape(CCD)[0]
    sigma = 1./np.sqrt(np.asarray(CCD).ravel())
    
    sig = sparse.spdiags(sigma.T, 0, n,n)
    #print(np.shape(B_matrix.T), np.shape(sig))
    P = B_matrix.T.dot(sig)
    x_old = np.zeros((np.shape(B_matrix)[1],1))
    x = x_old.copy()
    #print(CCD)
    diff=1.
    diff0=0. 
    Morozov=0.1*np.median(sigma)  #Value of tau=1. for the Moon
    print("morozov = ", Morozov)
    
    while (diff>Morozov) and (np.abs(diff-diff0)>eta_stop) and iter < 5000:
        diff0=diff
        iter+=1
        M = CCD - B_matrix.dot(x_old)
        x = x_old + w*P.dot(M)
        x_old = x
        diff=np.linalg.norm(M)/np.linalg.norm(x)
    print('diff',diff,'delta_diff',np.abs(diff-diff0),iter)
    
    
    return x
    

def generate_spectrum(test_data, Bscale, Bzero, eta_stop, B_matrix, Meth):
    lambd_max = int(end_index*pix - ini_index)
    lambd_min = int(ini_index*pix - ini_index)
    B_matrix = B_matrix[:, lambd_min:lambd_max]
    CCD = []
    for i in range(ini_index, end_index+1):
        for j in range(thickness):
            CCD.append(((test_data[i, rint(left_lane(i)+j)]-Bzero)/Bscale))
    # CCD = (CCD - np.min(CCD))/np.max(CCD)
    # plt.figure(50)
    # plt.plot(CCD)
    # plt.show()
    CCD = np.matrix(CCD).T
    cPickle.dump(CCD, open(r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\second_membre_CCD', 'wb'))
    print(np.shape(CCD))
    
    
    if Meth == 'LW':
        Spectrum = landweber(CCD, eta_stop, B_matrix)
    
    
    return Spectrum
    

def get_spectrum(test, origin, ord):
    global nbr_lanes
    global ini_index    # first index of the window of ThAr to be processed
    global end_index    # last index of the window of ThAr to be processed
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global order        # Considered order
    global lane         # Considered lane within the order
    global lenX
    global pix
    global Bzero        # BZERO of the fits considered
    global Bscale       # BSCALE of the fits considered
    global fBzero       # BZERO of the flat fits considered
    global fBscale      # BSCALE of the flat fits considered
    Bzero = 0
    Bscale = 1
    fBzero = 32768
    fBscale = 1
    order = ord
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
    
    # polys_data_file      = dic["n polys env"]             
    # thickness_data_file = dic["n lanes thic"]
    # polys_data_file      = dic["Polys envelope file"]             
    # thickness_data_file = dic["Lane thickness file"]
    polys_data_file      = dic["a polys env"]             
    thickness_data_file = dic["a lanes thic"]
    
    #-------------------------------------------------
    
    if test == 'ThAr':
        test_file       = dic["ThAr fts file"]
    elif test == 'test':
        test_file       = dic["test file"]
    elif test == 'FP':
        test_file       = dic["FP fts file"]
    elif test =='flat':
        test_file       =dic["flat file"]
    elif test == 'star':
        test_file       =dic["star file"]
    flat_file           =dic["flat file"]
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    pix                 = float(dic["pix/arturo"])
    art                 = 1./pix
    
    B_matrix_name = r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\Bmatrix\Bmat_'+origin+'_ord'+str(order)+r'_lane'+str(lane)+r'.npz'
    dic["B matrix"] = B_matrix_name
    dic["B matrix origin"] = origin
    dic["order"] = str(order)
    with open(path, 'w') as file:
        for value in dic.items():
            line = value[0] + " : " + value[1] + "\n"
            file.write(line)
            
            
    B_matrix_path       = dic["B matrix"]

    B_matrix            = sparse.load_npz(B_matrix_path).tocsr()
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    thickness_float = thickness_data[nbr_lanes * order + lane -1]
    thickness =rint(thickness_float)
    print(thickness)
    
    polys_data     = cPickle.load(open(polys_data_file, 'r'))
    left_lane = polys_data[order]       # A polynomial giving the upper enveloppe of the lane
    left_lane[0] += thickness_float * (lane - 1)
    if test == 'test':
        test_data = cPickle.load(open(test_file, 'r'))
    else:   
        image_file = pyfits.open(test_file)
        image_header = image_file[0].header
        #Bzero = image_header['BZERO']
        #Bscale = image_header['BSCALE']
        test_data = image_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
        image_file.close()
    
    image_f_file = pyfits.open(flat_file)
    image_f_header = image_f_file[0].header
    #print(image_f_header)
    # fBzero = image_f_header['BZERO']
    # fBscale = image_f_header['BSCALE']
    flat_data = image_f_file[0].data.astype(np.float32) # Data of the ccd of the flat fts file
    image_f_file.close()
    
    lenX = test_data.shape[0]
    min = np.min(test_data)
    test_data = test_data-min
    #------------flat-field handling---------
    
    
    flat_field = np.asarray(generate_spectrum(flat_data, fBscale, fBzero, 1e-3, B_matrix, 'LW')).ravel()
    plt.figure(102)
    plt.plot(flat_field)
    plt.show()
    
    Spectrum = np.asarray(generate_spectrum(test_data, Bscale, Bzero, 1e-2, B_matrix, 'LW')).ravel()
    plt.figure(101)
    plt.plot(Spectrum)
    plt.show()
    
    Spectrum = Spectrum/flat_field
    plt.figure(103)
    plt.plot(Spectrum)
    plt.show()
    #-----------------------------------------
    #-------lane pos check------------
    plt.matshow(np.array([flat_data[i,int(np.floor(left_lane(i))):int(np.floor(left_lane(i)))+thickness] for i in range(4612)]),aspect = 'auto')
    plt.show()
    #--------------------------------
    pickle_name = r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\ThAr_based_spec'+ '_OR'+str(order)+'_LA'+str(lane)+'.p'
    cPickle.dump(Spectrum, open(pickle_name, 'wb'))
    
get_spectrum('ThAr', 'simu', 14)