import pyfits
import cPickle
import numpy as np
import collections
from scipy import sparse
import matplotlib.pyplot as plt

# Method used to get the closest integer of a real number
def rint(x): return int(round(x))

def generate_spectrum(test_data, B_matrix):
    
   
    
    # B_matrix[B_matrix.nonzero()] = 1
    # plt.matshow(B_matrix.toarray(), aspect = 'auto')
    # plt.show()
    # 
    lambd_max = int(end_index*pix - ini_index)
    lambd_min = int(ini_index*pix - ini_index)
    
    B_matrix = B_matrix[:, lambd_min:lambd_max]
    
    CCD = []
    for i in range(ini_index, end_index):
        for j in range(thickness):
            CCD.append(test_data[i, rint(left_lane(i)+j)])
    CCD = np.matrix(CCD).T
    print(len(CCD))
    
    
 # We use the method of Normal Equations to retrieve the spectrum of the lane
    C = np.dot(B_matrix.T, B_matrix)            # C = At * A
    C = C + np.identity(C.shape[0]) * 0.000001  # We get rid of the potential null values on the diagonal

    G = np.linalg.cholesky(C)   # Cholesky factorization
    d = (B_matrix.T).dot(CCD)
    y = np.linalg.solve(G, d)
    X = np.linalg.solve(G.T, y)
    spectrum = (list(np.array(X.T.tolist()[0])))
    
    return spectrum
    
def get_spectrum(test):
    global nbr_lanes
    global ini_index    # first index of the window of ThAr to be processed
    global end_index    # last index of the window of ThAr to be processed
    global thickness    # thickness of the considered lane
    global left_lane    # polynomial giving the upper enveloppe of the lane
    global order        # Considered order
    global lane         # Considered lane within the order
    global lenX
    global pix

    path = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt'
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
    
    polys_data_file      = dic["Polys envelope file"]             
    thickness_data_file = dic["Lane thickness file"]
    if test == 'ThAr':
        test_file       = dic["ThAr fts file"]
    elif test == 'test':
        test_file       = dic["test file"]
    elif test == 'FP':
        test_file       = dic["FP fts file"]
    elif test =='flat':
        test_file       =dic["flat file"]
    order               = int(dic["order"])
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    pix                 = float(dic["pix/arturo"])
    art                 = 1./pix
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
        test_data = image_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
        image_file.close()
    lenX = test_data.shape[0]
    min = np.min(test_data)
    test_data = test_data-min
    Spectrum = generate_spectrum(test_data, B_matrix)
    pickle_name = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_based_spec'+ '_OR'+str(order)+'_LA'+str(lane)+'.p'
    cPickle.dump(Spectrum, open(pickle_name, 'wb'))
    
get_spectrum('test')