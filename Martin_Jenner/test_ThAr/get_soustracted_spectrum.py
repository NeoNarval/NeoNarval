import pyfits
import cPickle
import numpy as np
import collections
from scipy import sparse
import matplotlib.pyplot as plt

# Method used to get the closest integer of a real number
def rint(x): return int(round(x))

def generate_spectrum(test_data, B_matrix):
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
    
def get_soustracted_spectrum(data_file, matrix_file):
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
    ThAr_file           = dic["ThAr fts file"]
    FP_file             = dic["FP fts file"]
    test_file          = dic["test file"]
    flat_file           = dic["flat file"]
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
    
    test_data = cPickle.load(open(test_file, 'r'))
    
    image_file = pyfits.open(ThAr_file)
    ThAr_data = image_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
    image_file.close()
    
    image_file = pyfits.open(flat_file)
    flat_data = image_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
    image_file.close()
    
    image_file = pyfits.open(FP_file)
    FP_data = image_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
    image_file.close()
    
    
    if data_file == 'ThAr':
        Spectrum = generate_spectrum(ThAr_data, B_matrix)
    elif data_file =='FP':
        Spectrum = generate_spectrum(FP_data, B_matrix)
    elif data_file =='test':
        Spectrum = generate_spectrum(test_data, B_matrix)
    elif data_file =='flat':
        Spectrum = generate_spectrum(flat_data, B_matrix)
        
    
    if matrix_file =='ThAr':
        inf_envelope = generate_spectrum(ThAr_data, B_matrix)
    elif matrix_file =='FP':
        inf_envelope = generate_spectrum(FP_data, B_matrix)
    elif matrix_file =='test':
        inf_envelope = generate_spectrum(test_data, B_matrix)
    elif matrix_file =='flat':
        inf_envelope = generate_spectrum(flat_data, B_matrix)
        
    s=0
    k=0
    s_lim = np.mean(Spectrum)
    print(s_lim)
    while s<=s_lim:
        s = Spectrum[k]
        k+=1
    print(s)
    e=0
    k=0
    e_lim = np.mean(inf_envelope)
    print(e_lim)
    while e<=e_lim:
        e = inf_envelope[k]
        k+=1
    print(e)
    print(inf_envelope)
    # plt.figure(101)
    # plt.plot(Spectrum/s)
    # plt.show()
    # plt.figure(102)
    # plt.plot(inf_envelope/e)
    # plt.show()
    sous_Spectrum = [(Spectrum[i]/s)-(inf_envelope[i]/e) for i in range(len(Spectrum))]
    pickle_name = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_based_spec'+ '_OR'+str(order)+'_LA'+str(lane)+'.p'
    cPickle.dump(sous_Spectrum, open(pickle_name, 'wb'))

#argument : ('data_file', 'matrix_file')
get_soustracted_spectrum('ThAr','flat' )