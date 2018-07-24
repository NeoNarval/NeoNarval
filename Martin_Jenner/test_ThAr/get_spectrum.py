import pyfits
import cPickle
import numpy as np
import collections
from scipy import sparse
import matplotlib.pyplot as plt



# Method used to get the closest integer of a real number
def rint(x): return int(round(x))

def ichol(A):
    G = np.copy(A)
    n =np.shape(G)[0]
    for k in range(0,n):
        G[k,k] = np.sqrt(G[k,k])
        for i in range(k+1,n):
            if (G[i,k] != 0):
                G[i,k] = G[i,k]/G[k,k]            
            

        for j in range(k+1,n):
            for i in range(j,n):
                if (G[i,j]!=0):
                    G[i,j] = G[i,j]-G[i,k]*G[j,k]  

    for i in range(0,n):
        for j in range(i+1,n):
            G[i,j] = 0
    return(G)
    
def conjugate_grad_prec(A,Y):
    iter = 100
    plt.figure(200)

    # ========== Precondtioned conjugate gratient method ==========
    # https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
    # Pas besoin de preconditionner donc on peut directement utiliser C1 et b1
    n = np.shape(Y)[0]
    sig = sparse.spdiags(np.sqrt(Y.T), 0, n,n) 
    C=((A.T).dot(sig)).dot(A)
    b=A.T.dot(Y)
    
    # Calcul de la matrice M de preconditionnement
    
    G = ichol(C)
    
    M = np.dot(G.T,G)
    DM = np.diag(M)
    #DM=np.zeros_like(M) 
    #np.fill_diagonal(DM,1./np.diag(M))
    np.fill_diagonal(M,DM)
    
    # Initialisation :
    x0 = np.zeros((np.shape(C)[0],1))
    r0 = b - np.dot(C,x0)
    z0 = np.dot(M,r0)
    p0 = r0
    k = 0
    
    xk = x0
    pk = p0
    rk = r0
    zk = z0
    print(np.shape(xk),np.shape(pk),np.shape(rk),np.shape(zk))
    print(np.shape(C.dot(pk)))
    # Iteration
    for i in range(iter):
        
        ak = float(np.dot(zk.T,rk) / np.dot(pk.T, C.dot(pk) ))
        new_xk = xk + ak*pk
        new_rk = rk - ak*np.dot(C,pk)

        new_zk = np.dot(M.T,new_rk)
        bk = float(np.dot(new_zk.T,new_rk) / np.dot(zk.T,rk))
        new_pk = new_zk + bk*pk

        xk = new_xk
        rk = new_rk
        pk = new_pk
        zk = new_zk

    # Plotting the results
    new_xk = new_xk/np.max(new_xk)		
    plt.plot(new_xk,color='blue')
    plt.show()

    return new_xk


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
    cPickle.dump(CCD, open(r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\second_membre_CCD', 'wb'))
    print(np.shape(CCD))
    
    
 # # We use the method of Normal Equations to retrieve the spectrum of the lane
    n = np.shape(CCD)[0]
    sig = sparse.spdiags(np.sqrt(CCD.T), 0, n,n)      
    C = (B_matrix.T).dot(sig)
    C = C.dot(B_matrix)
    C = C.toarray()
    print(np.linalg.cond(C))
    k = ((B_matrix.T).dot(sig)).dot(CCD)
    
#--------------------linalg_solve------------------
    # x_lin = np.linalg.solve(C,k)
    # 
    # # plt.figure()
    # # plt.plot(x_lin)
    # # plt.show()
    # spectrum = (list(np.array(x_lin.T.tolist()[0])))
# # #--------------------SSOR Method-----------------------
    omega = 1
    iter = 50
    diag = np.diag(C)
    D = np.zeros_like(C)
    inv_D  = np.zeros_like(C)
    np.fill_diagonal(D, diag)
    np.fill_diagonal(inv_D, 1./diag)
    
    L = np.dot(inv_D,-np.tril(C))
    U = np.dot(inv_D,-np.triu(C))
    I = np.identity(C.shape[0])
    
    
    x_simple = inv_D.dot(k)
    # plt.close(100)
    # plt.figure(100)
    # plt.plot(x_simple)
    # plt.show()
    
    Alpha = I-omega*L
    Beta = (1-omega)*I + omega*U
    #x_old = np.linalg.solve(C,k)
    
    #print(np.linalg.cond(Alpha))
    
    x_old = x_simple
    
    for i in range(iter):
        x = np.linalg.solve(Alpha, Beta.dot(x_old)+ omega*inv_D.dot(k))
        #x = np.linalg.pinv(Alpha).dot(Beta).dot(x_old) + omega
        x_old = x
        #plt.cla()
    plt.figure()
    plt.plot(x)
    plt.title("SSOR iteration {0}".format(iter))
    plt.show()

    spectrum = (list(np.array(x.T.tolist()[0])))
#-----------------------------------------------------

#---------------Lucas Method--------------------
    # spectrum = (list(np.array(conjugate_grad_prec(B_matrix.toarray(),CCD).T.tolist()[0])))
#-----------------------------------------------

#-------------------Cholesky -simple- Method----------
    # G = np.linalg.cholesky(C) # Cholesky factorization
    # #plt.matshow(G, aspect ='auto')
    # #plt.matshow(np.linalg.inv(C)*B_matrix.T, aspect='auto')
    # #plt.show()
    # y = np.linalg.solve(G, k)
    # X = np.linalg.solve(G.T, y)
    # spectrum = (list(np.array(X.T.tolist()[0])))
#-----------------------------------------------------

# #---------------Landweber Method----------------------
#  for i in range(iter):
#     x = np.linalg.solve(Alpha, Beta.dot(x_old)+ omega*inv_D.dot(k))
#     #x = np.linalg.pinv(Alpha).dot(Beta).dot(x_old) + omega
#     x_old = x
#-----------------------------------------------------

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
    elif test == 'star':
        test_file       =dic["star file"]
    flat_file           =dic["flat file"]
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
    
    image_f_file = pyfits.open(flat_file)
    flat_data = image_f_file[0].data.astype(np.float32) # Data of the ccd of the star fts file
    image_f_file.close()
    
    lenX = test_data.shape[0]
    min = np.min(test_data)
    test_data = test_data-min
    #------------flat-field handling---------
    Spectrum = generate_spectrum(test_data, B_matrix)
    flat_field = generate_spectrum(flat_data, B_matrix)
    for i in range(len(Spectrum)):
        Spectrum[i] = Spectrum[i]/flat_field[i]
    plt.figure()
    plt.plot(Spectrum)
    plt.show()
    #-----------------------------------------
    pickle_name = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_based_spec'+ '_OR'+str(order)+'_LA'+str(lane)+'.p'
    cPickle.dump(Spectrum, open(pickle_name, 'wb'))
    
get_spectrum('ThAr')