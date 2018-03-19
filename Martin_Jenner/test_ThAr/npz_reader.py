import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt


path = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\Bmatrix\Bmat_ord14_lane1.npz'
def npz_reader(path):
    B_matrix = sparse.load_npz(path).tolil()    
    print(B_matrix.getnnz())
    n,m =B_matrix.get_shape()
    print(n,m)
    B_matrix.toarray()
    plt.matshow(B_matrix)
    plt.show()
  
npz_reader(path)   
    