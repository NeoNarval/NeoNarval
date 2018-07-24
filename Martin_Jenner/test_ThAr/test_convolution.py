import pyfits
import cPickle
import numpy as np
import collections
from scipy import sparse
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
test = 'ThAr'

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
        
#----------load fts file----------------- not used yet
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
#-----------------------------------------------------
    
order = int(dic['order'])
ini_index = int(dic['initial index'])
end_index = int(dic['final index'])
B_matrix_path       = dic["B matrix"]
B_matrix            = sparse.load_npz(B_matrix_path).tocsr()

spectrum_file = r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\DRS\DATA\th_calibered_updated.fits"
image_file = pyfits.open(spectrum_file)
spectrum = image_file[1].data['intensity_lane1']
wavelength = image_file[1].data['wavelength_lane1']
image_file.close()

Y_CCD = cPickle.load(open(r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\second_membre_CCD", 'r'))

# spectrum = np.flip(spectrum,0)
# spectrum = spectrum[order*4612 + ini_index:order*4612 + end_index]
# plt.figure(1)
# # plt.plot(wavelength[-((order-1)*4612+end_index):-((order-1)*4612+ini_index)],spectrum[-((order-1)*4612+end_index):-((order-1)*4612+ini_index)])
# plt.plot(spectrum)
# plt.show()
# 
# Y = B_matrix.dot(spectrum)
# 
# plt.figure(2)
# plt.plot(Y)
# plt.show()
# 
# plt.figure(3)
# plt.plot(Y_CCD)
# plt.show()
# 
# C = (B_matrix.T).dot(B_matrix)
# S = spsolve(C,(B_matrix.T).dot(Y))
# 
# plt.figure(1)
# plt.plot(S)
# plt.show()

X = (B_matrix.T).dot(Y_CCD)
plt.figure(1)
plt.plot(X)
plt.show()

