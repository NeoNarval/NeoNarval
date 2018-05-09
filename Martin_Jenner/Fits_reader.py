import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit
import collections

"""
read .fits file encoding a spectrum : 1 column for the wavelength another for the intensity

    : path : string with the relative path toward the file
"""
path =r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\th_calibered.fits"
path_th=r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\DRS\FILES\Narval_20180313_181059_th1.fts"
path_fp=r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\DRS\FILES\Narval_20180313_181342_fp1.fts"

def order_gen(lambd, intensity):
    L = len(lambd)
    k = 1
    order_l = []
    order_i = []
    current_order_l = [lambd[0]]
    current_order_i = [intensity[0]]
    for k in range (L):
        if lambd[k-1]<lambd[k] and k<L-2:
            current_order_l.append(lambd[k])
            current_order_i.append(intensity[k])
        else :
            order_l.append(list(current_order_l))
            order_i.append(list(current_order_i))
            current_order_l = [lambd[k]]
            current_order_i = [intensity[k]]
            
    return order_l, order_i
    
def plot_all(lambd, intensity):
    plt.plot(lambd, intensity)
    plt.xlabel(r'wavelength ($\AA$)')
    plt.ylabel('intensity')
    plt.show()
    
def plot_orders(lambd, intensity, n):
    
    order_l, order_i = order_gen(lambd, intensity)
    plt.figure()
    plt.plot(order_l[n], order_i[n])
    
    plt.show()
        
def plot_mat(path):
    image_file = pyfits.open(path)
    ThAr_data = image_file[0].data # Data of the ThAr fts file
    image_file.close()
    
    plt.matshow(ThAr_data, aspect = 'auto')
    plt.show()
    

def fits_reader(path):
    l = pyfits.open(path)
    a = l[1].data
    l.close()
    lambd = a['wavelength_lane1']
    #lambd = np.linspace(500,1000,500)
    intensity = a['intensity_lane1']
    #i2=intensity[500:1000]
    return lambd, intensity
 
# lambd, intensity = fits_reader(path)
# order_l, order_i = order_gen(lambd,intensity)
# print(len(order_l))
# #print(order_i[0])
# #plot_all(lambd, intensity)
# plot_orders(lambd, intensity, 23)

plot_mat(path_fp)