import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit
import collections

plt.close(50)
plt.close(51)
"""
read .fits file encoding a spectrum : 1 column for the wavelength another for the intensity

    : path : string with the relative path toward the file
"""
path =r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\th_calibered.fits"
path_th=r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18\Narval_20180313_181059_th0.fts"
path_fp=r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18\Narval_20180313_181342_fp0.fts"
path_fla=r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18\Narval_20180313_180054_10f.fts"

docpath = r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt"
dic = collections.OrderedDict()
with open(docpath, 'r') as file:
    content = file.readlines()
content = [line.strip() for line in content]
for line in content:
    param = line.split(" : ")
    if len(param) == 2:
        nom_param = param[0]
        value_param = param[1]
        dic[nom_param] = value_param
pix = float(dic["pix/arturo"])
lambd_ini = int(dic["initial index"])*pix
lambd_end = int(dic["final index"])*pix
order = dic["order"]
lane = dic["lane"]

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
    
    plt.figure(50)
    plt.imshow(ThAr_data, aspect = 'auto')
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



def check_pos_order():
    #plot_mat(path_fla)
    
    
    #----------fichiers arthur--------
    thickness_data_file = dic["Lane thickness file"]
    polys_data_file     = dic["Polys envelope file"]
    envelope_data_file  = dic["Lane envelope file"]
    polys_data     = cPickle.load(open(polys_data_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    envelope_data  = cPickle.load(open(envelope_data_file, 'r'))
    #---------------------------------
    #---------fichiers arturo---------
    order_ref = dic["order ref"]
    order_curv = dic["order curv"]
    order_ref_data = cPickle.load(open(order_ref, 'r'))
    order_curv_data = cPickle.load(open(order_curv, 'r'))
    #----------------------------------
    #---------fichiers recalc arthur---
    n_lanes_pos = dic["n lanes pos"]
    n_lanes_thic = dic["n lanes thic"]
    
    n_polys_env = dic["n polys env"]
    lanes_pos_data = cPickle.load(open(n_lanes_pos, 'r'))
    lanes_thic_data = cPickle.load(open(n_lanes_thic, 'r'))
    polys_env_data = cPickle.load(open(n_polys_env, 'r'))
    #----------------------------------
    # plt.figure(50)
    # c = np.array(order_ref_data['Profil'])
    # plt.plot(c/np.max(c))
    # plt.show()
    
 
    plt.figure(50)
    plt.plot(2000-lanes_pos_data)
    plt.show()

    l = pyfits.open(path_fla)
    a = l[0].data
    l.close()
    plt.figure(50)
    plt.imshow(np.flip(a.T,0), aspect = 'auto')
    plt.show()
    
plot_mat(path_th)







