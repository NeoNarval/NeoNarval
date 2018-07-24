import cPickle
import numpy.core.multiarray
import matplotlib.pyplot as plt
import numpy as np
import collections


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
path = r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\TEMP_\ThAr_based_spec_OR"+order+"_LA"+lane+".p"
print(path)
def read_pickle(path):
    f = open(path, 'r')
    data = cPickle.load(f)
    f.close()
    abs = np.linspace(lambd_ini,lambd_end,len(data))
    # print(len(data))
    # print(len(abs))
    plt.figure()
    plt.plot(abs, data, color = 'green')
    plt.xlabel("Wavelength in arturo")
    plt.show()
    
    return data
def read_data():
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
    n_order_pos = dic["n order pos"]
    n_order_thic = dic["n order thic"]
    n_polys_env = dic["n polys env"]
    lanes_pos_data = cPickle.load(open(n_lanes_pos, 'r'))
    lanes_thic_data = cPickle.load(open(n_lanes_thic, 'r'))
    order_pos_data = cPickle.load(open(n_order_pos, 'r'))
    order_thic_data = cPickle.load(open(n_order_thic, 'r'))
    polys_env_data = cPickle.load(open(n_polys_env, 'r'))
    #-----------------------------------
    
    print("arturo order ref Profil : ", np.shape(order_ref_data['Profil']))
    # plt.figure()
    # plt.plot(list(order_ref_data['Orders']),'o')
    # plt.show()
    # print("arturo order ref Beams : ", order_ref_data['Beams'])
    # print("arturo order ref Order : ", order_ref_data['Orders'])
    # print("arturo order curv: ", order_curv_data['Coefs'])
    # plt.figure()
    # plt.plot(order_curv_data['map'])
    # plt.show()
    # print("arturo order ref : ", order_ref_data)
    # print(np.shape(order_curv_data['map']))
    
    # print("polynome :",polys_data)
    # print("lanes thickness :",thickness_data)
    # print("length lanes thickness :", len(thickness_data))
    # print("enveloppe: ", envelope_data)
    # print(np.shape(envelope_data))
    # 
    # print('lanes pos : ', lanes_pos_data)
    # plt.figure()
    # plt.plot(lanes_pos_data)
    # plt.show()
    # print('lanes thic : ',lanes_thic_data)
    # print('order pos : ', order_pos_data)
    # print('order thic : ', order_thic_data)
    print('polys env : ', polys_env_data)
read_data()