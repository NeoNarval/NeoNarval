import cPickle
import numpy.core.multiarray
import matplotlib.pyplot as plt
import numpy as np
import collections


docpath = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt"
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
path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_based_spec_OR"+order+"_LA"+lane+".p"
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
    
read_pickle(path)