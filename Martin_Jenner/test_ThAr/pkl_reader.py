import cPickle
import numpy.core.multiarray
import matplotlib.pyplot as plt

path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\ThAr_based_spec_OR14_LA1.p"
lambd_ini = 0
lambd_end = 1000
order = 14
def read_pickle(path):
    f = open(path, 'r')
    data = cPickle.load(f)
    f.close()
    abs = np.linspace(lambd_ini,lambd_end,len(data))
    # print(len(data))
    # print(len(abs))
    plt.figure
    plt.plot(abs, data)
    plt.xlabel("Wavelength in arturo")
    plt.show()
    
    return data
    
    
print(read_pickle(path))