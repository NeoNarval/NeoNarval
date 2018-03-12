import cPickle
import numpy.core.multiarray

path = r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Lucas_Herbert\Documents\ThAr_Line_Cal.pkl"
def read_pickle(path):
    f = open(path, 'r')
    data = cPickle.load(f)
    f.close()
    return data
    
    
print(read_pickle(path))