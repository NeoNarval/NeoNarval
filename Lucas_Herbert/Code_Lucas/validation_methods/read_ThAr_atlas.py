import pyfits
import numpy as np
import cPickle
import pickle

""" 
This function reads the file containing all the wavelengths of the spikes of the Thorium Argon spectrum and returns a python list of those wavelengths. We will use it as a reference to compare our computed spectra.
"""


def read_ThAr_atlas():

    path = '/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/thar_UVES_MM090311.dat'
    file = open(path)    
    wavelengths = []
    
    for i in range(0,3005):
        
        line = file.readline()
        line_lambda = float(line[12:23])
        wavelengths.append(line_lambda)
            
        
    for i in range(3006,3223):
        
        line = file.readline()
        line_lambda = float(line[11:23])
        wavelengths.append(line_lambda)
    
    file = open('/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/ThAr_Atlas.pkl','w')
    ThAr_pickle = pickle.dump(wavelength_slits,file)
    file.close()
    return(wavelengths)


    

