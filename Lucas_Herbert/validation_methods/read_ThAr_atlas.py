#!/usr/bin/env python
# -*- coding: utf-8 -*-


import pyfits
import numpy as np
import cPickle
import pickle

""" 
This function reads the file containing all the wavelengths of the spikes of the Thorium Argon spectrum and returns a python list of those wavelengths. We will use it as a reference to compare our computed wavelengths.
Input :
None.
Output : 
- wavelengths : list of the wavelengths contained in the atlas (in Angqtroms)
"""


def read_ThAr_Atlas():

    path = '/home/stagiaire/Documents/Données utiles/thar_UVES_MM090311.dat'
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
    
    file = open('/home/stagiaire//Documents/Données utiles/ThAr_Atlas.pkl','w')
    ThAr_pickle = pickle.dump(wavelengths,file)
    file.close()
    return(wavelengths)


    

