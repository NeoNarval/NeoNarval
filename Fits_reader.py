import pyfits
import matplotlib.pyplot as plt
from os import chdir

"""
read .fits file encoding a spectrum : 1 column for the wavelength another for the intensity

    : path : string with the relative path toward the file
"""
path =r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\Dépot_Git\DRS\DATA\th_calibered.fits"
def fits_reader(path):
    l = pyfits.open(path)
    a = l[1].data
    l.close()
    lambd = a['wavelength_lane1']
    intensity = a['intensity_lane1']
    plt.plot(lambd, intensity)
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Intensity')
    plt.show()
    
fits_reader(path)