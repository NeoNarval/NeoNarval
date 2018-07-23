# -*- coding: utf-8 -*-

"""

"""

import os
import gc
import utils
import pyfits
import tables
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as interp
from scipy.signal import medfilt


def planck(wav, T):
    wav = wav*10**(-3)  #Wavelength given in mmm then converted into m
    h = 6.626e-34
    c = 3.0e+8
    k = 1.38e-23
    a = 2.0*h*c**2
    b = h*c/(wav*k*T)
    intensity = a/ ( (wav**5) * (np.exp(b) - 1.0) ) /2004623760881.600098
    return intensity

def readQE(wv):
	f=open('/home/main/Documents/DRS/DATA/NARVAL_QE.txt','r')
	f.readline()
	wvl=[]
	QE=[]
	for line in f:
	    wvl.append(float(line.split()[0]))
	    QE.append(float(line.split()[1]))
	f.close()
	wvl=np.asarray(wvl)
	QE=np.asarray(QE)

	QEf=np.zeros(len(wv))
	for i,w in enumerate(wv):
	    idx=(np.abs(wvl-w)).argmin()
	    idx0=np.amax([idx-4,0])
	    idx1=np.amin([idx+4,len(wvl)-1])
	    q=interp.UnivariateSpline(wvl[idx0:idx1],QE[idx0:idx1])
	    QEf[i]=q(w)
	QEf=QEf/QEf.max()
	QEf=np.clip(QEf,0.03,1)
	return QEf/QEf.max()

file_Lune='/home/main/Documents/DRS/DATA/th_calibered.fits'
l=pyfits.open(file_Lune)
a=l[1].data
l.close()
QEf=readQE(a['wavelength_lane1']/10.)
f1=planck(a['wavelength_lane1']/1e7,6000)

f2=f1+planck(a['wavelength_lane1']/1e7,3400)
f2=f2/f2.max()
esc=np.amax(a['intensity_lane1'])
plt.plot(a['wavelength_lane1'],a['intensity_lane1']/esc,label='Moon')
plt.plot(a['wavelength_lane1'],QEf,label='QE CCD')
plt.plot(a['wavelength_lane1'],f2,label='Lamp spectrum')
plt.plot(a['wavelength_lane1'],f1)
#plt.plot(a['wavelength_lane1'],a['intensity_lane1'])
plt.legend()
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Normalized Intensity')
plt.show()


fS=planck(a['wavelength_lane1']/1e7,5777.)
#med=medfilt(a['intensity_lane1']/esc,2501)







