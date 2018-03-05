#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy.ma as ma
import pickle
#from Chelli_Module import *
import scipy.interpolate as interp
from scipy.optimize import curve_fit

def Voigt(x, *p):
    A, mu, sigma = p
    return 1-A*np.exp(-(x-mu)**2/(2.*sigma**2))


# On lit l'atlas solaire sur a
atlas="/Users/arturo/DeepStokes/SolarSpectrum.pkl"
f=open(atlas,'r')
a=pickle.load(f)
f.close()

Nwvl=len(a['Wavelength'])
disp=(np.amax(a['Wavelength'])-np.amin(a['Wavelength']))/Nwvl
cuantos=np.fix(20/disp).astype(int)  #cogemos unos 20A cada vez


rangos=np.arange(0,Nwvl,cuantos)
lineas=[]
plt.ion()
for i in rangos[100:202]:
	
	sp=a['Intensity'][i:i+cuantos]
	der=np.gradient(sp)
	tilts=np.roll(der,1)*der
	donde=np.where(tilts<0)
	#aqui=i+np.argmin(a['Intensity'][i:i+cuantos])
	#lineas.append(aqui)
	for aqui in donde[0][1:-2]:
		
		if (der[aqui]<der[aqui+1])  and (sp[aqui]<0.99):
			plt.cla()
			sp1=a['Intensity'][i+aqui-10:i+aqui+10]
			der2=np.gradient(np.gradient(sp1))
			tilts2=np.roll(der2,1)*der2
			donde2=np.where(tilts2<0)[0]
			if len(donde2)>1:
				dist=np.argsort((donde2-10)**2)
				d1=donde2[dist[0]]
				d2=donde2[dist[1]]
				Perf0=a['Intensity'][i+aqui]
				Perf1=a['Intensity'][i+aqui-10+d1]
				Perf2=a['Intensity'][i+aqui-10+d1]
				ProbableAmp=2*max([Perf1-Perf0,Perf2-Perf0])
				if d1<10 and d2>10 and ProbableAmp>0.05:
			
					p0 = [1., 0., 1.]

					coeff, var_matrix = curve_fit(Voigt, np.arange(-10+d1,-10+d2+1), sp1[d1:d2+1], p0=p0)
					dispersion=np.mean(np.gradient(a['Wavelength'][i+aqui-10:i+aqui+10]))
					lineas.append(a['Wavelength'][i+aqui]+coeff[1]*dispersion)
					
					
					
					plt.plot(a['Wavelength'][i+aqui-30:i+aqui+30],a['Intensity'][i+aqui-30:i+aqui+30])
					plt.plot(a['Wavelength'][i+aqui-30:i+aqui+30],Voigt(np.arange(-30,30),*coeff))
					plt.axvline(a['Wavelength'][i+aqui],linestyle='--',color='b')
					plt.axvline(a['Wavelength'][i+aqui-10+d1],linestyle='--',color='r')
					plt.axvline(a['Wavelength'][i+aqui-10+d2],linestyle='--',color='r')
					plt.axvline(a['Wavelength'][i+aqui]+coeff[1]*dispersion)
					plt.draw()
					
			
stop