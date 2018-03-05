#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy.ma as ma
import pyfits
from cursor import *
from rvm import *
from scipy import interpolate
from Chelli_Module import *
import scipy.interpolate as interp
from scipy.optimize import curve_fit

def Voigt(x, A, mu, sigma,cont):
    #A, mu, sigma,cont = p
    return cont-A*np.exp(-(x-mu)**2/(2.*sigma**2))





# On lit l'atlas solaire sur a
atlas="/Users/arturo/DeepStokes/SolarSpectrum.pkl"
f=open(atlas,'r')
atlas=pickle.load(f)
f.close()

#Leemos la observacion de Narval
file_Lune='luna_normalised.fts'
l=pyfits.open(file_Lune)
a=l[1].data
h0=l[0].header
h1=l[1].header
l.close()

#Leemos la mascara
file_mask='G2_ekdra'
f=open(file_mask,'r')
cuantos=int(f.readline())
lineas=np.zeros(cuantos)
for i,line in enumerate(f.readlines()):
	lineas[i]=10.*float(line.split()[0])
f.close()

Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))


#Ca IR: ordre 26=20+6
# O2 770nm: ordre 30=20+10

plt.ion()
error1=[]
error2=[]
errorD=[]
Res=[]
puntos=[]
for order in range(2,37):
	
	print('Order {0}'.format(20+order))
	w0=a['wavelength_lane1'][Orderlimit[order]]
	w1=a['wavelength_lane1'][Orderlimit[order-1]+1]
	
	
	
	donde=ma.masked_inside(atlas["Wavelength"],w0,w1)
	watlas=np.asarray(atlas["Wavelength"])[donde.mask]
	Fatlas=np.asarray(atlas["Intensity"])[donde.mask]
	p=interpolate.InterpolatedUnivariateSpline(watlas,Fatlas)
	
	
	
	
	
	dispersion=(w1-w0)/(Orderlimit[order-1]-Orderlimit[order]) #Angstroms/pixel
	Oprof=a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	Owv=a['wavelength_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	
	donde=np.where(Oprof>0)
	for i in np.arange(donde[0][0]+50,donde[0][-1]-100,10):
	
		d0=i
		d1=i+100
	
		cont=np.amax(Oprof[d0:d1])
		dispersion=(Owv[d1]-Owv[d0])/100.
		
		
		
		Atprof=p(Owv[d0:d1])
		
		
		v=Chelly(Atprof,Oprof[d0:d1]/cont,v0=0.1)
		
		
		# plt.plot(Owv[d0:d1],Oprof[d0:d1])
# 		plt.plot(Owv[d0:d1],Atprof)
# 		plt.draw()
		error1.append(dispersion*v)
		puntos.append(Owv[d0])
	plt.cla()
	plt.plot(puntos,error1,'b.')
	plt.title(r'Position error: {0:.3f}$\pm$ {1:.3f} $\AA$'.format(np.mean(error1),np.std(error1)))
	plt.ylabel('Error (A)')
	
	
	plt.draw()
	
	
	

