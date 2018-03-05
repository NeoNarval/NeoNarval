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
EWA=[]
EWN=[]
Res=[]
puntos=[]


fin=open("Lambda_Cal.pkl",'r')
Info=pickle.load(fin)
#{'Order':order,'Lines':puntos,'Error':error2}                                  
fin.close()
order0=Info['Order']+1
puntos=Info['Lines']
error2=Info['Error']
		


for order in range(order0,37):
	
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
	
	
	for i in np.arange(0,len(Owv)-300,300):
	
		d0=i
		d1=i+300
	
		cont=np.amax(Oprof[d0:d1])
		dispersion=(Owv[d1]-Owv[d0])/50.
		cuales=[]
		flag=True
		while flag:
			fig=plt.plot(Owv[d0:d1],Oprof[d0:d1])
			#plt.plot(Owv[d0:d1],p(Owv[d0:d1]))
			#plt.xlim(Owv[d0],Owv[d1])
			#for este in cuales:
			#	plt.axvline(este)
			
			c=Cursor(fig,Owv[d0:d1],Oprof[d0:d1],tipo='V')
			plt.show(fig)
			if c.xdata>0:
				cuales.append(c.xdata)
			else:
				flag=False
		plt.close()
		plt.ioff()
		plt.subplot(2,1,1)
		fig=plt.plot(Owv[d0:d1],Oprof[d0:d1],'b')
		plt.plot(Owv[d0:d1],p(Owv[d0:d1]),'g')	
		plt.title('Fitting...')
		for cual in cuales:
			aqui=np.argmin(abs(Owv-cual))
		
			#try:
			ext=7
			cont=np.amax(Oprof[aqui-ext:aqui+ext])
			if cont>0:
				p0 = [1., 0., 1.,cont]
				coeff, var_matrix = curve_fit(Voigt,np.arange(-ext,ext),Oprof[aqui-ext:aqui+ext], p0=p0,bounds=([0.1,-3,0.1,0.],[1.,3,3.,1.5]))
				aquiA=np.argmin(abs(watlas-cual))
				dispersion=(Owv[aqui+10]-Owv[aqui-10])/21.
				p0 = [1., 0., 1.,1.]
				coeffA, var_matrix = curve_fit(Voigt,np.arange(-20,20),Fatlas[aquiA-20:aquiA+20], p0=p0)
				dispAtlas=(watlas[aquiA+50]-watlas[aquiA])/50.
				DeltaS=dispersion*coeff[2]-dispAtlas*coeffA[2]
				#print(dispersion*coeff[2],dispAtlas*coeffA[2])
				if abs(coeff[1])<ext and DeltaS>0 and abs(DeltaS)<2:
					#plt.plot(watlas,Fatlas)	
					plt.plot(np.linspace(Owv[aqui-20],Owv[aqui+20],100),Voigt(np.linspace(-20,20,100),*coeff),color='r')
					
			
					plt.plot(watlas[aquiA-20:aquiA+20],Voigt(np.arange(-20,20),*coeffA),color='r')
					plt.axvline(watlas[aquiA]+dispAtlas*coeffA[1],color='r',linestyle='--')
					plt.axvline(Owv[aqui]+dispersion*coeff[1],color='g',linestyle='--')
					escala=1.#1000.*3e5/i
					error1.append((cual-(Owv[aqui]+dispersion*coeff[1]))*escala)
					error2.append((watlas[aquiA]+dispAtlas*coeffA[1]-(Owv[aqui]+dispersion*coeff[1]))*escala)
					errorD.append(DeltaS)
					EWN.append(coeff[0]*coeff[2]*dispersion*np.sqrt(2.*np.pi))
					EWA.append(coeffA[0]*coeffA[2]*dispAtlas*np.sqrt(2.*np.pi))
					#errorEW.append(EW-EWA)
					puntos.append(watlas[aquiA]+dispAtlas*coeffA[1])
					Res.append((dispersion*coeff[2])**2-(dispAtlas*coeffA[2])**2)
				else:
					print(coeff[1],DeltaS)
			#except:
			#	pass	
		plt.subplot(2,1,2)
		plt.plot(puntos,error2,'.')
		plt.title(r' {0:0.3f} $\pm$ {1:0.3f}'.format(np.mean(error2),np.std(error2)))
		plt.show()
	Info={'Order':order,'Lines':puntos,'Error':error2}                                  
	fout=open("Lambda_Cal.pkl",'w')
	pickle.dump(Info,fout)
	fout.close()
		
	# plt.cla()
# 	plt.plot(puntos,error1,'b.')
# 	plt.title(r'Position error: {0:.3f}$\pm$ {1:.3f} $\AA$'.format(dispersion*np.mean(error1),dispersion*np.std(error1)))
# 	plt.ylabel('Error (A)')
# 	
# 	
# 	plt.draw()
# 	
	
	

