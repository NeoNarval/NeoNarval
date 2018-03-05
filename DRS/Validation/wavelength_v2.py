#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy.ma as ma
import pyfits
from cursor import *
from rvm import *
from scipy import interpolate

#from Chelli_Module import *
import scipy.interpolate as interp
from scipy.optimize import curve_fit

def Voigt(x, A, mu, sigma,cont):
    #A, mu, sigma,cont = p
    return cont-A*np.exp(-(x-mu)**2/(2.*sigma**2))

def FitCont(lbd,flx):
#
#--- input is:
#          lbd..... wavelength grid
#          flx..... initial (unproperly normalized) flux values
#
#--- returns estimate of continuum level (to be reinterpolated onto original grid for correction)
#
    pp=np.polyfit(lbd,flx,8)
    fcp=pp[0]*lbd*lbd*lbd*lbd*lbd*lbd*lbd*lbd + pp[1]*lbd*lbd*lbd*lbd*lbd*lbd*lbd + pp[2]*lbd*lbd*lbd*lbd*lbd*lbd + pp[3]*lbd*lbd*lbd*lbd*lbd + pp[4]*lbd*lbd*lbd*lbd + pp[5]*lbd*lbd*lbd + pp[6]*lbd*lbd + pp[7]*lbd + pp[8]

    dif=flx-fcp
    mmf=np.mean(dif)
    ect=np.std(dif)
    crit= (dif > mmf-0.5*ect) & (dif < (mmf+3*ect) )

    lbd=lbd[crit]
    flx=flx[crit]

    return lbd,flx



# wv=np.arange(10)
# Basis=np.zeros((10,500*20))
# MemoriaWidth=np.zeros(500*20)
# MemoriaLambda=np.zeros(500*20)
# count=0
# for basisWidth in np.linspace(1,5,20):
#   for C in np.linspace(0,10,50):
#      Basis[:,count]=np.exp(-(wv-C)**2 / basisWidth**2)
#      MemoriaWidth[count]=basisWidth
#      MemoriaLambda[count]=C             
#      count=count+1      
# print Basis.shape



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

#plt.ion()
error1=[]
error2=[]
errorD=[]
EWN=[]
EWA=[]
Res=[]
puntos=[]
for order in range(1,37):
	
	print('Order {0}'.format(20+order))
	w0=a['wavelength_lane1'][Orderlimit[order]]
	w1=a['wavelength_lane1'][Orderlimit[order-1]+1]
	
	
	
	donde=ma.masked_inside(atlas["Wavelength"],w0,w1)
	watlas=np.asarray(atlas["Wavelength"])[donde.mask]
	Fatlas=np.asarray(atlas["Intensity"])[donde.mask]
	
	cuales=ma.masked_inside(lineas,w0,w1)
	
	
	
	
	dispersion=(w1-w0)/(Orderlimit[order-1]-Orderlimit[order]) #Angstroms/pixel
	Oprof=a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	Owv=a['wavelength_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]

	plt.subplot(2,1,1)
	plt.plot(Owv,Oprof,'b')
	plt.plot(watlas,Fatlas,color='g')
	
	for i in lineas[cuales.mask]:
		
		aqui=np.argmin(abs(Owv-i))
		ext=6
		
		try:
			cont=np.amax(Oprof[aqui-ext:aqui+ext])
			if cont>0:
				p0 = [1., 0., 1.,cont]
				coeff, var_matrix = curve_fit(Voigt,np.arange(-ext,ext),Oprof[aqui-ext:aqui+ext], p0=p0,bounds=([0.1,-3,0.1,0.],[1.,3,3.,1.5]))
				dispersion=(Owv[aqui+10]-Owv[aqui-10])/20.
				aquiA=np.argmin(abs(watlas-i))
				p0 = [1., 0., 1.,1.]
				coeffA, var_matrix = curve_fit(Voigt,np.arange(-15,15),Fatlas[aquiA-15:aquiA+15], p0=p0)
				dispAtlas=(watlas[aquiA+50]-watlas[aquiA])/50.
				DeltaS=dispersion*coeff[2]-dispAtlas*coeffA[2]
				#print(dispersion*coeff[2],dispAtlas*coeffA[2])
				if abs(coeff[1])<4 and DeltaS>0 and abs(DeltaS)<2:
					#plt.plot(watlas,Fatlas)	
					plt.plot(Owv[aqui-4:aqui+4],Voigt(np.arange(-4,4),*coeff),color='r')
			
				
					plt.plot(watlas[aquiA-15:aquiA+15],Voigt(np.arange(-15,15),*coeffA),color='r')
					
					#print(i,Owv[aqui]+dispersion*coeff[1],watlas[aquiA]+dispAtlas*coeffA[1])
					escala=1.#1000.*3e5/i
					error1.append((i-(Owv[aqui]+dispersion*coeff[1]))*escala)
					error2.append((watlas[aquiA]+dispAtlas*coeffA[1]-(Owv[aqui]+dispersion*coeff[1]))*escala)
					errorD.append(DeltaS)
					EWN.append(coeff[0]*coeff[2]*dispersion*np.sqrt(2.*np.pi))
					EWA.append(coeffA[0]*coeffA[2]*dispAtlas*np.sqrt(2.*np.pi))
					#errorEW.append(EW-EWA)
					puntos.append(watlas[aquiA]+dispAtlas*coeffA[1])
					Res.append((dispersion*coeff[2])**2-(dispAtlas*coeffA[2])**2)
			#plt.axvline(i)
		except:
			pass
	#plt.cla()
	plt.subplot(2,1,2)
	
	plt.plot(puntos,errorD,'b.')
	plt.plot(puntos,error2,'r.')
	plt.ylabel('Error (m/s)')
#	plt.ylim(-0.1,0.1)
	plt.title(np.std(error2))	
	plt.draw()
#	stop


print(np.mean(Res))


plt.ioff()
plt.cla()
plt.subplot(1,2,1)
plt.plot(puntos,error2,'b,')
plt.title(r'Position error: {0:.3f} $\pm$ {1:.3f} $\AA$'.format(np.mean(error2),np.std(error2)))
plt.ylabel('Error (A)')
plt.subplot(1,2,2)
plt.plot(EWN,EWA,'r,')
plt.ylabel('Atlas Eq. Width (A)')
plt.xlabel('Narval Eq. Width (A)')
plt.plot(np.arange(0,3),np.arange(0,3),'b-')

#plt.title(r'Resolution: {0:5.0f} $\pm$ {1:5.0f} $\AA$'.format(np.mean(puntos/np.sqrt(np.mean(Res))),np.std(puntos/np.sqrt(np.mean(Res)))))
# plt.plot(puntos,errorD,'r.')
# plt.title(r'Doppler width error: {0:.3f} $\pm$ {1:.3f} $\AA$'.format(np.mean(errorD),np.std(errorD)))

#	plt.ylim(-0.1,0.1)

#plt.ylabel('Error (A)')
plt.show()
plt.plot(errorD,error2,'.')
plt.show()
#	donde=np.where(Oprof>0)
	# plt.subplot(2,1,1)
# 	plt.plot(Owv,Oprof)
# 	plt.ylim(0,1)
# 	plt.subplot(2,1,2)
# 	plt.plot(watlas,Fatlas)
# 	plt.ylim(0,1)
	

