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

def Voigt(x, *p):
    A, mu, sigma,cont = p
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



wv=np.arange(50)
Basis=np.zeros((50,500*20))
MemoriaWidth=np.zeros(500*20)
MemoriaLambda=np.zeros(500*20)
count=0
for basisWidth in np.linspace(1,10,20):
  for C in np.linspace(0,50,500):
     Basis[:,count]=np.exp(-(wv-C)**2 / basisWidth**2)
     MemoriaWidth[count]=basisWidth
     MemoriaLambda[count]=C             
     count=count+1      
print Basis.shape



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



Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))


#Ca IR: ordre 26=20+6
# O2 770nm: ordre 30=20+10
plt.ion()

for order in range(1,31):
	aquiD=[]
	aquiA=[]
	error=[]
	print('Order {0}'.format(20+order))
	w0=a['wavelength_lane1'][Orderlimit[order]]
	w1=a['wavelength_lane1'][Orderlimit[order-1]+1]
	
	
	
	donde=ma.masked_inside(atlas["Wavelength"],w0,w1)
	watlas=np.asarray(atlas["Wavelength"])[donde.mask]
	Fatlas=np.asarray(atlas["Intensity"])[donde.mask]
	
	
	dispersion=(w1-w0)/(Orderlimit[order-1]-Orderlimit[order]) #Angstroms/pixel
	Oprof=a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	Owv=a['wavelength_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]

	donde=np.where(Oprof>0)
	# plt.subplot(2,1,1)
# 	plt.plot(Owv,Oprof)
# 	plt.ylim(0,1)
# 	plt.subplot(2,1,2)
# 	plt.plot(watlas,Fatlas)
# 	plt.ylim(0,1)
	
	for i in np.arange(donde[0][0]+50,donde[0][-1]-100,50):
	
		d0=i
		d1=i+50
	
		cont=np.amax(Oprof[d0:d1])
		dispersion=(Owv[d1]-Owv[d0])/50.
		pL =rvm(Basis,Oprof[d0:d1]-cont,noise=0.01)
		pL.iterateUntilConvergence()
		
		for j in pL.relevant:
			if (pL.wInferred[j]<-0.1):
				
				#plt.plot(Owv[d0:d1],cont+pL.wInferred[j]*Basis[:,j],color='r')
				aquiD.append(MemoriaLambda[j]*dispersion+Owv[d0])
				# plt.subplot(2,1,1)
# 				plt.axvline(aquiD[-1],color='r')
# 				
# 				#plt.text(aquiD[-1],0.5,'{0}'.format(aqui))
# 				plt.subplot(2,1,2)
# 				plt.axvline(aquiD[-1],color='r')
# 				plt.draw()
			
			
	
	for i in np.arange(0,len(watlas)-50,50):
		prof=Fatlas[i:i+50]-1
		dispersion=(watlas[i+50]-watlas[i])/50.
		pL =rvm(Basis,prof,noise=0.01)
		pL.iterateUntilConvergence()
		for j in pL.relevant:
			if (pL.wInferred[j]<-0.1):
				
				#plt.plot(watlas[i:i+50],1.+pL.wInferred[j]*Basis[:,j],color='r')
				aquiA.append(MemoriaLambda[j]*dispersion+watlas[i])
				#plt.text(aqui,0.5,'{0}'.format(aqui))
				# plt.subplot(2,1,1)
# 				plt.axvline(aquiA[-1],color='g')
# 				
# 				#plt.text(aquiA[-1],0.5,'{0}'.format(aqui))
# 				plt.subplot(2,1,2)
# 				plt.axvline(aquiA[-1],color='g')
# 				plt.draw()
	
	
	for w in aquiD:
		error.append(np.amin(abs(w-np.asarray(aquiA))))
	
	plt.plot(aquiD,error,'.')
	plt.xlim(4000,11000)
	plt.ylim(0,0.5)
	print('Ordre {0}, lambda {1}, moyenne {2}'.format(20+i,w0,np.mean(error)))
	plt.draw()		
stop		
	
	
Cont={'C0':C0,'C1':C1,'C2':C2}					
fout=open("Continuum.pkl",'w')
pickle.dump(Cont,fout)
fout.close()



col1=pyfits.Column(name="wavelength_lane1",format='E',array=a['wavelength_lane1'])
col2=pyfits.Column(name="intensity_lane1",format='E',array=a['intensity_lane1'])

cols=pyfits.ColDefs([col1,col2])
tbhdu2=pyfits.BinTableHDU.from_columns(cols)



tbhdu1=pyfits.PrimaryHDU(header=h0)

final_hdu=pyfits.HDUList([tbhdu1,tbhdu2])
final_hdu.writeto('luna_normalised.fts',clobber=True)			
stop	
	
