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



wv=np.arange(100)
Basis=np.zeros((100,300*20))
MemoriaWidth=np.zeros(300*20)
MemoriaLambda=np.zeros(300*20)
count=0
for basisWidth in np.linspace(1,10,20):
  for C in np.linspace(0,100,300):
     Basis[:,count]=np.exp(-(wv-C)**2 / basisWidth**2)
     MemoriaWidth[count]=basisWidth
     MemoriaLambda[count]=C             
     count=count+1      
print Basis.shape



# On lit l'atlas solaire sur a
# atlas="/Users/arturo/DeepStokes/SolarSpectrum.pkl"
# f=open(atlas,'r')
# atlas=pickle.load(f)
# f.close()

#Leemos la observacion de Narval
file_Lune='luna_calibered.fits'
l=pyfits.open(file_Lune)
a=l[1].data
h0=l[0].header
h1=l[1].header
l.close()

Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))

C0=[]
C1=[]
C2=[]
#Ca IR: ordre 26=20+6
# O2 770nm: ordre 30=20+10
plt.ion()
for order in range(1,36):
	print('Order {0}'.format(20+order))
	w0=a['wavelength_lane1'][Orderlimit[order]]
	w1=a['wavelength_lane1'][Orderlimit[order-1]+1]
	
	
	dispersion=(w1-w0)/(Orderlimit[order-1]-Orderlimit[order]) #Angstroms/pixel
	Oprof=a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]
	Owv=a['wavelength_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]


	calib0=Oprof
	
	m=np.median(calib0)
	estos=np.where(calib0-m>0)
	coeff0=np.polyfit(Owv[estos],calib0[estos],1)
	C0.append(coeff0)
	print(coeff0)
	p0 = np.poly1d(coeff0)
	#plt.plot(Owv,Oprof/p(Owv))
	#plt.plot(Owv,p(Owv))
	#plt.show()
	#a['intensity_lane1'][Orderlimit[i-1]+1:Orderlimit[i]]=a['intensity_lane1'][Orderlimit[i-1]+1:Orderlimit[i]]/p(Owv)
	Oprof=Oprof/p0(Owv)
	
	donde=np.where(Oprof>0)
		
	pesos=np.ones(len(Owv))
	pesos[0:donde[0][0]+50]=0.
	bkg=p0(Owv)
	
	for i in np.arange(donde[0][0]+50,donde[0][-1]-100,100):
	
		d0=i
		d1=i+100
	
		cont=np.amax(Oprof)
		
		pL =rvm(Basis,Oprof[d0:d1]-cont,noise=0.01)
		pL.iterateUntilConvergence()
		for j in pL.relevant:
			if (abs(pL.wInferred[j])>0.1):
				w3=MemoriaWidth[j]*np.sqrt(np.log(4.))
				ini=np.fix(MemoriaLambda[j]-w3).astype(int)
				fin=np.fix(MemoriaLambda[j]+w3).astype(int)
				pesos[d0+ini:d0+fin]=0.
			
			
	cuantos=len(donde[0])-151#donde[0][0]+100,donde[0][-1]-50
	aqui=np.arange(len(Owv))-(donde[0][0]+50)
	coeff=np.polyfit(np.arange(cuantos),10+Oprof[donde[0][0]+50:donde[0][-1]-100],2,w=pesos[donde[0][0]+50:donde[0][-1]-100])
	C1.append(coeff)
	print(coeff)
	p = np.poly1d(coeff)
		
	Oprof=Oprof/(p(aqui)-10)
	bkg=bkg*(p(aqui)-10)
	
	
	lbd=np.copy(Owv)
	flx=np.copy(Oprof)
	for i in np.arange(10):
		lbd,flx=FitCont(lbd,flx)
	flin=interpolate.interp1d(lbd,flx,kind='linear',bounds_error=False,fill_value=1.0)
	
	coeff2=np.polyfit(lbd,flx,2)
	C2.append(coeff2)
	fcp=flin(Owv)
	p2 = np.poly1d(coeff2)
	
	plt.subplot(2,1,1)
	plt.title('Lines of order {0}'.format(20+order))
	plt.plot(Owv,calib0,color='blue')
	plt.plot(Owv,p0(Owv)*(p(aqui)-10)*fcp,color='red')
	plt.plot(Owv,p0(Owv)*(p(aqui)-10)*p2(Owv),color='green')
	plt.text(np.mean(Owv),np.max(calib0),'{0}'.format(20+order))
	plt.xlim(4000,11000)
	
	plt.subplot(2,1,2)
	plt.plot(Owv,Oprof/p2(Owv),color='blue')
	plt.text(np.mean(Owv),1.1,'{0}'.format(20+order))
	plt.xlim(4000,11000)
	plt.draw()
	a['intensity_lane1'][Orderlimit[order-1]+1:Orderlimit[order]]=Oprof/p2(Owv)
		
	
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
	
