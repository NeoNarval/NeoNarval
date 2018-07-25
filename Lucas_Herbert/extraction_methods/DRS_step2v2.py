#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
from extraction_methods.Chelli_Module import *
from scipy.optimize import curve_fit
from scipy import sparse
import time
import extraction_methods.gaussianfit2D as gf2D
from pylab import *
from mpl_toolkits.mplot3d import Axes3D



def Vfunction(x,a,p,d):
		gauss=a * np.exp(-(x-p)**2/d)
		return(-np.gradient(gauss)) 


workdir='./temp'
CCDsize=4640
ref=2320 
ancho=11

dir='extraction_methods/Extraction_files/13mar18/'
#dir='extraction_methods/Extraction_files/Moon/'

# Get FLAT done
files=sorted([f for f in os.listdir(dir) if f.endswith('fp0.fts')])
#f='/Users/arturo/NeoNARVAL/Espadons/1772472f.fits'
if len(files)>0:
	img=0.
	for f in files:
		a=pyfits.open(dir+f)
		img=img+a[0].data
		a.close()

	fp=img/len(files)
else:
	files=sorted([f for f in os.listdir(dir) if f.endswith('fp1.fts')])
	a=pyfits.open(dir+files[0])
	fp=a[0].data
	a.close()	

f = pyfits.open('extraction_methods/Extraction_files/temp/Narval_20180313_180054_f10.fts')
flat = f[0].data
f.close()


#fp = fp/flat


f=open('extraction_methods/Extraction_files/temp/order_reference.pkl','r')
a=pickle.load(f) 
f.close()
perfref0=a['Profil']
interorden=a['Orders']
interpol=a['Beams']

f=open("extraction_methods/Extraction_files//temp/Order_Curvature.pkl","r")
a=pickle.load(f) 
f.close()
mapa=a['map']
Norders=len(interorden)-2

if (False): 	
	
	FPmap=[]
	FPpics=[]
	for i,order in enumerate(interorden[0:Norders]):
		perfil=np.zeros(4640)
		#img=np.zeros((50,4640))
		o0=order
		o1=interorden[i+1]
		ov=interpol[i]
		for j,desp in enumerate(mapa[i,:]):
			perfil[j]=np.sum(fp[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)])
			#img[0:ov+1,j]=fp[j,o0+np.floor(desp):o0+ov+np.floor(desp)+1]
		
		der=np.gradient(perfil)
		#plt.plot(der)
		picos=[]
		for j in range(10,CCDsize-11):
			j0=max((j-50),0)
			j1=min((j+50),CCDsize-1)
			if (der[j]*der[j-1]<0) and (perfil[j]>np.mean(perfil[j0:j1])):
		
				try:
					popt, pcov = curve_fit(Vfunction, np.arange(20),der[j-10:j+10],p0=[perfil[j],10,4] )
					if (1<popt[2]<9) and (popt[0]<0):
						picos.append(j-10+popt[1])
						#plt.axvline(j-10+popt[1],linestyle='--')					
	# 					plt.plot(j-10+np.arange(20), Vfunction(np.arange(20), *popt),color='r')
				except:
					pass
		print("Order {0}, {1} peaks found".format(i,len(picos)))
		FPpics.append(picos)
		FPmap.append(perfil)
# 		plt.subplot(2,1,1)
# 		plt.title('Order {0}'.format(i))
# 		plt.plot(der)
# 		for i in picos:
# 			plt.axvline(i,linestyle='--')
# 		plt.subplot(2,1,2)
# 		x=np.arange(len(picos))
# 		c,var=np.polyfit(x,picos,3,cov=True)
# 		# plt.plot(picos,'o')
# # 		plt.plot(c[0]*x**3+c[1]*x**2+c[2]*x+c[3])
# 		plt.plot(picos-(c[0]*x**3+c[1]*x**2+c[2]*x+c[3]))
# 		plt.show()

	
		fout=open("extraction_methods/Extraction_files/temp/FP_map_new.pkl","w")
		pickle.dump({'picos':FPpics,'perfiles':FPmap,'info':"Positions of the FP peaks, order by order"},fout)
		fout.close()
else:
	fout=open("extraction_methods/Extraction_files/temp/FP_map_new.pkl","r")
	a=pickle.load(fout)
	FPpics=a['picos']
	FPmap=a['perfiles']
	fout.close()


for i in np.arange(0,Norders-1):
	t0=time.time()
	o0=interorden[i]
	o1=interorden[i+1]
	ov=interpol[i]
	
	picos=FPpics[i]
	arturos=len(FPmap[i])
	
	jini=5
	jfin=len(picos)-5
	
	
	

	#modif	
	#A=np.zeros( ((ov+1)*CCDsize,arturos))

	
	A2 = sparse.lil_matrix((arturos,(ov+1)*CCDsize)) # A2 = A.T
	
	for j in np.arange(jini,jfin+1):
		cual=int(np.floor(picos[j]))
		desp=mapa[i,cual]

		B1=fp[int(cual-ancho):int(cual+ancho+1),int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		B1 -= np.min(B1)
		B1 = B1/np.max(B1) 
		cualn=int(np.floor(picos[j+1]))
		desp=mapa[i,cualn]
		B2=fp[int(cualn-ancho):int(cualn+ancho+1),int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		B2 -= np.min(B2)
		B2 = B2/np.max(B2)
		delta=cualn-cual
		for cual2 in np.arange(cual,cual+delta+1):
			alfa=(cual2-cual)/(delta+1)
			#A[np.floor((cual2-ancho)*(ov+1)):np.floor((cual2+ancho+1)*(ov+1)),cual2]=(1.-alfa)*np.ravel(B1)+alfa*np.ravel(B2)
			#A[int((cual2-ancho)*(ov+1)):int((cual2+ancho+1)*(ov+1)),int(cual2)]=(1.-alfa)*np.ravel(B1)+alfa*np.ravel(B2)
			A2.rows[cual2] = np.arange(((cual2-ancho)*(ov+1)).astype(int),((cual2+ancho+1)*(ov+1)).astype(int))
			A2.data[cual2] = (1.-alfa)*np.ravel(B1)+alfa*np.ravel(B2)
			
	cual0=np.floor(picos[jini])
	cualf=np.floor(picos[jfin+1])
	#A2=A[(cual0-ancho)*(ov+1):cualf*(ov+1),cual0:cualf]
	


	nombre='extraction_methods/Extraction_files/REF/Amatrix_order'+str(i)+'_lane1_ancho'+str(ancho)+'_test.pkl'
	fout=open(nombre,'w')
		
	As= A2.tocsr().transpose()
		
	# COmpute the largest singular value....it will hint on the relaxation parameter for Landweber's
	U,l,v=sparse.linalg.svds(A,1)
	pickle.dump({'Amatrix':As,'Limits':(jini,jfin),'Width':ancho,'Lambda_max':l},fout)
	fout.close()
	print(nombre)
	print("Tiempo de calculo para este orden: {0}secs".format(time.time()-t0))
	