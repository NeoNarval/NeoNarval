#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
from Chelli_Module import *
from scipy.optimize import curve_fit
from scipy import sparse


def Vfunction(x,a,p,d):
		gauss=a * np.exp(-(x-p)**2/d)
		return(-np.gradient(gauss)) 


workdir='./temp'
CCDsize=4640
ref=2320 
ancho=1

dir='13mar18/'

# Get FLAT done
files=sorted([f for f in os.listdir(dir) if f.endswith('fp0.fts')])
#f='/Users/arturo/NeoNARVAL/Espadons/1772472f.fits'

img=0.
for f in files:
	a=pyfits.open(dir+f)
	img=img+a[0].data
	a.close()

fp=img/len(files)

f=open('./REF/oder_reference.pkl','r')
a=pickle.load(f) 
f.close()
perfref0=a['Profil']
interorden=a['Orders']
interpol=a['Beams']

f=open("./REF/Order_Curvature.pkl","r")
a=pickle.load(f) 
f.close()
mapa=a['map']
Norders=len(interorden)

if (False): 	
	
	FPmap=[]
	FPpics=[]
	for i,order in enumerate(interorden[0:Norders-1]):
		perfil=np.zeros(4640)
		#img=np.zeros((50,4640))
		o0=int(order)
		o1=int(interorden[i+1])
		ov=int(interpol[i])
		
		for j,desp in enumerate(mapa[i,:]):

			perfil[j]=np.sum(fp[j,o0+int(np.floor(desp)):o0+ov+int(np.floor(desp))+1])
			#img[0:ov+1,j]=fp[j,o0+np.floor(desp):o0+ov+np.floor(desp)+1]
		
		der=np.gradient(perfil)
		#plt.plot(der)
		picos=[]
		for j in range(50,CCDsize-50):
			if (der[j]*der[j-1]<0) and (perfil[j]>np.mean(perfil[j-50:j+50])):
		
				try:
					popt, pcov = curve_fit(Vfunction, np.arange(20),der[j-10:j+10],p0=[perfil[j],10,4] )
					if 1<popt[2]<9:
						picos.append(j-10+popt[1])
						# plt.axvline(j-10+popt[1],linestyle='--')					
	# 					plt.plot(j-10+np.arange(20), Vfunction(np.arange(20), *popt),color='r')
				except:
					pass
		print("Order {0}, {1} peaks found".format(i,len(picos)))
		FPpics.append(picos)
		FPmap.append(perfil)
		
		fout=open("./REF/FP_map.pkl","w")
		pickle.dump({'picos':FPpics,'perfiles':FPmap,'info':"Positions of the FP peaks, order by order"},fout)
		fout.close()
else:

	fout=open("./REF/FP_map.pkl","r")
	a=pickle.load(fout)
	FPpics=a['picos']
	FPmap=a['perfiles']
	fout.close()


	# plt.plot(picos)
# 	plt.title('Order {0}'.format(i))
# 	plt.show()

for i in np.arange(0,12):#Norders):
	
	o0=int(interorden[i])
	o1=int(interorden[i+1])
	ov=int(interpol[i])
	
	picos=FPpics[i]
	arturos=len(FPmap[i])
	jini=75
	jfin=81
	
	
	

		
	A=np.zeros( ((ov+1)*CCDsize,arturos))

	for j in np.arange(jini,jfin+1):
		cual=int(picos[j])
		desp=mapa[i,cual]

		B1=fp[int(cual-ancho):int(cual+ancho+1),int(o0+np.floor(desp)) : int( o0+ov+np.floor(desp)+1 )]
		B1=B1/np.sum(B1)
		cualn=int(picos[j+1])
		desp=mapa[i,cualn]
		B2=fp[int(cualn-ancho):int(cualn+ancho+1),int(o0+np.floor(desp)):int(o0+ov+np.floor(desp))+1]
		B2=B2/np.sum(B2)
		delta=picos[j+1]-cual
		for cual2 in np.arange(cual,cual+delta+1):
			alfa=(cual2-cual)/(delta+1)
			A[int((cual2-ancho))*(ov+1):int((cual2+ancho+1))*(ov+1),int(cual2)]=(1.-alfa)*np.ravel(B1)+alfa*np.ravel(B2)


	cual0=np.floor(picos[jini])
	cualf=np.floor(picos[jfin+1])
	A2=A[int((cual0-ancho))*(ov+1):int(cualf)*(ov+1),int(cual0):int(cualf)]
	# plt.imshow(A[20000:23000,:],aspect='auto',interpolation='nearest')
# 	plt.colorbar()
# 	plt.title('Order {0}'.format(i))
# 	#plt.clim(200,400)
# 	plt.show()
#  	


	#print(np.linalg.cond(A2))
	C=np.dot(A2.T,A2)
	C=C/np.amax(C)
# 	
# 	for j in np.arange(0,(ov+1)):
# 		C[j,j]=1.
	# for j in np.arange(np.floor(picos[-2])+1,arturos):
# 		Cs[j,j]=1.
# 	stop
	G=np.linalg.cholesky(C)
	# plt.imshow(G,interpolation='nearest')
# 	plt.show()

	nombre='./REF/Amatrix_order{0}_lane1.pkl'.format(i)
	fout=open(nombre,'w')
		
	sG = sparse.csr_matrix(G) 
	As= sparse.csr_matrix(A2) 
		
	pickle.dump({'Amatrix':As,'Gmatrix':sG},fout)
	fout.close()
	print(nombre)
 	
	