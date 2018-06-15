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
import numpy.ma as ma


workdir='./temp'
CCDsize=4640
ref=2320 
ancho=5

# On lit l'atlas solaire sur a
# atlas="/Users/arturo/DeepStokes/SolarSpectrum.pkl"
# f=open(atlas,'r')
# sol=pickle.load(f)
# f.close()


f=open('extraction_methods/Extraction_files/REF/order_reference.pkl','r')
a=pickle.load(f) 
f.close()
perfref0=a['Profil']
interorden=a['Orders']
interpol=a['Beams']

f=open("extraction_methods/Extraction_files/REF/Order_Curvature.pkl","r")
a=pickle.load(f) 
f.close()
mapa=a['map']
Norders=len(interorden)


calibfile='calibration_methods/Calibration_files/th_calibered.fits'
f=pyfits.open(calibfile)
a=f[1].data
Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))
wvl=a['wavelength_lane1']


flatfile="extraction_methods/Extraction_files/temp/Narval_20180313_180054_f10.fts"
f=pyfits.open(flatfile)
flat=f[0].data
f.close()

#datafile="extraction_methods/Extraction_files/Moon/Narval_20170614_020541_st1.fts"
datafile="extraction_methods/Extraction_files/13mar18/Narval_20180313_181059_th0.fts"
f=pyfits.open(datafile)
img=f[0].data
f.close()

fout=open("extraction_methods/Extraction_files/REF/FP_map.pkl","r")
a=pickle.load(fout)
FPpics=a['picos']
fout.close()


sp=[]
calw=[]
SL=[]
WL=[]
plt.ion()
all_Xshow = []

for orden in range(0,35):
	t0=time.time()
	print("Empezamos orden {0}".format(60-orden))
	nombre='extraction_methods/Extraction_files/REF/Amatrix_order{0}_lane1.pkl'.format(orden)

	o0=interorden[orden]
	o1=interorden[orden+1]
	ov=interpol[orden]


	
	picos=FPpics[orden]
	# jini=5
	# jfin=len(picos)-10

	nombre='extraction_methods/Extraction_files/REF/Amatrix_order{0}_lane1.pkl'.format(orden)
	f=open(nombre,'r')
	matriz=pickle.load(f)
	Acompleta=matriz['Amatrix']#.toarray()
	jini,jfin=matriz['Limits']
	ancho=matriz['Width']
	lambda_max=0.33#matriz['Lambda_max']
	f.close()
# 	print("Todo leido en {0}s".format(time.time()-t0)) #0.3sec

	#for jini in np.arange(55,60,30):#np.arange(5,len(picos)-5,50):#np.arange(55,60,30):#np.arange(5,len(picos)-35,30):
		
# 		jfin=jini+50
# 		if jfin>len(picos)-6:
# 			jfin=len(picos)-6
# 		print("		FP {0} a {1}".format(jini,jfin))
	
	Y=np.zeros(CCDsize*(ov+1))
	FY=np.zeros(CCDsize*(ov+1))
	arturos=np.shape(mapa)[1]
	Souv=np.zeros((arturos,ov+1))

	for j in np.arange(0,arturos):
		desp=mapa[orden,j]
		B=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		Y[int((ov+1)*j):int((ov+1)*(j+1))]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		FY[int((ov+1)*j):int((ov+1)*(j+1))]=flat[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		Souv[j,:]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]

	cual0=np.floor(picos[jini+2])
	cualf=np.floor(picos[jfin-2])
	
	Y2=Y[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
	Sigma=1./np.sqrt(Y2)
	FY2=FY[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
	FSigma=1./np.sqrt(FY2)
	A=Acompleta[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1)),int(cual0):int(cualf)]#[(cual0-ancho)*(ov+1):cualf*(ov+1),cual0:cualf]

# 		print("Matrices llenas en {0}".format(time.time()-t0)) #+0.4sec

	angstroms=wvl[Orderlimit[orden]:Orderlimit[orden+1]]
	w=angstroms[int(cual0):int(cualf)]
	# Je vais chercher dans l'atlas le morceaux de spectre choisi
	# donde=ma.masked_inside(sol["Wavelength"],w[0],w[-1])
	# Solw=np.asarray(sol["Wavelength"])[donde.mask]
	# SolI=np.asarray(sol["Intensity"])[donde.mask]


# 		print("LibreEsprit encontrado en {0}s".format(time.time()-t0)) #+0.005 sec


	# AT=sparse.csc_matrix.tocsr(A.T)
	Th2=sparse.csr_matrix.dot(A.T, Sigma*Y2)  #np.inner(A.T,np.sqrt(Y2)*Y2))

# 		print("Th2 calculado en {0}s".format(time.time()-t0)) #0.001s

	Sig2=sparse.lil_matrix((Y2.shape[0],Y2.shape[0]),dtype=np.float)
	Sig2.setdiag(Sigma)

	#Get the diagonal of C ....my way
	D1=np.array(sparse.csr_matrix.sum(Sig2*A.multiply(A),0)).ravel()
	# and solve directly Jacobi like
	X=Th2/D1

	print("		Jacobi resuelto en {0}s".format(time.time()-t0)) #0.02secs


	Fh2=sparse.csr_matrix.dot(A.T, FSigma*FY2) #np.inner(A.T,np.dot(FSigma,FY2))
# 	

	Sig2.setdiag(FSigma)
# 	
	FD1=np.asarray(sparse.csr_matrix.sum(Sig2*A.multiply(A),0)).ravel()
	FX=Fh2/FD1

	FX=FX/np.mean(FX)

	X0=X/FX
	X0=X0/np.amax(X0)

	print("Lo mismo para el flat en {0}s".format(time.time()-t0)) #+10sec, 0.1sec sin pinv




	wrelax=0.5   # <2./lambda_max**2
	X=X0.copy()
	diff=1.
	diff0=0.
	
	
	Morozov=.75*np.median(Sigma)  #Value of tau=1. for the Moon
	
	while (diff>Morozov) and (np.abs(diff-diff0)>1e-10) :
		diff0=diff
		r=np.sqrt(Sigma)*sparse.csc_matrix.dot(A,X)-np.sqrt(Sigma)*Y2
		rF=np.sqrt(FSigma)*sparse.csc_matrix.dot(A,FX)-np.sqrt(FSigma)*FY2
	
		
		X=X-wrelax*sparse.csr_matrix.dot(A.T,r)
		FX=FX-wrelax*sparse.csr_matrix.dot(A.T,rF)
		diff=np.linalg.norm(r)/np.linalg.norm(X)
		
		
	print(diff,Morozov)

	Xshow=X/(FX+1e-40)
	Xshow=Xshow/np.amax(Xshow[10:-10])
	print("		Landweber converge en {0}s".format(time.time()-t0)) #+1sec
	plt.figure(1)
	plt.plot(Xshow,color = 'blue')
	plt.show()
	# if len(SolI)>0:
	# 	slope=(1.-np.amin(Xshow[10:-10]))/(np.amax(SolI)-np.amin(SolI))
	# 	off=1.-slope*np.amax(SolI)
	# 	SL.extend(SolI*slope+off)
	# 	WL.extend(Solw)
	# 
	all_Xshow.extend(Xshow)
	sp.extend(FX)#Xshow)	
	calw.extend(w)
	#plt.cla()
	# asi=np.argsort(calw)
# 	plt.plot(np.asarray(calw)[asi],np.asarray(sp)[asi])
	#plt.plot(WL,SL,'g')
	#plt.plot(calw,sp,'r')
	
	#plt.ylim(0,1)
	# for i in np.arange(jini,jfin+1):
# 		plt.axvline(w[picos[i]-cual0+7],linestyle='--')
	#plt.draw()
plt.figure(3)
plt.plot(all_Xshow,color = 'blue')
plt.show()


