#! /usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
import extraction_methods.Chelli_Module as chlm
from scipy.optimize import curve_fit
from scipy import sparse
import time as time
import numpy.ma as ma

CCDsize=4640
ref=2320 


""" Fonctions """

# Factorisation incomplete de Choleski
def ichol(A):
	
	beg = time.time()
	G = np.copy(A)
	n = np.shape(G)[0]

	for k in range(0,n):
		G[k,k] = np.sqrt(G[k,k])
		for i in range(k+1,n):
			if (G[i,k] != 0):
				G[i,k] = G[i,k]/G[k,k]            
			

		for j in range(k+1,n):
			for i in range(j,n):
				if (G[i,j]!=0):
					G[i,j] = G[i,j]-G[i,k]*G[j,k]  

	for i in range(0,n):
		for j in range(i+1,n):
			G[i,j] = 0
	end = time.time()
	print(" ===== ichol TIME = "+str(end-beg)+" =====")
	return(G)



			
""" archives """

# def conjugate_grad_prec_sparsev2(A,Y,x0_file):
# 	
# 	"""
# 	Cette fonction implemente l'algorithme du gradient conjugue preconditionne, permettant de resoudre un systeme lineaire du type Ax = Y , meme mal conditionne en un temps optimal. Cependant, cela exige A auto-adjointe c'est pourquoi on va commencer par transformer A en A.T*A pour obtenir C auto-adjointe. De plus dans cette version on effctuera les gros calculs en format sparse pour gagner du temps.
# 	Input : 
# 	- A : matrice du systeme en format sparse
# 	- Y : vecteur second membre en format sparse
# 	- x0_file : fichier ou se situe le vecteur d'initialisation
# 	Ouput :
# 	- result : solutions approchees du systeme : xk et xnaif
# 	"""
# 
# 	# ========== Precondtioned conjugate gratient method ==========
# 	# https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
# 	beg = time.time()
# 	
# 	lenY = np.shape(Y)[1]
# 	# All the following matrix are in sparse format!
# 	sigma = 1 / np.sqrt( Y.toarray() )
# 	#sigma = [ sigma[i][0] for i in range(len(sigma)) ]
# 	sigma = sparse.spdiags(sigma,0,lenY,lenY)
# 	sigma = sparse.csr_matrix(sigma)
# 
# 	C = sigma.dot(A)
# 	C = A.T.dot(C)
# 	b = A.T.dot(sigma)
# 	b = b.dot(Y.T)
# 
# 	# Calcul de la solution simple avec Jacobi
# 	Cjacobi = np.zeros_like(C.toarray())
# 	np.fill_diagonal(Cjacobi,1/np.diag(C.toarray()))
# 	xnaif = np.dot(Cjacobi,b.toarray())
# 	xnaifn = xnaif/np.max(xnaif)
# 	
# 	print("==================== CONDITION ==============")
# 	# On preconditionne avec Jacobi pour aller plus vite en temps de calcul
# 	M = Cjacobi #np.eye(np.shape(C)[0])
# 	#M = M.T.dot(M)
# 	#M = np.linalg.inv(M)
# 	M = sparse.csr_matrix(M)
# 	#print(np.linalg.cond( (M.dot(C)).toarray() ))
# 	
# 	#================================================
# 	# Facteur de relaxation? 
# 	omega = 10
# 	#================================================
# 	
# 	# Initialisation :
# 	f = open(x0_file,'r')
# 	x0 = pickle.load(f)
# 	f.close()
# 	x0 = x0 / np.max(x0)
# 	plt.figure(12)
# 	plt.plot(xnaifn,'black')
# 	plt.plot(x0,'orange')
# 	plt.show()
# 	
# 	x0 = sparse.csr_matrix(x0).T
# 	print(np.shape(C))
# 	print(np.shape(x0))
# 	
# 	r0 = b - omega*C.dot(x0)
# 	z0 = M.dot(r0)
# 	p0 = r0
# 
# 	xk = x0
# 	pk = p0
# 	rk = r0
# 	zk = z0
# 	#ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xk)) + np.linalg.norm(b) ) 
# 	# e0 = 1.0
# 	end = time.time()
# 	print("============== COMPUTATION TIME ================")
# 	print(end-beg) 
# 	
# 	# ========== Iteration ==========
# 	beg = time.time()	
# 	
# 	# Criteres d'arret
# 	
# 	iter = 1000 # pour eviter d'iterer dans le vide
# 	
# 	Morozov=0.75*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
# 	print("Morozov : "+ str(Morozov))
# 	
# 	diff = 1 # initialisation
# 	new_diff = 2
# 	k = 0
# 	
# 	#and (diff>Morozov)
# 	while ( ( np.abs(new_diff - diff) > 1e-10) and (k < iter) and (diff>Morozov)):
# 		
# 		diff = np.linalg.norm(rk.toarray()) / (np.linalg.norm(xk.toarray())+1e-20) #The "+1e-20" stands for the first case, where xk = 0 and the division fataly means nothing
# 
# 		ak =  (rk.T.dot(zk)).toarray()[0][0] /  (pk.T.dot(C).dot(pk)).toarray()[0][0]  # les produits scalaires sont mal geres par sparse d'ou l'ecriture etrange
# 		new_xk = xk + omega*ak*pk
# 		new_rk = rk - omega*ak*C.dot(pk)
# 		new_zk = M.dot(new_rk)
# 		bk = (new_zk.T.dot(new_rk)).toarray()[0][0] / (zk.T.dot(rk) ).toarray()[0][0]
# 		new_pk = new_zk + omega*bk*pk
# 		xk = new_xk
# 		rk = new_rk
# 		pk = new_pk
# 		zk = new_zk
# 		k += 1
# 		
# 		new_diff = np.linalg.norm(rk.toarray()) / np.linalg.norm(xk.toarray())
# 		
# 	end = time.time()
# 	
# 	xk = list(xk.T.toarray()[0]) # conversion du format sparse vers un array classique exploitable
# 	if (k == iter):
# 		print("Maximum iteration has been reached")
# 	print("Iterations number : "+str(k))
# 	print("========= ITERATION TIME ========")
# 	print(end-beg) 
# 	# Plotting the results
# 	xkn = xk/np.max(xk)	# on normalise xk
# 	plt.figure(1)	
# 	plt.plot(xkn,color='green')
# 	plt.plot(xnaifn,color='black')
# 	plt.show()
# 	return(xk,xnaif)

def conjugate_grad_prec_sparse(A,Y):
	
	"""
	Cette fonction implemente l'algorithme du gradient conjugue preconditionne, permettant de resoudre un systeme lineaire du type Ax = Y , meme mal conditionne en un temps optimal. Cependant, cela exige A auto-adjointe c'est pourquoi on va commencer par transformer A en A.T*A pour obtenir C auto-adjointe. De plus dans cette version on effctuera les gros calculs en format sparse pour gagner du temps.
	Input : 
	- A : matrice du systeme en format sparse
	- Y : vecteur second membre en format sparse
	Ouput :
	- result : solutions approchees du systeme : [xk,xnaif]
	"""

	# ========== Precondtioned conjugate gratient method ==========
	# https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
	beg = time.time()
	
	lenY = np.shape(Y)[1]
	# All the following matrix are in sparse format!
	sigma = 1 / np.sqrt( Y.toarray() )
	#sigma = [ sigma[i][0] for i in range(len(sigma)) ]
	sigma = sparse.spdiags(sigma,0,lenY,lenY)
	sigma = sparse.csr_matrix(sigma)
	print( np.shape(sigma))
	print(np.shape(A))
	C = sigma.dot(A)
	C = A.T.dot(C)
	b = A.T.dot(sigma)
	b = b.dot(Y.T)

	# Calcul de la solution simple avec Jacobi
	Cjacobi = np.zeros_like(C.toarray())
	np.fill_diagonal(Cjacobi,1/np.diag(C.toarray()))
	xnaif = np.dot(Cjacobi,b.toarray())
	xnaifn = xnaif/np.max(xnaif)
	
	print("==================== CONDITION ==============")
	# On preconditionne avec Jacobi pour aller plus vite en temps de calcul
	M = Cjacobi #np.eye(np.shape(C)[0])
	#M = M.T.dot(M)
	#M = np.linalg.inv(M)
	M = sparse.csr_matrix(M)
	#print(np.linalg.cond( (M.dot(C)).toarray() ))
	
	#================================================
	# Facteur de relaxation? 
	omega = 0.01
	#================================================
	
	# Initialisation :
	
	x0 = sparse.csr_matrix(xnaif)    #sparse.csr_matrix(xnaif)    #sparse.csr_matrix( np.zeros(np.shape(C)[0]) ).T
	# print(np.shape(C))
	# print(np.shape(x0))
	r0 = b - omega*C.dot(x0)
	z0 = M.dot(r0)
	p0 = r0

	xk = x0
	pk = p0
	rk = r0
	zk = z0
	#ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xkfout = open("extraction_methods/Extraction_files/temp/x0thar_ancho0_order_"+str(orden),'w')
	# e0 = 1.0
	end = time.time()
	print("============== COMPUTATION TIME ================")
	print(end-beg) 
	
	# ========== Iteration ==========
	beg = time.time()	
	
	# Criteres d'arret
	
	iter = 1000 # pour eviter d'iterer dans le vide
	
	Morozov=0.75*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
	print("Morozov : "+ str(Morozov))
	
	diff = 1 # initialisation
	new_diff = 2
	k = 0
	
	#and (diff>Morozov)
	while ( ( np.abs(new_diff - diff) > 1e-10) and (k < iter) and (diff>Morozov)):
		
		diff = np.linalg.norm(rk.toarray()) / (np.linalg.norm(xk.toarray())+1e-20) #The "+1e-20" stands for the first case, where xk = 0 and the division fataly means nothing

		ak =  (rk.T.dot(zk)).toarray()[0][0] /  (pk.T.dot(C).dot(pk)).toarray()[0][0]  # les produits scalaires sont mal geres par sparse d'ou l'ecriture etrange
		new_xk = xk + omega*ak*pk
		new_rk = rk - omega*ak*C.dot(pk)
		new_zk = M.dot(new_rk)
		bk = (new_zk.T.dot(new_rk)).toarray()[0][0] / (zk.T.dot(rk) ).toarray()[0][0]
		new_pk = new_zk + omega*bk*pk
		xk = new_xk
		rk = new_rk
		pk = new_pk
		zk = new_zk
		k += 1
		
		new_diff = np.linalg.norm(rk.toarray()) / np.linalg.norm(xk.toarray())
		
	end = time.time()
	
	xk = list(xk.T.toarray()[0]) # conversion du format sparse vers un array classique exploitable
	if (k == iter):
		print("Maximum iteration has been reached")
	print("Iterations number : "+str(k))
	print("========= ITERATION TIME ========")
	print(end-beg) 
	# Plotting the results
	xkn = xk/np.max(xk)	# on normalise xk
	plt.figure(1)	
	plt.plot(xkn,color='green')
	plt.plot(xnaifn,color='black')
	plt.show()
	return(xk,xnaif)



""" Traitement des donnees """



f=open('extraction_methods/Extraction_files/REF/oder_reference.pkl','r')
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

datafile = "extraction_methods/Extraction_files/"

datafileFlat = datafile + "temp/simu00_20180725_141934_fp0.fts"
datafileStar = datafile + "13mar18/Narval_20180314_003846_st0.fts"
datafileThar = datafile + "13mar18/Narval_20180313_181059_th0.fts"
f=pyfits.open(datafileThar)
img=f[0].data
f.close()
f = pyfits.open(datafileFlat)
img2=f[0].data
f.close()

# Normalisation du second membre
#img = img/img2
# ######################


# On lit l'atlas solaire
atlas="extraction_methods/Extraction_files/SolarSpectrum.pkl"
f=open(atlas,'r')
sol=pickle.load(f)
f.close()

# Fichier utilise pour la calibration en longueurs d'onde
calibfile='calibration_methods/Calibration_files/th_calibered.fits'
f=pyfits.open(calibfile)
a=f[1].data
Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))
wvl=a['wavelength_lane1']


ancho_list = ['8'] #['0_flat_edge','5_flat_edge','8_flat_edge','11_flat_edge','11_divided_by_flat'] #'3','5','8','15_flat_edge'
graphs3, graphs4 , ancho_spectra, ancho_flat, ancho_thar_brut = [], [], [], [], []

for ancho_str in ancho_list :
	
	print("################ ancho : " + str(ancho_str) + " ################")
	
	to = time.time()
	all_thar = []
	all_w = []
	all_Solw = []
	all_SolI = []
	all_flat = []
	all_thar_brut = []
	#	===== Selection de l'ordre =======
	
	for orden in range(10,20) :
		
		print("=================== ORDER : "+str(orden)+" ====================")
		
		o0=interorden[orden]
		o1=interorden[orden+1]
		ov=interpol[orden]
		
		fout=open("extraction_methods/Extraction_files/REF/FP_map_new.pkl","r")
		a=pickle.load(fout)
		FPpics=a['picos']
		fout.close()
		
		picos=FPpics[orden]
		
		Y=np.zeros(CCDsize*(ov+1))
		YFlat = np.zeros(CCDsize*(ov+1))
		arturos=np.shape(mapa)[1]
		Souv=np.zeros((arturos,ov+1))
		nombre='extraction_methods/Extraction_files/REF/Amatrix_order'+str(orden)+'_lane1_ancho'+str(ancho_str)+'_simu.pkl'
		f=open(nombre,'r')
		matriz=pickle.load(f)
		Acompleta=matriz['Amatrix']#.toarray()
		f.close()
		jini,jfin=matriz['Limits']
		ancho=matriz['Width']
		lambda_max=0.33 #matriz['Lambda_max']
		f.close()
	
		
		cual0=int(np.floor(picos[jini+2]))
		cualf=int(np.floor(picos[jfin-2]))
		
		for j in np.arange(0,arturos):
			desp=mapa[orden,j]
			B=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
			BFlat = img2[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
			Y[(ov+1)*j:(ov+1)*(j+1)]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
			YFlat[(ov+1)*j:(ov+1)*(j+1)]=img2[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
			Souv[j,:]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
		
		
		A=Acompleta[(cual0-ancho)*(ov+1):(cualf+ancho+1)*(ov+1),cual0:cualf]#[(cual0-ancho)*(ov+1):cualf*(ov+1),cual0:cualf]
		#Y2=Y[int((cual0-ancho)*(ov+1)):int(cualf*(ov+1))]
		Y2=Y[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
		Y2Flat=YFlat[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
				
		Y3 = sparse.csr_matrix(Y2)
		Y3Flat = sparse.csr_matrix(Y2Flat)
		
		# Calibration
		angstroms=wvl[Orderlimit[orden]:Orderlimit[orden+1]]
		w=angstroms[int(cual0):int(cualf)]
		# Je vais chercher dans l'atlas le morceaux de spectre choisi
		donde=ma.masked_inside(sol["Wavelength"],w[0],w[-1])
		Solw=np.asarray(sol["Wavelength"])[donde.mask]
		SolI=np.asarray(sol["Intensity"])[donde.mask]

		
	
		""" Calculs """
# DEBUT COMMENTAIRE
		a = time.time()
		plt.figure(1)
		print("ooooo Computing spectrum ooooo")
		thar_result = conjugate_grad_prec_sparse(A,Y3) #"extraction_methods/Extraction_files/temp/x0thar_ancho0_order_"+str(orden))
		print("ooooo Computing flat ooooo")
		#flat_result = conjugate_grad_prec_sparse(A,Y3Flat) #"extraction_methods/Extraction_files/temp/x0flat_ancho0_order_"+str(orden))
		thar = thar_result[0]
		#tharnaif = thar_result[1]
		#flat = flat_result[0]
		#flatnaif = flat_result[1]
		thar = thar / np.max(thar)
		#flat = flat / np.max(flat)
		#thar2 = thar/flat 
		#thar2 = thar2 - np.min(thar2)
		#thar2 = thar2 / np.max(thar2)
		#plt.plot(thar2,color = 'green')
		#tharnaif2 = tharnaif/flatnaif
		#tharnaif2 = tharnaif2 / np.max(tharnaif2)
		#plt.plot(tharnaif2,color = 'black')
		#plt.figure(2)
		#plt.plot(thar2,color = 'red')
		#plt.plot(Solw,SolI,'black')
		#plt.show()
		
		fout = open("extraction_methods/Extraction_files/reduced/reduced_thar_"+str(orden)+"_ancho"+str(ancho_str)+".pkl",'w')
		pickle.dump(thar_result[0],fout)
		fout.close()
		all_thar_brut.extend(thar)
		#all_thar.extend(thar2)
		all_w.extend(w)
		#all_flat.extend(flat)
		# all_Solw.extend(Solw)
		# all_SolI.extend(SolI)
		
		print("_______________ TEMPS TOTAL ORDRE ________________")
		print(time.time() - a)
			
	# plt.figure(3)
	# graph3, = plt.plot(all_w,all_thar,label='Ancho ='+str(ancho_str))
	# plt.legend(handles=[graph3])
	# plt.plot(all_Solw,all_SolI,color='black')
	# plt.show()
	plt.figure(4)
	graph4, = plt.plot(all_thar_brut,label='Ancho ='+str(ancho_str))
	plt.legend(handles=[graph4])
	plt.show()
	
	print("=_=_=_=_=_=_=_= TEMPS TOTAL GLOBAL =_=_=_=_=_=_=_=_=")
	print(time.time()-to)
	
	
	graphs3.append(graph3)
	graphs4.append(graph4)
	ancho_spectra.append(all_thar)
	#ancho_flat.append(all_flat)
	ancho_thar_brut.append(all_thar_brut)

plt.figure(10)
for i in range(len(ancho_list)) :
	plt.plot(ancho_spectra[i])
	#plt.plot(ancho_flat[i],color='black')
	#plt.plot(ancho_thar_brut[i])
plt.legend([g for g in graphs4],[ancho for ancho in ancho_list])
plt.show()

plt.figure(11)
for i in range(len(ancho_list)) :
	plt.plot(all_w, ancho_spectra[i])
plt.legend([g for g in graphs3],[ancho for ancho in ancho_list])
plt.show()
# FIN COMMENTAIRE












##########################################################################################
def conjugate_grad_prec_sparse_version_Martin(A,Y):
	
	"""
	Cette fonction implemente l'algorithme du gradient conjugue preconditionne, permettant de resoudre un systeme lineaire du type Ax = Y , meme mal conditionne en un temps optimal. Cependant, cela exige A auto-adjointe c'est pourquoi on va commencer par transformer A en A.T*A pour obtenir C auto-adjointe. De plus dans cette version on effctuera les gros calculs en format sparse pour gagner du temps.
	Input : 
	- A : matrice du systeme en format sparse
	- Y : vecteur second membre en format sparse
	Ouput :
	- result : solutions approchees du systeme : [xk,xnaif]
	"""

	# ========== Precondtioned conjugate gratient method ==========
	# https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
	beg = time.time()
	
	lenY = np.shape(Y)[1]
	# All the following matrix are in sparse format!
	sigma = 1 / np.sqrt( Y.toarray() )
	#sigma = [ sigma[i][0] for i in range(len(sigma)) ]
	sigma = sparse.spdiags(sigma,0,lenY,lenY)
	sigma = sparse.csr_matrix(sigma)

	C = sigma.dot(A)
	C = A.T.dot(C)
	b = A.T.dot(sigma)
	b = b.dot(Y.T)

	# Calcul de la solution simple avec Jacobi
	Cjacobi = np.zeros_like(C.toarray())
	np.fill_diagonal(Cjacobi,1/np.diag(C.toarray()))
	xnaif = np.dot(Cjacobi,b.toarray())
	xnaifn = xnaif/np.max(xnaif)
	
	print("==================== CONDITION ==============")
	# On preconditionne avec Jacobi pour aller plus vite en temps de calcul
	M = Cjacobi #np.eye(np.shape(C)[0])
	#M = M.T.dot(M)
	#M = np.linalg.inv(M)
	M = sparse.csr_matrix(M)
	#print(np.linalg.cond( (M.dot(C)).toarray() ))
	
	#================================================
	# Facteur de relaxation? 
	omega = 0.001
	#================================================
	
	# Initialisation :
	
	x0 = sparse.csr_matrix(xnaif)    #sparse.csr_matrix(xnaif)    #sparse.csr_matrix( np.zeros(np.shape(C)[0]) ).T
	# print(np.shape(C))
	# print(np.shape(x0))
	r0 = b - omega*C.dot(x0)
	z0 = M.dot(r0)
	p0 = r0

	xk = x0
	pk = p0
	rk = r0
	zk = z0
	#ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xkfout = open("extraction_methods/Extraction_files/temp/x0thar_ancho0_order_"+str(orden),'w')
	# e0 = 1.0
	end = time.time()
	print("============== COMPUTATION TIME ================")
	print(end-beg) 
	
	# ========== Iteration ==========
	beg = time.time()	
	
	# Criteres d'arret
	
	iter = 1000 # pour eviter d'iterer dans le vide
	
	Morozov=1e-2*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
	print("Morozov : "+ str(Morozov))
	
	diff = 1 # initialisation
	new_diff = 2
	k = 0
	
	#and (diff>Morozov)
	while ( ( np.abs(new_diff - diff) > 1e-12) and (k < iter) and (diff>Morozov)):
		
		diff = np.linalg.norm(rk.toarray()) / (np.linalg.norm(xk.toarray())+1e-20) #The "+1e-20" stands for the first case, where xk = 0 and the division fataly means nothing

		ak =  (rk.T.dot(zk)).toarray()[0][0] /  (pk.T.dot(C).dot(pk)).toarray()[0][0]  # les produits scalaires sont mal geres par sparse d'ou l'ecriture etrange
		new_xk = xk + omega*ak*pk
		new_rk = rk - omega*ak*C.dot(pk)
		new_zk = M.dot(new_rk)
		bk = (new_zk.T.dot(new_rk)).toarray()[0][0] / (zk.T.dot(rk) ).toarray()[0][0]
		new_pk = new_zk + omega*bk*pk
		xk = new_xk
		rk = new_rk
		pk = new_pk
		zk = new_zk
		k += 1
		
		new_diff = np.linalg.norm(rk.toarray()) / np.linalg.norm(xk.toarray())
		
	end = time.time()
	
	xk = list(xk.T.toarray()[0]) # conversion du format sparse vers un array classique exploitable
	if (k == iter):
		print("Maximum iteration has been reached")
	print("Iterations number : "+str(k))
	print("========= ITERATION TIME ========")
	print(end-beg) 
	# Plotting the results
	xkn = xk/np.max(xk)	# on normalise xk
	plt.figure(1)	
	plt.plot(xkn,color='green')
	plt.plot(xnaifn,color='black')
	plt.show()
	return(xk,xnaif)


# f = open("extraction_methods/Extraction_files/Matrices_Martin/Bmat_test_dirak_ord14_lane1.npz",'r')
# Am = sparse.load_npz(f)
# f.close()
# 
# f = open("extraction_methods/Extraction_files/Matrices_Martin/second_membre_ThAr_1306-3306")
# bm = pickle.load(f)
# f.close()
# bm = np.array(bm)
# bm = [bm[i][0] for i in range(len(bm))]
# bm = sparse.csr_matrix(bm)
# 
# 
# th_result = conjugate_grad_prec_sparse_version_Martin(Am,bm)
# th_naif = th_result[1] 
# th = th_result[0]
# th_naif = th_naif/np.max(th_naif)
# th = th / np.max(th)
# plt.figure(33)
# plt.plot(th,'red')
# plt.show()




	
	
	
	