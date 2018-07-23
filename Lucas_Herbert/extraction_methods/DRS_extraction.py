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


workdir='./temp'
CCDsize=4640
ref=2320 
ancho=1

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

datafileFlat = datafile + "temp/Narval_20180313_180054_f10.fts"
datafileStar = datafile + "13mar18/Narval_20180314_003846_st0.fts"
datafileThar = datafile + "13mar18/Narval_20180313_181059_th0.fts"
f=pyfits.open(datafileThar)
img=f[0].data
f.close()
f = pyfits.open(datafileFlat)
img2=f[0].data
f.close()

# ============ Ordre =====================================================
orden= 2
#=========================================================================
to = time.time()
all_thar = []
for orden in range(0,35) :
	
	print("=================== ORDER : "+str(orden)+" ====================")
	nombre='extraction_methods/Extraction_files/REF/Amatrix_order{0}_lane1.pkl'.format(orden)
	
	o0=interorden[orden]
	o1=interorden[orden+1]
	ov=interpol[orden]
	
	fout=open("extraction_methods/Extraction_files/REF/FP_map.pkl","r")
	a=pickle.load(fout)
	FPpics=a['picos']
	fout.close()
	
	picos=FPpics[orden]
	
	Y=np.zeros(CCDsize*(ov+1))
	YFlat = np.zeros(CCDsize*(ov+1))
	arturos=np.shape(mapa)[1]
	Souv=np.zeros((arturos,ov+1))
	
	nombre='extraction_methods/Extraction_files/REF/Amatrix_order{0}_lane1.pkl'.format(orden)
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
		
	
	f = open("extraction_methods/Extraction_files/second_membre_CCD_ThAr_ord15",'r')
	bm = pickle.load(f) 
	bm = sparse.csr_matrix(bm)
	bm = bm.toarray()
	bm = [bm[i][0] for i in range(len(bm))]
	bm = sparse.csr_matrix(bm)
	f.close()
	f = open("extraction_methods/Extraction_files/second_membre_CCD_flat_ord15",'r')
	bmflat = pickle.load(f) 
	bmflat = sparse.csr_matrix(bmflat)
	bmflat = bmflat.toarray()
	bmflat = [bmflat[i][0] for i in range(len(bmflat))]
	bmflat = sparse.csr_matrix(bmflat)
	f.close()
	#f = "extraction_methods/Extraction_files/Bmat_ThAr_ord14_lane1.npz"
	f = "extraction_methods/Extraction_files/Bmat_ThAr_ord15_lane1.npz"
	Am = sparse.load_npz(f)
	
	
	
	
	def conjugate_grad_prec(A,Y):
		
		"""
		Cette fonction implemente l'algorithme du gradient conjugue preconditionne, permettant de resoudre un systeme lineaire du type Ax = Y , meme mal conditionne en un temps optimal. Cependant, cela exige A auto-adjointe c'est pourquoi on va commencer par transformer A en A.T*A pour obtenir C auto-adjointe. On utilisera la décomposition Cholesky incomptlete pour le preconditionnement.
		Input : 
		- A : matrice du systeme
		- Y : vecteur second membre
		Ouput :
		- new_xk : solution approchee du systeme
		"""
	
	
		# ========== Precondtioned conjugate gratient method ==========
		# https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
		sigma = np.sqrt( Y )
		sigma = np.diag(sigma)
		
		# Pas besoin de preconditionner donc on peut directement utiliser C1 et b1
		
		beg = time.time()
		C = np.dot(sigma,A)
		C = np.dot(A.T,C)
		b = np.dot(A.T,np.sqrt(Y)*Y)
		
		# Calcul de la matrice M de preconditionnement
		G = ichol(C)
		M = np.dot(G.T,G)
		DM=np.zeros_like(M) 
		np.fill_diagonal(DM,np.diag(M))
		#M = DM
		M = np.linalg.inv(M)
		print(" ==================== CONDITION ==============")
		print(np.linalg.cond(np.dot(M,C)) )
		
		# Calcul de la solution simple avec Jacobi
		Cjacob = np.zeros_like(C)
		np.fill_diagonal(Cjacob,1/np.diag(C))
		xnaif = np.dot(Cjacob,b)
		xnaif = xnaif/np.max(xnaif)
		
		# Initialisation :
		x0 = np.zeros(np.shape(C)[0])
		r0 = b - np.dot(C,x0)
		z0 = np.dot(M,r0)
		p0 = r0
		
		xk = x0
		pk = p0
		rk = r0
		zk = z0
		ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xk)) + np.linalg.norm(b) ) # e0 = 1.0
		end = time.time()
		print("============== COMPUTATION TIME ================")
		print(end-beg) 
		
		# Iteration
		beg = time.time()	
		while( ek > 0.95) :
	
			ak = np.dot(rk.T,zk) / np.dot(pk.T, np.dot(C,pk) )
			new_xk = xk + ak*pk
			new_rk = rk - ak*np.dot(C,pk)
			new_zk = np.dot(M,new_rk)
			bk = np.dot(new_zk.T,new_rk) / np.dot(zk.T,rk)
			new_pk = new_zk + bk*pk
			xk = new_xk
			rk = new_rk
			pk = new_pk
			zk = new_zk
			
			ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xk)) + np.linalg.norm(b) ) 
			#print("ek :",ek)
			error = np.linalg.norm(rk)/np.linalg.norm(b)
			#print("error :",error)
		end = time.time()
		print("========= ITERATION TIME ========")
		print(end-beg) 
		# Plotting the results
		new_xkf = new_xk/np.max(new_xk)	
		plt.figure(1)	
		plt.plot(new_xkf,color='purple')
		plt.plot(xnaif,color='black')
		plt.show()
		return(xk)
	
	def biconjugate_grad(A,Y):
		""" /!\ NON FONCTIONNEL, NECESSITE A CARREE /!\
		Voir ici :  https://en.wikipedia.org/wiki/Biconjugate_gradient_method
		Cette fonction implemente l'algorithme du gradient biconjugue, permettant de resoudre un systeme lineaire du type Ax = Y sans exiger que A soit auto-adjointe (dans notre cas, que A.T = A). 
		Input : 
		- A : matrice du systeme
		- Y : vecteur second membre
		Ouput :
		- new_xk : solution approchee du systeme
		"""
		# Initialisation
		b = np.zeros((len(Y2),1))
		for i in range(len(Y2)):
			b[i][0] = Y2[i]
		x0 = np.zeros((len(A[0]),1))
		x0_ = x0.T
		b_ = b.T
		r0 = b - A.dot(x0)
		r0_ = b_ - x0_.dot(A.T)
		p0 = r0
		p0_ = r0_
		
		xk = x0
		xk_ = x0_
		rk = r0
		rk_ = r0_
		pk = p0
		pk_ = p0_
		
		
		# Iteration
		beg = time.time()
		for k in range(20):
				alphak = (rk_.dot(rk)) / ( (pk_.dot(A)).dot(pk))
				new_xk = xk + alphak*pk
				new_xk_ = xk_ + alphak*pk_
				new_rk = rk - alphak*A.dot(pk)
				new_rk_ = rk_ - alphak*pk_.dot(A)
				betak = (new_rk_.dot(new_rk)) / (rk_.dot(rk))
				new_pk = new_rk + betak*pk
				new_pk_ = new_rk_ + betak*pk_
				
				xk = new_xk
				xk_ = new_xk_
				rk = new_rk
				rk_ = new_rk_
				pk = new_pk
				pk_ = new_pk_
			
		end = time.time()
		print("========== ITERATION TIME ==========")
		print(end-beg)
		
		xkf = xk / np.max(xk)
		plt.figure(2)
		plt.plot(xkf,color = 'green')
		plt.show()
		return(xk)
	
	
	def conjugate_grad_prec_sparse(A,Y):
		
		"""
		Cette fonction implemente l'algorithme du gradient conjugue preconditionne, permettant de resoudre un systeme lineaire du type Ax = Y , meme mal conditionne en un temps optimal. Cependant, cela exige A auto-adjointe c'est pourquoi on va commencer par transformer A en A.T*A pour obtenir C auto-adjointe. De plus dans cette version on effctuera les gros calculs en format sparse pour gagner du temps.
		Input : 
		- A : matrice du systeme en format sparse
		- Y : vecteur second membre en format sparse
		Ouput :
		- new_xk : solution approchee du systeme
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
		omega = 0.0005
		#================================================
		
		# Initialisation :
		
		x0 = sparse.csr_matrix(xnaif)    #sparse.csr_matrix(xnaif)    #sparse.csr_matrix( np.zeros(np.shape(C)[0]) ).T
		r0 = b - omega*C.dot(x0)
		z0 = M.dot(r0)
		p0 = r0
	
		xk = x0
		pk = p0
		rk = r0
		zk = z0
		#ek = np.linalg.norm(rk) / ( np.linalg.norm(C.dot(xk)) + np.linalg.norm(b) ) 
		# e0 = 1.0
		end = time.time()
		print("============== COMPUTATION TIME ================")
		print(end-beg) 
		
		# ========== Iteration ==========
		beg = time.time()	
		
		# Criteres d'arret
		
		iter = 500 # pour eviter d'iterer dans le vide
		
		Morozov=0.1*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
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
	
	a = time.time()
	plt.figure(1)
	
	thar_result = conjugate_grad_prec_sparse(A,Y3)
	flat_result = conjugate_grad_prec_sparse(A,Y3Flat)
	thar = thar_result[0]
	tharnaif = thar_result[1]
	flat = flat_result[0]
	flatnaif = flat_result[1]
	thar = thar / np.max(thar)
	flat = flat / np.max(flat)
	thar2 = thar/(flat) 
	thar2 = thar2 - np.min(thar2)
	thar2 = thar2 / np.max(thar2)
	plt.plot(thar2,color = 'green')
	tharnaif2 = tharnaif/flatnaif
	tharnaif2 = tharnaif2 / np.max(tharnaif2)
	plt.plot(tharnaif2,color = 'black')
	plt.figure(2)
	plt.plot(thar2,color = 'red')
	#plt.plot(tharnaif2,color = 'black')
	plt.show()
	
	
	all_thar.extend(thar2)
	
	print("_______________ TEMPS TOTAL ORDRE ________________")
	print(time.time() - a)
		
plt.figure(3)
plt.plot(all_thar,color = 'red')
plt.show()
#plt.close(1)
plt.close(2)

print("=_=_=_=_=_=_=_= TEMPS TOTAL GLOBAL =_=_=_=_=_=_=_=_=")
print(time.time()-to)