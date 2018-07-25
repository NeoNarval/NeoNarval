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
#import extraction_methods.DRS_extraction_Lucas as DRS_Lucas

for k in range(100) :
    plt.close(k)

CCDsize=4640
ref=2320 
global order_len 
order_len = 4612

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

datafileSimu = datafile + "simu00_20180724_074450_fp0.fts"
datafileFlat = datafile + "temp/Narval_20180313_180054_f10.fts"
datafileStar = datafile + "13mar18/Narval_20180314_003846_st0.fts"
datafileThar = datafile + "13mar18/Narval_20180313_181059_th0.fts"
datafileFP = datafile + "13mar18/Narval_20180313_181342_fp0.fts"
f=pyfits.open(datafileThar)
img=f[0].data

f.close()
f = pyfits.open(datafileFlat)
imgFlat=f[0].data
f.close()
f=pyfits.open(datafileSimu)
imgFP=f[0].data
f.close()


fout=open("extraction_methods/Extraction_files/REF/FP_map_new.pkl","r")
a=pickle.load(fout)
FPpics=a['picos']
fout.close()


# On lit l'atlas solaire
atlas="extraction_methods/Extraction_files/SolarSpectrum.pkl"
f=open(atlas,'r')
sol=pickle.load(f)
f.close()

# Fichier utilisé pour la calibration en longueurs d'onde
calibfile='calibration_methods/Calibration_files/th_calibered.fits'
f=pyfits.open(calibfile)
a=f[1].data
Orderlimit=np.concatenate(([0],np.where(np.diff(a['wavelength_lane1'])<0)[0],[len(a['wavelength_lane1'])-1]))
wvl=a['wavelength_lane1']
spectrum = a['intensity_lane1']

# Code  pour afficher le vecteur CCD en matrice 2D

ep_voie = interpol # liste des epaisseurs de chaque voie 
ord_centres= interorden # liste des ordonnees du centre de chaque ordre (axe Y allant de 0 a 2098)
delta_courb = mapa # matrice des listes (par ordre) des delta_pixel de courbure de l'ordre par rapport à la colonne centrale
M = np.zeros((4612,2098)) # matrice representant le CCD complet

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
    C = sigma.dot(A)
    C = A.T.dot(C)
    b = A.T.dot(sigma)
    b = b.dot(Y.T)
    
    # Calcul de la solution simple avec Jacobi
    Cjacobi = np.zeros_like(C.toarray())
    np.fill_diagonal(Cjacobi,1/np.diag(C.toarray()))
    xnaif = np.dot(Cjacobi,b.toarray())
    
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
    
    Morozov=0.7*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
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
    
    xk = xk.T.toarray()[0] # conversion du format sparse vers un array classique exploitable
    if (k == iter):
        print("Maximum iteration has been reached")
    print("Iterations number : "+str(k))
    print("========= ITERATION TIME ========")
    print(end-beg) 
 
    return(xk,xnaif)


def print_CCD(order,Y,beg,end):
    """
    Fonction qui print le CCD sous forme 2D pour un ordre et un vecteur de pixels donne.
    Inputs :
    - order : numero de l'ordre a tracer
    - Y : vecteur CCD a tracer
    """
    y0 = ord_centres[order]
    delta_courb_order = delta_courb[order]
    epaisseur = ep_voie[order]
    
    error = 0
    
    # On va remplir M avec l'odre selectionne. 
    for x in range(beg,end):
        # x est l'asbcisse de la colonne a afficher le long de l'ordre donne
        for i in range(epaisseur):
            # i est le numero en ordonee du pixel dans la colonne d'abscisse j
            # on calcule y, son indice d'ordonnee dans la matrice globale
            
            y = int(round(y0 + i + delta_courb_order[x]))
            #M[x,y] = Y[x*epaisseur+i]
            try :
                M[x,y] = Y[x*epaisseur+i]
            except :
                error += 1
                break
    plt.imshow(M,aspect='auto')
    plt.show()
    print("Nombre d'erreurs dans le remplissage du CCD : "+str(error))

    
    
# 
# 
# # Code pour multiplier le spectre choisi par A et obtenir le vecteur CCD correspondant
# ancho_str = '0_flat_edge'
# CCD_brut_all = []
# CCD_reduced_all = []
# 
# 
# for orden in range(10,20):
#     
#     
#     #print("Ordre : "+str(orden))
#     o0=interorden[orden]
#     o1=interorden[orden+1]
#     ov=interpol[orden]
#     
#     fout=open("extraction_methods/Extraction_files/REF/FP_map_new.pkl","r")
#     a=pickle.load(fout)
#     FPpics=a['picos']
#     fout.close()
#     
#     picos=FPpics[orden]
#     
#     Y=np.zeros(CCDsize*(ov+1))
#     YFlat = np.zeros(CCDsize*(ov+1))
#     arturos=np.shape(mapa)[1]
#     Souv=np.zeros((arturos,ov+1))
#     
#     nombre='extraction_methods/Extraction_files/REF/Amatrix_order'+str(orden)+'_lane1_ancho'+str(ancho_str)+'.pkl'
#     f=open(nombre,'r')
#     matriz=pickle.load(f)
#     Acompleta=matriz['Amatrix']#.toarray()
#     f.close()
#     jini,jfin=matriz['Limits']
#     ancho=matriz['Width']
#     lambda_max=0.33 #matriz['Lambda_max']
#     f.close()
# 
#     
#     cual0=int(np.floor(picos[jini+2]))
#     cualf=int(np.floor(picos[jfin-2]))
#     
#     
#     Yccd = np.zeros(order_len*ov)
#     Yccd_flat = np.zeros(order_len*ov)
#     for j in np.arange(0,arturos):
#         desp=mapa[orden,j]
#         B=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
#         BFlat = img2[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
#         Y[(ov+1)*j:(ov+1)*(j+1)]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
#         YFlat[(ov+1)*j:(ov+1)*(j+1)]=img2[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
#         Souv[j,:]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp)+1)]
#         Yccd[ov*j:ov*(j+1)]=img[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp))]
#         Yccd_flat[ov*j:ov*(j+1)]=img2[j,int(o0+np.floor(desp)):int(o0+ov+np.floor(desp))]
#     
#     A=Acompleta[(cual0-ancho)*(ov+1):(cualf+ancho+1)*(ov+1),cual0:cualf]#[(cual0-ancho)*(ov+1):cualf*(ov+1),cual0:cualf]
#     #Y2=Y[int((cual0-ancho)*(ov+1)):int(cualf*(ov+1))]
#     Y2=Y[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
#     Y2Flat=YFlat[int((cual0-ancho)*(ov+1)):int((cualf+ancho+1)*(ov+1))]
#     #Y3 = sparse.csr_matrix(Y2)
#     #Y3Flat = sparse.csr_matrix(Y2Flat)
#     
#     Yccd = Yccd/np.max(Yccd)
#     Yccd_flat = Yccd_flat/np.max(Yccd_flat)
#     Yccd2 = Yccd / Yccd_flat
#     #plt.figure(1)
#     #print_CCD(orden,Yccd)
#     # plt.figure(2)
#     # print_CCD(orden,Yccd_flat)
#     # plt.figure(3)
#     # print(np.shape(Yccd2))
#     # print_CCD(orden,Yccd2,cual0,cualf)
#     
#     # f = open("extraction_methods/Extraction_files/reduced/reduced_thar_"+str(orden)+"_ancho"+str(ancho_str)+".pkl",'r')
#     # thar_reduced = np.array(pickle.load(f))
#     # f.close()
#     
#     spectrum_order = spectrum[orden*order_len : (orden+1)*order_len]
#     spectrum_order = spectrum_order[cual0:cualf]
#     spectrum_order /= np.max(spectrum_order)
#     spectrum_order = sparse.csr_matrix(spectrum_order)
#     
#     # CCD reduit : on multiplie le spectre reduit par A pour retrouver le vecteur pixe CCD    
#     CCD_reversed = A.dot(spectrum_order.T)
#     #plt.figure(2)
#     CCD_reversed = CCD_reversed.toarray()
#     # print(np.shape(CCD_reversed))
#     # print_CCD(orden,CCD_reversed,cual0,cualf)
#     # thar_reduced = sparse.csr_matrix(thar_reduced)
#     # ccd_reduced = A.dot(thar_reduced.T)
#     # ccd_reduced = ccd_reduced - np.min(ccd_reduced)
#     # ccd_reduced = ccd_reduced.toarray() / np.max(ccd_reduced.toarray())
#     # ccd_brut = Y2 - np.min(Y2)
#     # ccd_brut = ccd_brut / np.max(ccd_brut)
#     # CCD_brut_all.extend(ccd_brut)
#     # CCD_reduced_all.extend(ccd_reduced)
# 

    

# Code pour construire un CCD en 1D a partir du CCD en 2D

ep_voie = interpol # liste des epaisseurs de chaque voie 
ord_centres= interorden # liste des ordonnees du centre de chaque ordre (axe Y allant de 0 a 2098)
delta_courb = mapa # matrice des listes (par ordre) des delta_pixel de courbure de l'ordre par rapport à la colonne centrale

def build_Yccd(ccd2d,order):
    """ Cette fonction construit le vecteur ccd second membre a partir du fits ouvert et pour l'ordre donné.
    Inputs : 
    - ccd2d : image du fits
    - order : numero de l'ordre a traiter
    Outputs :
    - Yccd : vecteur ccd en 1D
    """

    y0 = ord_centres[order]
    delta_courb_order = delta_courb[order]
    epaisseur = ep_voie[order]
    picos=FPpics[order]
    beg = int(round(picos[0]))
    end = int(round(picos[-1]))
    Yccd = np.zeros(epaisseur*order_len)
    
    for x in range(beg,end):
        
        for i in range(epaisseur):
            
            y = int(round(y0 + delta_courb_order[x] + i))
            
            Yccd[x*epaisseur+i] = ccd2d[x,y]
    return(Yccd)
    

def build_Yccd_1row(ccd2d,order):
    """ Cette fonction construit le vecteur ccd second membre a partir du fits ouvert et pour l'ordre donné, en sommant toutes les lignes de la voie pour en obtenir qu'une seule.
    Inputs : 
    - ccd2d : image du fits
    - order : numero de l'ordre a traiter
    Outputs :
    - Yccd : vecteur ccd en 1D
    """

    y0 = ord_centres[order]
    delta_courb_order = delta_courb[order]
    epaisseur = ep_voie[order]
    picos=FPpics[order]
    beg = int(round(picos[0]))
    end = int(round(picos[-1]))
    Yccd = np.zeros(order_len)
    
    for x in range(beg,end):
        
        Yccd[x] = np.sum([ccd2d[x,int(round(y0 + delta_courb_order[x] + i))] for i in range(epaisseur)])
        
    return(Yccd)

# Code pour remplir A 

fout=open("extraction_methods/Extraction_files/REF/FP_map_new.pkl","r")
a=pickle.load(fout)
FPpics=a['picos']
fout.close()


def fill_A(ccd2d,order,largeur):
    """ Cette fonction construit la matrice de convolution A pour un order donne, sans interpoler entre les raies de FP.
    Inputs : 
    - ccd2d la matrice originelle du fits
    - order : l'ordre considéré
    - largeur : demi largeur de la fenetre de construction de A
    Outputs:
    - A : la matrice construite sans interpolation pour l'odre donne
    """
    y0 = ord_centres[order]
    delta_courb_order = delta_courb[order]
    epaisseur = ep_voie[order]
    picos=FPpics[order]
    Yccd = build_Yccd(ccd2d,order)
    M = np.zeros((order_len,len(Yccd)))
    M = sparse.lil_matrix(M)
    for k in range(len(picos)):
        x = int(round(picos[k]))
        beg = (x-largeur)*epaisseur
        end = (x+largeur)*epaisseur
        M.rows[x] = np.arange(beg,end)
        M.data[x] = Yccd[beg:end] - np.min(Yccd[beg:end])
    return(M)


def interpol_A(Abrut,order,largeur):
    """ Cette fonction prend une matrice A remplie uniquement aux pics du FP/Thar et interpole entre les pics pour compléter la "diagonale".
    Inputs :
    - Abrut : matrice issue de fill_A, en format sparse csr matrix
    - order : numero de l'ordre a interpoller
    - epaisseur : largeur de la fenetre consideree autour du pic de FP (voir fill_A)
    Outputs :
    - Ainterp : matrice remplie par interpollation 
    """
    # Il est plus pratique pour travailler en format sparse.lil_matrix d'utiliser A.T, on retransposera apres
    Ainterp = np.zeros(np.shape(Abrut))
    Ainterp = sparse.lil_matrix(Ainterp)
    # On va commencer par retrouver les colonnes deja remplies de la matrice pour pouvoir ebsuite interpoler entre elles. Par construction, on peut reutiliser la liste des pics de FP : picos
    picos=FPpics[order]
    epaisseur = ep_voie[order]
    for k in range(len(picos)-1):
        xk1 = int(round(picos[k]))
        xk2 = int(round(picos[k+1]))
        Axk1 = Abrut.toarray()[xk1,(xk1-largeur)*epaisseur:(xk1+largeur)*epaisseur]
        Axk2 = Abrut.toarray()[xk2,(xk2-largeur)*epaisseur:(xk2+largeur)*epaisseur]
        for x in range(xk1,xk2):
            alpha = ( x - xk1) / (xk2 - xk1)
            beg = (x-largeur)*epaisseur 
            end = (x+largeur)*epaisseur
            Ainterp.rows[x] = np.arange(beg,end)
            Ainterp.data[x] = alpha*Axk2 + (1-alpha)*Axk1
    Ainterp = Ainterp.tocsr().T
    return(Ainterp)
    

def fill_A_1row(ccd2d,order,largeur):
    """ Cette fonction construit la matrice de convolution A pour un order donne, sans interpoler entre les raies de FP, en utilisant qu'une ligne pour chaque ordre (on va sommer toutes les lignes).
    Inputs : 
    - ccd2d la matrice originelle du fits
    - order : l'ordre considéré
    - largeur : demi largeur de la fenetre de construction de A
    Outputs:
    - A : la matrice construite sans interpolation pour l'odre donne
    """
    y0 = ord_centres[order]
    delta_courb_order = delta_courb[order]
    picos=FPpics[order]
    Yccd = build_Yccd_1row(ccd2d,order)
    M = np.zeros((order_len,len(Yccd)))
    M = sparse.lil_matrix(M)
    for k in range(len(picos)):
        x = int(round(picos[k]))
        beg = x-largeur
        end = x+largeur
        M.rows[x] = np.arange(beg,end)
        M.data[x] = Yccd[beg:end] - np.min(Yccd[beg:end])
    return(M)
    
def interpol_A_1row(Abrut,order,largeur):
    """ Cette fonction prend une matrice A remplie uniquement aux pics du FP/Thar et interpole entre les pics pour compléter la "diagonale".
    Inputs :
    - Abrut : matrice issue de fill_A, en format sparse csr matrix
    - order : numero de l'ordre a interpoller
    - epaisseur : largeur de la fenetre consideree autour du pic de FP (voir fill_A)
    Outputs :
    - Ainterp : matrice remplie par interpollation 
    """
    # Il est plus pratique pour travailler en format sparse.lil_matrix d'utiliser A.T, on retransposera apres
    Ainterp = np.zeros(np.shape(Abrut))
    Ainterp = sparse.lil_matrix(Ainterp)
    # On va commencer par retrouver les colonnes deja remplies de la matrice pour pouvoir ebsuite interpoler entre elles. Par construction, on peut reutiliser la liste des pics de FP : picos
    picos=FPpics[order]
    for k in range(len(picos)-1):
        xk1 = int(round(picos[k]))
        xk2 = int(round(picos[k+1]))
        Axk1 = Abrut.toarray()[xk1,xk1-largeur:xk1+largeur]
        Axk2 = Abrut.toarray()[xk2,xk2-largeur:xk2+largeur]
        for x in range(xk1,xk2):
            alpha = ( x - xk1) / (xk2 - xk1)
            beg = x-largeur
            end = x+largeur
            Ainterp.rows[x] = np.arange(beg,end)
            Ainterp.data[x] = alpha*Axk2 + (1-alpha)*Axk1
    Ainterp = Ainterp.tocsr().T
    return(Ainterp)

# imgFlat = imgFlat/np.max(imgFlat)
# img = img/imgFlat
# 
# order =15
# ini = int(round(FPpics[order][0]))
# end = int(round(FPpics[order][-1]))
# 
# Yccdflat_1row = build_Yccd_1row(imgFlat,order)[ini:end]
# Yccd1row = build_Yccd_1row(img,order)[ini:end]
# Yccd = build_Yccd(img,order)[ini:end]
# plt.figure(1)
# Yccd1rown = Yccd1row / np.max(Yccd1row)
# #plt.plot(Yccd,'black')
# plt.plot(Yccd1rown,'blue')
# plt.show()
# 
# A = fill_A_1row(imgFP,order,10)
# # plt.figure(3)
# # plt.imshow(A.toarray(),aspect='auto')
# A = interpol_A_1row(A,order,10)
# A = A.toarray()[ini:end,ini:end]
# A = sparse.csr_matrix(A)
# plt.figure(2)
# plt.imshow(A.toarray(),aspect='auto')
# plt.show()
# 
# Yccd1row = sparse.csr_matrix(Yccd1row)
# Yccdflat_1row = sparse.csr_matrix(Yccdflat_1row)
# Thar_sp1row = conjugate_grad_prec_sparse(A,Yccd1row)[0]
# Flat_sp1row = conjugate_grad_prec_sparse(A,Yccdflat_1row)[0]
# #Thar_sp1row = Thar_sp1row/np.max(Thar_sp1row)
# #Flat_sp1row = Flat_sp1row/np.max(Flat_sp1row)
# plt.figure(1)
# #plt.plot(Thar_sp1row,'red')
# #plt.plot(Flat_sp1row,'black')
# sp = Thar_sp1row/(Flat_sp1row+1e-12)
# #sp = [ k[0] for k in sp]
# for k in range(len(sp)) :
#     if np.isnan(sp[k]) :
#         sp[k] = 0
# sp = sp / np.max(sp)
# plt.plot(sp,'purple')
# plt.show()



# Construction et enregistrement de la matrice A
ancho = 8
spectra = []
flats = []
for order in range(14,15):
    if True :
        print("###### ORDER : "+str(order)+" ######")
        t0 = time.time()
        A = fill_A(imgFP,order,ancho)
        print("Temps de construction de A : ",time.time()-t0)
        
        plt.figure(1)
        B = A.toarray()[2000:3000,:]
        plt.imshow(B,aspect='auto')
        plt.show()
        
        t1 = time.time()
        Ainterp = interpol_A(A,order,ancho)
        print("Temps d'interpolation de A :",time.time()-t1)
        
        plt.figure(2)
        C = Ainterp.toarray()[:,2000:3000]
        plt.imshow(C,aspect='auto')
        plt.show()
        
       
        picos=FPpics[order]
        jini=5
        jfin=len(picos)-5
        f = 'extraction_methods/Extraction_files/REF/Amatrix_order'+str(order)+'_lane1_ancho'+str(ancho)+'_new_test'
        f = open(f,'w')
        pickle.dump({'Amatrix':Ainterp,'Limits':(jini,jfin),'Width':10,'Lambda_max':1},f)
        f.close()
        
        

        # test de reduction en utilisant A construite au dessus
        plt.figure(3)
        YccdThar = build_Yccd(img,order)
        plt.plot(YccdThar,'red')
        plt.show()
        YccdThar = sparse.csr_matrix(YccdThar)
        YccdFlat = build_Yccd(imgFlat,order)
        for k in range(len(YccdFlat)) :
            if YccdFlat[k] != 0 :
                YccdFlat[k] = YccdFlat[k] - 32768.0  # On retire l'offset decrit dans les headers du fit : B0 = 2**15
        plt.plot(YccdFlat,'black')
        plt.show()
        YccdFlat = sparse.csr_matrix(YccdFlat)
        
        f = 'extraction_methods/Extraction_files/REF/Amatrix_order'+str(order)+'_lane1_ancho'+str(ancho)+'_new_test'
        f = open(f,'r')
        Adata = pickle.load(f)
        f.close()
        A = Adata['Amatrix']
        
        Thar_spectrum = conjugate_grad_prec_sparse(A,YccdThar)[1]
        #Thar_spectrum = Thar_spectrum/np.max(Thar_spectrum)
        Flat_spectrum = conjugate_grad_prec_sparse(A,YccdFlat)[1]
        #Flat_spectrum = Flat_spectrum/np.max(Flat_spectrum)
        #Thar_normalized_spectrum = Thar_spectrum / Flat_spectrum
        plt.figure(4)
        plt.plot(Thar_spectrum,'blue')
        plt.plot(Flat_spectrum,'black')
        #plt.plot(Thar_normalized_spectrum,'red')
        plt.show()
        
        
        Thar_spectrum2 = conjugate_grad_prec_sparse(A,YccdThar)[0]
        Flat_spectrum2 = conjugate_grad_prec_sparse(A,YccdFlat)[0]
        #Thar_normalized_spectrum = Thar_spectrum / Flat_spectrum
        plt.figure(5)
        plt.plot(Thar_spectrum2,'blue')
        plt.plot(Flat_spectrum2,'black')
        #plt.plot(Thar_normalized_spectrum,'red')
        plt.show()
        
        spectra.extend(Thar_spectrum)
        flats.extend(Flat_spectrum)
    else :
        print("Problem with order "+str(order))
plt.figure(6)
plt.plot(flats,'black')
plt.plot(spectra,'red')
plt.show()
    
