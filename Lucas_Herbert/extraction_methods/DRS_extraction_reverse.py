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

datafileSimu = datafile + "Simulation_Antony_FP.fts"
datafileFlat = datafile + "Simulation_Antony_Flat.fts"
datafileStar = datafile + "13mar18/Narval_20180314_003846_st0.fts"
datafileThar = datafile + "13mar18/Narval_20180313_181059_th0.fts"
datafileFP = datafile + "Simulation_Antony_FP.fts"
f=pyfits.open(datafileThar)
img=f[0].data
f.close()
f = pyfits.open(datafileFlat)
imgFlat=f[0].data
f.close()
f=pyfits.open(datafileFP)
imgFP=f[0].data
f.close()
f=pyfits.open(datafileSimu)
imgSimu=f[0].data
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


""" 
Ce code retrouve la localisation de chaque ordre dans l'image CCD 2D et la renvoie sous forme de liste.
"""
# 
# # On commence par retrouver la position des ordres en trouvant l'ordonnee du pixel de la colonne centrale de chaque ordre grâce au code qui suit.
# 
# x_central = int(order_len/2)
# 
# colonne_centrale = [] # liste des pixels de la colonne centrale
# 
# for y in range(np.shape(imgFlat)[1]) :
#     colonne_centrale.append(imgFlat[x_central,y])
# 
# ordre_centres = [58,89,119,150,181,212,244,276,310,344,378,414,451,489,528,568,609,650,695,740,786,834,883,934,987,1041,1096,1154,1213,1274,1338,1403,1470,1540,1612,1687,1763,1843,1926,2012]
# epaisseurs = [15 for i in range(len(ordre_centres))]
# 
# f = open("extraction_methods/Extraction_files/REF/reference_simu",'w')
# pickle.dump(ordre_centres,f)
# f.close()
# 
# f = open("extraction_methods/Extraction_files/REF/epaisseurs",'w')
# pickle.dump(epaisseurs,f)
# f.close()
# 
# 
# plt.figure(1)
# plt.plot(colonne_centrale)
# for i in range(len(ordre_centres)) :
#     plt.axvline(ordre_centres[i],color='red')
#     y = ordre_centres[i] + epaisseurs[i]
#     plt.axvline(y,color='green')
# plt.show()
# 
# # Puis il nous faut trouver la forme de chaque ordre : le decalage en y suivant l'indice x. Pour cela on va utilise la colonne centrale du flat que l'on va intercorreler avec la meme colonne decalee en indices (abscisses). Le shift obtenu sera notre decalage. Pour chaque ordre on retiendra donc le decalage pour chaque indice nous donnant une liste de decalages de longueur order_len.

def correlator(a,b):
    
    """
    We can also find a big shift with the help of a cross correlation function. This function will implement the finding of a shift thanks to the cross correlation.
    Inputs :
    - a : array or list
    - b : array or list (a with a shift which we want to comute)
    Output :
    - ind : the shift between a and b
    """
    plt.plot(a,'black')
    plt.plot(b,'red')
    plt.show()
    
    r = np.correlate(a,b,"full")
    
    max = r[0]
    ind = 0
    
    for k in range(len(r)):
        
        if r[k] > max :
            max= r[k]
            ind = k
            
            
    return( - (-len(a) + ind + 1) )
    
    
def shifter(Sr, shift):
    """
    The following function is used to test chelli algorithm. It takes a rerefrence spectrum and shift it from the given delta.
    Inputs :
    - Sr : reference spectrum
    - delta : shift to compute
    Outputs :
    - S : shifted spectrum
    """
    fSr = np.fft.rfft(Sr)
    iList = np.arange(len(fSr))
    k = -2j*np.pi*iList*1.0/len(Sr)*shift
    fS = np.exp(k)*fSr
    S = np.fft.irfft(fS)
    return(S)
    


def chelli_shift(Sr, S, step, v=1):
    """
        Function to derive the Doppler shift between two spectrales rays, based on Chelli(2000). The arguments are :
            Sr : spectrum of reference
            S : shift spectrum
            step : spatial discretization step
            v : initial value of the shift
        Return the shift in pixels, the error
    """
    # I hat
    I = np.fft.rfft(Sr)*np.conjugate(np.fft.rfft(S))
    lenI = len(I)

    # List of Q, the quantity to minimize
    Qinit = 0.
    Qfinal = 1.

    cpt = 0
    cpt_max = 2000
    # Stop condition on Q and the cpt_max
    condition = (abs(Qfinal-Qinit) > 1e-20) and cpt < cpt_max 
    # Qi,Qf = [Qinit],[Qfinal]
    iList = np.arange(lenI)
    k = -2j*np.pi*iList*v/(2*step*lenI)

    C = I*np.exp(k)

    while (condition):
        # Sigma2 represents the variance of Im(C). We use a constant value,
        # this one can be modified to matches Chelli's method
        Sigma2 = 0.5 + np.abs(C)
        # Derive of Deltav to add to v
        Num = np.sum(iList*C.real*C.imag / Sigma2**2)
        Den = np.sum((iList*C.real)**2 / Sigma2**2)

        DeltaV = Num / Den
        sigma_noise = 1/np.sum(iList**2 * np.abs(C))
        v += DeltaV

        k = -2j * np.pi * iList * v /(2*step*lenI)
        C = I * np.exp(k)
        # Update of Q
        Qinit = Qfinal
        Qfinal = np.sum(C.imag**2 / Sigma2)

        while (Qfinal >= Qinit and abs(DeltaV) > 1e-20):
            v -= DeltaV
            DeltaV *= 0.1
            v += DeltaV  # we change the step for V if Q increase
            k = -2j*np.pi*iList*v/(2*step*lenI)
            C = I*np.exp(k)

            Qtest = np.sum(C.imag**2 / Sigma2)
            if (Qtest > Qfinal):
                v -= DeltaV
                DeltaV *= 10
                v += DeltaV
                break
            else:
                Qfinal = Qtest

        # Update of the conditionJe joins à ce mail un CV et une lettre de motivation expliquant ma situation et mes disponibilités.
        condition = (abs(Qfinal-Qinit) > 1e-20) and cpt < cpt_max

        cpt += 1

    if (cpt == cpt_max):             # Case of non convergence
        return (-1, 0, 0)
        print(" /!\ Chelli has not converged! ")
    return (v, Qfinal, sigma_noise)  # Case of convergence


def find_shift(Sr,S):

    """
    The following function will find the shift between two spectrum but only with an error of 1 pixel : it is a way to find a big shift. Since chelli shift is only working on less than one pixel shifts we need to first find the shift with a precision of one pixel, then apply chelli shift. This function is made for this objective.
    Input :
    - Sr : reference spectrum.
    - S : shifted spectrum.
    Output :
    - shift : the shift in pixel with a precision of one pixel.
    """    
    shift = 0
    shiftj = 0
    # Initiating the least square indicator with an initial value
    Q = np.sum( [ (Sr[i]-S[i])**2 for i in range(len(Sr)) ] )
    
    # We are going to improve the accuracy step by step, dividing it by 10 in each loop... We start between -110 and 110, then -11 to 11, the -1.1 to 1.1 so we don't miss anything.
    for precision in range(1,6):
        
        jlist = [ i*10**-(precision-1) for i in range(-11,12) ]
        for j in jlist :
            
            Srj = shifter(Sr,shift+j)
            # Computing the current least square indicator for the current j shift    
            Qj = np.sum( [ (Srj[i]-S[i])**2 for i in range(len(Srj)) ] )
            # Updating only if the lq indicator is better ad recording the best found shift until here...
            if  Qj < Q :
                Q = Qj
                shiftj = j
                
        # updating the shift by adding the new shift we have just computed    
        shift += shiftj

    return(shift)

"""
Code ayant deja ete execute calculant la position des ordres, voies et pics.
"""
# Maintenant, on peut trouver le decalage etre deux colonnes du ccd grace aux fonctions ci dessus on va donc pouvoir retrouver la courbure d'un ordre. 

# pour chaque ordre on va selectionner le morceau de colonne lui correspondant sur la colonne centrale et trouver le shift entre ce morceau de colonne et celui de chaque colonne du ccd pour voir comment chaque colonne est decalee. Pour chaque colonne, pour chaque ordre, on retiendra ce shift nous donnant une liste de 4612 shifts pour chaque ordre, sur laquelle sera fitte un polynome de degre 2 ou 3 qui nous donnera la forme precise de la courbure de l'ordre. Ces resultats permettront ensuite de retrouver nos ordres et poursuivre le travail en construisant nos matrices, etc.
# 
# for order in range(0,32) :
#     # code pour un ordre : par exemple l'ordre 10 :
#     #order = 10
#     ref = ordre_centres[order]
#     ep = epaisseurs[order]
#     marge = 20
#     colonne_centrale_ordre = colonne_centrale[ref - marge : ref + 2*ep + marge]
#     
#     courbure = [] # liste des shifts 
#     # Detection des shifts a gauche du pixel central
#     shift = 0
#     for x in range(0,-int(order_len/2),-1):
#         colonne_x = imgFlat[x_central+x, int(shift) + ref-marge : int(shift) + ref +2*ep+marge]
#         plt.figure(2)
#         plt.plot(colonne_centrale_ordre,'black')
#         plt.plot(colonne_x,'red')
#         plt.show()
#         if ( np.abs(find_shift(colonne_centrale_ordre,colonne_x)) < ep ) :
#             shift += find_shift(colonne_centrale_ordre,colonne_x) 
#         courbure.insert(0,shift)
#     
#     # Detection des shifts a droite du pixel central
#     shift = 0
#     for x in range(0,int(order_len/2)):
#         colonne_x = imgFlat[x_central+x,int(shift) + ref-marge : int(shift) + ref +2*ep+marge]
#         plt.figure(2)
#         plt.plot(colonne_centrale_ordre,'black')
#         plt.plot(colonne_x,'red')
#         plt.show()
#         if ( np.abs(find_shift(colonne_centrale_ordre,colonne_x)) <= ep ) :
#             shift += find_shift(colonne_centrale_ordre,colonne_x) 
#         courbure.append(shift)
#     plt.figure(3)
#     plt.plot(courbure,'black')
#     plt.show()
#     
#     # Nous avons la liste des shifts pout un ordre donne : on va interpoler avec un polynome 
#     
#     indices = [i for i in range(4612) ]
#     
#     coeffs = np.polyfit(indices,courbure,3)
#     
#     pol_courb = np.poly1d(coeffs)
#     courbure_interp = pol_courb(indices)
#     plt.figure(3)
#     plt.plot(courbure_interp,'red')
#     plt.show()
#     
#     # Cette liste est constituee de floats, or on va l'utiliser pour se balader dans une matrice, donc on doit tout transformer en indices entiers.
#     courbure_ordre = [ int(round(shift)) for shift in courbure_interp ]
#     
#     # Il nous reste a enregistrer tous ces parametres
#     f = open("extraction_methods/Extraction_files/REF/courbure_ordre_"+str(order),'w')
#     pickle.dump(courbure_ordre,f)
#     f.close()
#     print("Order "+str(order)+" : recorded!")
#     
#     
#     # Il nous faut maintenant faire le mapping des psf de fabry perot sur le CCD FP (imgFP). On sait ou chercher les ordres donc pour chaque ordre on va retouver sommer toutes les lignes et detecter automatiquement tous les pics d'intensite qui correspondent aux raies du fabry perot. 
#     
#     y0 = ref
#     epaisseur = epaisseurs[order]
#     
#     Y_fp = []
#     
#     for x in range(4612) :
#             
#         value = np.sum([ imgFP[x, ref + courbure_ordre[x] + i] for i in range(epaisseur) ])
#         
#         Y_fp.insert(x,value)
#     
#     plt.figure(4)
#     plt.plot(Y_fp,'green')
#     plt.show()
#     
#     # Nous devons maintenant detecter les pics de fp et les localiser dans une liste qu'on reutilisera plus tard.
#     
#     pics_FP = [] #liste des indices des pics de fp
#     
#     # On va detecter ces pics mathematiquement en parcourant la liste
#     
#     for i in range(10,4602) :
#         
#         local_mean = np.mean([Y_fp[k] for k in range(i-10,i+10) ])
#         
#         if ( Y_fp[i-1] < Y_fp[i] and Y_fp[i+1] < Y_fp[i] and Y_fp[i] > local_mean ) :
#             pics_FP.append(i)
#     
#     plt.figure(4)
#     for pic in pics_FP :
#         plt.axvline(pic)
#     plt.show()
#     
#     # Enregistrement de notre liste de pics :
#     f = open("extraction_methods/Extraction_files/REF/pics_FP_order_"+str(order),'w')
#     pickle.dump(pics_FP,f)
#     f.close()

# courbures_simu = []
# pics_FP_simu = []
# 
# for order in range(0,32):
#     
#     fcourb = open("extraction_methods/Extraction_files/REF/courbure_ordre_"+str(order),'r')
#     courbure_ordre = pickle.load(fcourb)
#     fcourb.close()
#     courbures_simu.append(courbure_ordre)
# 
#     fpics = open("extraction_methods/Extraction_files/REF/pics_FP_order_"+str(order),'r')
#     pics_ordre = pickle.load(fpics)
#     fpics.close()
#     pics_FP_simu.append(pics_ordre)
# 
# f = open("extraction_methods/Extraction_files/REF/courbures_simu",'w')
# pickle.dump(courbures_simu,f)
# f.close()
# 
# f = open("extraction_methods/Extraction_files/REF/pics_FP_simu",'w')
# pickle.dump(pics_FP_simu,f)
# f.close()


f = open("extraction_methods/Extraction_files/REF/courbures_simu",'r')
courbures_simu = pickle.load(f)
f.close()

f = open("extraction_methods/Extraction_files/REF/pics_FP_simu",'r')
pics_FP_simu = pickle.load(f)
f.close()

f = open("extraction_methods/Extraction_files/REF/epaisseurs",'r')
epaisseurs = pickle.load(f)
f.close()

f = open("extraction_methods/Extraction_files/REF/reference_simu",'r')
ref_ordres = pickle.load(f)
f.close()







# Code  pour afficher le vecteur CCD en matrice 2D
# 
# ep_voie = epaisseurs #interpol # liste des epaisseurs de chaque voie 
# ord_centres= ordre_centres#interorden # liste des ordonnees du centre de chaque ordre (axe Y allant de 0 a 2098)
# delta_courb = mapa # matrice des listes (par ordre) des delta_pixel de courbure de l'ordre par rapport à la colonne centrale

ep_voie = epaisseurs
ord_centres = ref_ordres
delta_courb = courbures_simu
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

    
    

# Code pour construire un CCD en 1D a partir du CCD en 2D

# ep_voie = interpol # liste des epaisseurs de chaque voie 
# ord_centres= interorden # liste des ordonnees du centre de chaque ordre (axe Y allant de 0 a 2098)
# delta_courb = mapa # matrice des listes (par ordre) des delta_pixel de courbure de l'ordre par rapport à la colonne centrale

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

# fout=open("extraction_methods/Extraction_files/REF/FP_map_new.pkl","r")
# a=pickle.load(fout)
# FPpics=a['picos']
# fout.close()
FPpics = pics_FP_simu

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
    #M = np.zeros((order_len,len(Yccd)))
    M = sparse.lil_matrix((order_len,len(Yccd)))
    for k in range(len(picos)):
        x = int(round(picos[k]))
        beg = (x-largeur)*epaisseur
        end = (x+largeur)*epaisseur
        M.rows[x] = np.arange(beg,end)
        M.data[x] = Yccd[beg:end] #- np.min(Yccd[beg:end])
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
    #M = np.zeros((order_len,len(Yccd)))
    M = sparse.lil_matrix((order_len,len(Yccd)))
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





# Construction et enregistrement de la matrice A
ancho = 3
for order in range(10,11):
    if True :
        print("###### ORDER : "+str(order)+" ######")
        t0 = time.time()
        A = fill_A(imgSimu,order,ancho)
        print("Temps de construction de A : ",time.time()-t0)
         
        # plt.figure(1)
        # B = A.toarray()[2000:3000,:]
        # plt.imshow(B,aspect='auto')
        # plt.show()
        # 
        t1 = time.time()
        Ainterp = interpol_A(A,order,ancho)
        print("Temps d'interpolation de A :",time.time()-t1)
        
        # plt.figure(2)
        # C = Ainterp.toarray()[:,2000:3000]
        # plt.imshow(C,aspect='auto')
        # plt.show()
        # 
        
        picos=FPpics[order]
        jini=5
        jfin=len(picos)-5
        f = 'extraction_methods/Extraction_files/REF/Amatrix_order'+str(order)+'_lane1_ancho'+str(ancho)+'_test_simu'
        f = open(f,'w')
        pickle.dump({'Amatrix':Ainterp,'Limits':(jini,jfin),'Width':10,'Lambda_max':1},f)
        f.close()
        
        

        # test de reduction en utilisant A construite au dessus
        plt.figure(3)
        YccdSimu = build_Yccd(imgSimu,order)
        plt.plot(YccdSimu,'red')
        plt.show()
        YccdSimu = sparse.csr_matrix(YccdSimu)
        
        f = 'extraction_methods/Extraction_files/REF/Amatrix_order'+str(order)+'_lane1_ancho'+str(ancho)+'_test_simu'
        f = open(f,'r')
        Adata = pickle.load(f)
        f.close()
        A = Adata['Amatrix']
        
        Simu_spectrum = conjugate_grad_prec_sparse(A,YccdSimu)[1]
        #Simu_spectrum = Simu_spectrum/np.max(Simu_spectrum)
        plt.figure(4)
        plt.plot(Simu_spectrum,'red')
        plt.show()
        
        
        # Simu_spectrum2 = conjugate_grad_prec_sparse(A,YccdSimu)[0]
        # plt.figure(5)
        # plt.plot(Simu_spectrum2,'blue')
        # plt.show()

    else :
        print("Problem with order "+str(order))

