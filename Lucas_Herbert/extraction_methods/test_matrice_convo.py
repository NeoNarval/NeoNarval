#! /usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
from scipy.optimize import curve_fit
from scipy import sparse
import time as time
import numpy.ma as ma
import scipy.optimize
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

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
   
   
def clean() :
    for i in range(1000):
        plt.close(i)
        
def plotter2D(Z,n) :
    """
    Plot Z in 3D in the plt.figure(n). 
    """
    fig = plt.figure(n)
    ax = Axes3D(fig)
    X= np.arange(np.shape(Z)[1])
    Y = np.arange(np.shape(Z)[0])
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')
    plt.show()
    

def gaussian_fit(lambdas,data):
    """
    This function finds the gaussian function which is the closest from our data. It's a fitting algorithm which uses a least squares method. It takes an an input the necessary data for the fit : the wavelengths, lambdas, and their intensities (or whatever we need to fit), data. It returns the centre and width of the fitted gaussian and also the report of the fit (saying if it is good enough or not). It also plots the fitted gaussians and their centers so that we can "visualy check" if the fit is good enough.
    Inputs :
    - lambdas : list of wavelengths
    - data : list of intensities (or whatever) associated to the wavelengths list.
    Outputs :
    - lambda_centre : float, wavelength of the center of the gaussian which has been fitted to our data and lambdas. It represents the actual center of the spike.
    - lambda_width : float, width of the fitted gaussian. It represents the width of the spike.
    - report : string, the report contains all the informations about the fit : the chi_square, the nuber of iterations before converging to the solution, etc.
    """

    y = np.copy(data) - np.min(data)
    X = lambdas
    def gaussian(x,cen,amp,wid):
        if (wid != 0) :
            return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
        else : 
            wid = 1e-10
            return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
    
    # computing the initial gaussian in wavemengths (before the optimization of the fit)
    naive_center = float(np.sum(lambdas*y))/np.sum(y)
    naive_width = np.sqrt(abs((np.sum((lambdas-naive_center)**2*y)/np.sum(y))))
    naive_ampl = np.max(y)
    naive_gaussian = [ gaussian(x,naive_center,naive_ampl,naive_width) for x in lambdas ]
    
    
    lambda_centre = naive_center
    lambda_width = naive_width
    

    try :
        # we use the lmfit algorithm to improve our fit's precision
        gaussian_model = lmfit.Model(gaussian)
        params = gaussian_model.make_params(cen=naive_center,amp=naive_ampl,wid=naive_width)
        result = gaussian_model.fit(y,params,x=X) 
        # printing the best gaussian fit
        best_gaussian_fit = result.best_fit
        best_cen = result.best_values['cen']
        best_wid = result.best_values['wid']
        best_amp = result.best_values['amp']
        best_params_gaussian = [ gaussian(x,best_cen,best_amp,best_wid) + np.min(data) for x in lambdas ]
        plt.plot(lambdas, best_params_gaussian, 'b--', color='purple')
        plt.plot(best_cen,best_amp,'.',color='purple')
        plt.show()
        #computed_centre = float(np.sum(X*best_gaussian_fit))/np.sum(best_gaussian_fit) 
        #plt.plot(lambdas, best_gaussian_fit, 'b--' ,color='purple')
        #plt.axvline(computed_centre, color='purple')
        lambda_centre = best_cen
        lambda_width = best_wid
        # we need the report data to improve our understanding of the results
        report = result.fit_report()
         
        
    except : 
        report = "Computation failed for this spike : default data = naive fit"
        pass
        
    return(lambda_centre,lambda_width,report)


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
    # print( np.shape(sigma))
    # print(np.shape(A))
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
    M = np.eye(np.shape(C)[0]) #Cjacobi #np.eye(np.shape(C)[0])
    #M = M.T.dot(M)
    #M = np.linalg.inv(M)
    M = sparse.csr_matrix(M)
    #print(np.linalg.cond( (M.dot(C)).toarray() ))
    
    #================================================
    # Facteur de relaxation? 
    omega = 0.5
    #================================================
    
    # Initialisation :
    
    x0 = sparse.csr_matrix( np.zeros(np.shape(C)[0]) ).T    #sparse.csr_matrix(xnaif)    #sparse.csr_matrix( np.zeros(np.shape(C)[0]) ).T
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
    
    Morozov=0.0001*np.median(1 / np.sqrt( Y.toarray() )) # critere d'arret de Morozov
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
    # xkn = xk/np.max(xk)	# on normalise xk
    # plt.figure(1)	
    # plt.plot(xkn,color='green')
    # plt.plot(xnaifn,color='black')
    # plt.show()
    return(xk,xnaif)


# on fait de la place
#clean()



# parametres 

#  generation des ccd avant ajout psf

ccd1 = np.zeros((20,100))
ccd2 = np.zeros((20,100))
ccdFP = np.zeros((20,100))
ccdth = np.zeros((20,100))
ccdth_red = np.zeros((20,100))




cen_psf_cons = 50

# allons chercher une psf fabry perot sur un fits :
f = pyfits.open("extraction_methods/Extraction_files/FP_test.fts")
img_FP = f[0].data
local_FP = img_FP[2632:2677,546:562] - np.min(img_FP[2632:2677,546:562])
local_FP = local_FP 

# allons chercher une raie ThAr sur un fits :
f = pyfits.open("extraction_methods/Extraction_files/th0_test")
img_th = f[0].data
local_th = img_th[1782:1790,791:807] - np.min(img_th[1782:1790,791:807])
local_th = local_th 


def ajout_img_ccd(ccd,posX,posY,img) :
    larg = np.shape(img)[0]
    haut  = np.shape(img)[1]
    for x in range(posX - int(larg/2), posX + int(larg/2) ):
        for y in range(posY-int(haut/2),posY+int(haut/2)):
            ccd[y,x] += img[x - posX- int(larg/2) ,y - posY - int(haut/2)]
            
ajout_img_ccd(ccdFP,50,10,local_FP)
plotter2D(ccdFP,10)



ajout_img_ccd(ccdth,10,10,local_th)
ajout_img_ccd(ccdth,90,10,local_th)
plotter2D(ccdth,70)

ajout_img_ccd(ccdth_red,50,10,local_FP)
ajout_img_ccd(ccdth_red,53,10,local_th)
plotter2D(ccdth_red,700)


# version naive de la reduction
somme_yFP = []
for i in range(100) :
    somme_yFP.append(np.mean([ccdFP[j,i] for j in range(20)]))
# plt.figure(100)
# plt.title("Somme sur Y des intensites ccdFP")
# plt.plot(np.arange(100),somme_yFP)
# plt.show()


somme_yth = []
for i in range(100) :
    somme_yth.append(np.mean([ccdth[j,i] for j in range(20)]))
plt.figure(71)
plt.title("Somme sur Y des intensites ccdth")
plt.plot(np.arange(100),somme_yth)
plt.show()

# somme yth_red reduit par la matrice Ath
somme_yth_red = []
for i in range(100) :
    somme_yth_red.append(ccdth_red[10,i])
plt.figure(711)
plt.title("Ligne 10 du ccdth_red")
plot = somme_yth_red / np.max(somme_yth_red)
plt.plot(np.arange(80),plot[10:90],'blue')
plt.show()



# detection des indices des pics_FP sur somme_yFP
pics_FP = []
for i in range(1,len(somme_yFP)-1):
    if somme_yFP[i-1]<somme_yFP[i] and somme_yFP[i+1]<somme_yFP[i] :
        pics_FP.append(i)
        
# detection des indices des pics_th sur somme_yth
pics_th = []
for i in range(1,len(somme_yth)-1):
    if somme_yth[i-1]<somme_yth[i] and somme_yth[i+1]<somme_yth[i] :
        pics_th.append(i)
        

# construction du vecteur YccdFP
YccdFP = np.zeros(20*100)
for x in range(100):
    for y in range(20) :
        YccdFP[20*x + y] = ccdFP[y,x]

            
# construction du vecteur Yccdth
Yccdth = np.zeros(20*100)
for x in range(100):
    for y in range(20) :
        Yccdth[20*x + y] = ccdth[y,x]   
        
# construction du vecteur Yccdth_red a reduire
Yccdth_red = np.zeros(20*100)
for x in range(100):
    for y in range(20) :
        Yccdth_red[20*x + y] = ccdth_red[y,x]
        
# construction de la matrice AFP :
AFP = np.zeros((20*100,100))
fen_X = 5
for x in pics_FP:
    AFP[20*(x-fen_X) : 20*(x+fen_X) ,x] = YccdFP[20*(x-fen_X) : 20*(x+fen_X)] - np.min( YccdFP[20*(x-fen_X) : 20*(x+fen_X)] )
# plt.figure(14)
# plt.imshow(AFP,aspect='auto')
# plt.show()

# construction de la matrice Ath :
Ath = np.zeros((20*100,100))
fen_X = 5
for x in pics_th:
    Ath[20*(x-fen_X) : 20*(x+fen_X) ,x] = Yccdth[20*(x-fen_X) : 20*(x+fen_X)] - np.min( Yccdth[20*(x-fen_X) : 20*(x+fen_X)] )
plt.figure(72)
plt.imshow(Ath,aspect='auto')
plt.show()


# interpolation entre les raies de  fabry perots de AFP
pics_FP = [35,64]
epaisseur = 20
largeur = fen_X   
x1 = pics_FP[0]
x2 = pics_FP[1]
for x in range(x1+1,x2):
    alpha = (x-x1)/(x2-x1)
    beg = (x-largeur)*epaisseur 
    end = (x+largeur)*epaisseur
    for y in range(beg,end) :
        yx1 = y - (x-x1)*epaisseur
        yx2 = y + (x2-x)*epaisseur 
        AFP[y,x] = alpha*AFP[yx2,x2] + (1-alpha)*AFP[yx1,x1]
# plt.figure(15)
# plt.imshow(AFP,aspect='auto')
# plt.show()

# interpolation entre les raies de  fabry perots de Ath
epaisseur = 20
largeur = fen_X   
x1 = pics_th[0]
x2 = pics_th[1]
for x in range(x1+1,x2):
    alpha = (x-x1)/(x2-x1)
    beg = (x-largeur)*epaisseur 
    end = (x+largeur)*epaisseur
    for y in range(beg,end) :
        yx1 = y - (x-x1)*epaisseur
        yx2 = y + (x2-x)*epaisseur 
        Ath[y,x] = alpha*Ath[yx2,x2] + (1-alpha)*Ath[yx1,x1]
plt.figure(73)
plt.imshow(Ath,aspect='auto')
plt.show()


# test de reduction
Ath_test = Ath[200:1800,10:90]
Cth_test = Ath_test.T.dot(Ath_test)
Yccdth_red_test = Ath_test.T.dot(Yccdth_red[200:1800])
spectre_th_red = np.linalg.solve(Cth_test,Yccdth_red_test)
plt.figure(74)
plt.plot(spectre_th_red,'red')
plt.show()

clean()


# # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # ## # # # # # # ## # # # #

def ajout_psf_gaussienne(ccd,sigma_psf,cen_psf,larg_psf,hau_psf,ampl_psf) :
    
    #creation psf
    psf = np.zeros((hau_psf,larg_psf))
    ampl = ampl_psf
    largeur = sigma_psf
    cen = int(larg_psf/2)
    # construction du profil 1D de la psf
    x = [i for i in range(larg_psf) ]
    y = [ ampl*np.exp( -(i - cen)**2/(2*largeur**2) ) for i in range(larg_psf)]
    
    # construction de la psf 2D
    for i in range(hau_psf) :
        psf[i] = y
        
    # ajout psf
    for x in range(cen_psf-int(larg_psf/2),cen_psf+int(larg_psf/2)):
        for y in range(5,5+hau_psf):
            ccd[y,x] += psf[y-5,x-cen_psf-int(larg_psf/2)] - np.min(psf)

def ajout_psf_gaussienne_inclinee(ccd,sigma_psf,cen_psf,larg_psf,hau_psf,ampl_psf) :
    #cree une psf gaussienne 2D comme la fonction ci dessus mais en deux parties avec l'une decalee par rapport a l'autre comme si cette psf etait inclinee.
    #creation psf
    psf = np.zeros((hau_psf,larg_psf))
    ampl = ampl_psf
    largeur = sigma_psf
    cen = int(larg_psf/2)
    # construction du profil 1D de la psf
    x = [i for i in range(larg_psf) ]
    y = [ ampl*np.exp( -(i - cen)**2/(2*largeur**2) ) for i in range(larg_psf)]
    
    # construction de la psf 2D
    for i in range(int(hau_psf/2)) :
        psf[i] = y
        
    y_shifted = np.zeros(len(y))
    for i in range(len(y)):
        if i == len(y)-1 :
            y_shifted[i] = y[0]
        else :
            y_shifted[i] = y[i+1]
                
    for i in range(int(hau_psf/2),hau_psf) :
        psf[i] = y_shifted
    
    # ajout psf
    for x in range(cen_psf-int(larg_psf/2),cen_psf+int(larg_psf/2)):
        for y in range(5,5+hau_psf):
            ccd[y,x] += psf[y-5,x-cen_psf-int(larg_psf/2)] - np.min(psf)




# # # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # ## # # # # # # ## # # # #
# tests gaussiens ! 
ajout_psf_gaussienne(ccd1,0.1,cen_psf_cons+1,21,10,1)
plotter2D(ccd1,2)

local_th = local_th / 10000
local_FP = local_FP / 1000


ajout_img_ccd(ccd2,15,10,local_th)
ajout_psf_gaussienne_inclinee(ccd2,0.8,90,21,8,1)
ajout_img_ccd(ccd2,65,10,local_FP)
plotter2D(ccd2,-2)

# construction du vecteur ccd1 qui fera donc 20*100 utile pour construire A
Yccd1 = np.zeros(20*100)

for x in range(100) :
    for y in range(20) :
        Yccd1[20*x + y] = ccd1[y,x]
# construction vecteur Yccd2
Yccd2 = np.zeros(20*100)

for x in range(100) :
    for y in range(20) :
        Yccd2[20*x + y] = ccd2[y,x]

# version naive de la reduction
somme_y = []
for i in range(100) :
    somme_y.append(np.mean([ccd2[j,i] for j in range(20)]))
plt.figure(3)
plt.title("Somme sur Y des intensites")
plt.plot(somme_y)
plt.show()

ligne_10 = []
for i in range(100) :
    ligne_10.append(ccd2[10,i])
plt.figure(3)
plt.title("Somme sur Y des intensites")
plt.plot(ligne_10,'green')
plt.show()

# remplissage de la matrice de convolution pour notre longueur d'onde specifique a la psf fixee
A = np.zeros((20*100,100))

for x in range(100):
    for y in range(20*100):
        A[y,x] = shifter(Yccd1,20*(x-cen_psf_cons))[y]

# plt.figure(4)
# plt.imshow(A,aspect='auto')
# plt.show()


# resolution du systeme
C = A.T.dot(A)
Y = A.T.dot(Yccd2)
# resolution python 
spectre = np.linalg.solve(C,Y)
plt.figure(3)
plt.title("Spectre reduit")
plt.plot(spectre,'red')
plt.show()

# resolution gradient conjugue
# Asparse = sparse.csr_matrix(A)
# Ysparse = sparse.csr_matrix(Yccd2)
# spectre_grad = conjugate_grad_prec_sparse(Asparse,Ysparse)
# plt.figure(6)
# plt.title("Spectre reduit avec la methode du gradient conjugue")
# plt.plot(spectre_grad[1],'black')
# plt.plot(spectre_grad[0],'r--')
# plt.show()


begpic = 85
endpic = 95
lambdas = [ i for i in range(begpic,endpic) ]
pic_ligne_10 = ligne_10[begpic:endpic]
pic_red = spectre[begpic:endpic]
result_ligne_10 = gaussian_fit(lambdas,pic_ligne_10)
result_red = gaussian_fit(lambdas,pic_red)
larg_ligne_10 = result_ligne_10[1]
larg_red = result_red[1]
print(larg_ligne_10,larg_red)

