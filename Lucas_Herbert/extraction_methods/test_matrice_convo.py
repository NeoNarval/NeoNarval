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

# on fait de la place
clean()



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
ajout_img_ccd(ccdth_red,35,10,local_th)
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
    somme_yth_red.append(np.mean([ccdth_red[j,i] for j in range(20)]))
plt.figure(711)
plt.title("Somme sur Y des intensites ccdth_red")
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


##
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

ajout_psf_gaussienne(ccd1,5,cen_psf_cons,21,8,1)
plotter2D(ccd1,2)


ajout_psf_gaussienne(ccd2,5,60,21,8,1)
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
plt.plot(np.arange(100),somme_y)
plt.show()



# remplissage de la matrice de convolution pour notre longueur d'onde specifique a la psf fixee
A = np.zeros((20*100,100))

for x in range(100):
    for y in range(20*100):
        A[y,x] = shifter(Yccd1,20*(x-cen_psf_cons))[y]

plt.figure(4)
plt.imshow(A,aspect='auto')
plt.show()

# resolution dy systeme

C = A.T.dot(A)
Y = A.T.dot(Yccd2)
spectre = np.linalg.solve(C,Y)
plt.figure(5)
plt.title("Spectre reduit")
plt.plot(spectre)
plt.show()

