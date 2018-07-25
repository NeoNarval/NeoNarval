#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math
from csv import reader
import pyfits
import methods as meth
from multiprocessing import Pool
import gc
import scipy.sparse
import settings
import time
# import matplotlib.pyplot as plt
# import cPickle 
# IMPORT DU DARK

tps1 = time.clock()
type_image = settings.type_image
hdulist = pyfits.open('../DATA/Narval_20161003_155730_10b.fts')  #Narval_20161028_053601_10b.fts
dark = hdulist[0].data.T
offset_lin = np.linspace(400, 50, dark.shape[0])
for i, line in enumerate(dark):
    line = line + offset_lin[i]
    dark[i, :] = line

# IMPORT DE LA PSF
file='../DATA/IMA_45Langst5019_S2.csv'
f=open(file,'r')
a=reader(f,delimiter=';')
PSF=np.zeros((320,320))
i=0
for row in a:
    PSF[i,:]=np.asarray(map(float,row))
    i=i+1
f.close()

# CALCUL DE LA PSF SOUS-ECHANTILLONNEE
rapport = 14 # rapport entre la taille d'un pixel image et un pixel PSF
rapport_PSF = 1 # mouais
taillex_pxlimg = 6
tailley_pxlimg = 14
PSF = PSF[160 - tailley_pxlimg/2 * rapport:160 + taillex_pxlimg/2 * rapport,\
 160 - taillex_pxlimg/2 * rapport:160 + taillex_pxlimg/2 * rapport] #valeurs originelles 70:247, 145:178
angle = 9
PSF2 = np.zeros((tailley_pxlimg * rapport_PSF, taillex_pxlimg * rapport_PSF))
X = np.linspace(0, PSF.shape[1] - rapport/rapport_PSF, taillex_pxlimg * rapport_PSF)
Y = np.linspace(0, PSF.shape[0] - rapport/rapport_PSF, tailley_pxlimg * rapport_PSF)

if angle != 0:
    for i in xrange(len(PSF)):
        PSF[i, :] = meth.desplaza_perfil(PSF[i, :], np.tan(np.radians(angle)) * ((len(PSF)/2) -i))
for j, x in enumerate(X):
    for i, y in enumerate(Y):
        somme = 0
        for plus_i in xrange(rapport/rapport_PSF):
            for plus_j in xrange(rapport / rapport_PSF):
                somme += PSF[int(y) + plus_i, int(x) + plus_j]
        PSF2[i, j] = somme / (rapport**2)

# plt.imshow(PSF2)
# plt.show()
# exit()
# CREATION IMAGE2 ET IMAGE
num_ligne = 4096 ##2098
num_colonnes = 4096 #4612
image = np.zeros([num_ligne, num_colonnes])
shift_i = ((PSF2.shape[0])/rapport_PSF + 4) / 2 # decalage image/image2 en i
shift_j = ((PSF2.shape[1])/rapport_PSF + 4) / 2 # decalage image/image2 en j
image2 = np.zeros((num_ligne + shift_i * 10, num_colonnes + shift_j * 10))


def ajouter_PSF(matrice, x, y, facteur):
    PSF3 = PSF2 * facteur / 4
    def infl(cam, PSF): 
        return (min(cam + 1.0, PSF + 1.0/rapport_PSF) - float(max(cam, PSF))) * rapport_PSF
    centre1 = [0] * 2
    for i in xrange(PSF3.shape[0]):
        centre1[0] = y + (- PSF3.shape[0] / 2.0 + i) / rapport_PSF + shift_i
        for j in xrange(PSF3.shape[1]):
            centre1[1] = x + ( - PSF3.shape[1] /2.0 + j)/rapport_PSF  + shift_j
            i1 = int(centre1[0])
            j1 = int(centre1[1])
            inflx0 = infl(i1, centre1[0])
            inflx1 = infl(i1 + 1, centre1[0])
            infly0 = infl(j1, centre1[1])
            infly1 = infl(j1 + 1, centre1[1])
            matrice[i1, j1] += PSF3[i, j] *  inflx0 * infly0
            matrice[i1 + 1, j1] += PSF3[i, j] * inflx1 * infly0
            matrice[i1, j1 + 1] += PSF3[i, j] * inflx0 * infly1
            matrice[i1 + 1, j1 + 1] += PSF3[i, j] * inflx1 * infly1



# CHARACTERISTIQUES DU SYSTEME OPTIQUE
f1 = 1500 # focale 1 en mm
f2 = 388 # focale avant la camera CCD en mm 388
alpha = math.radians(63.477) # angle d'arrivee sur le reseau
gamma = math.radians(0.6) # inclinaison du reseau selon l'axe horizontal
chi = math.radians(56.34) # angle entre les deux prismes
G = 79 # nombre de traits par mm pour le reseau
A = math.radians(34.5)
pixel = 0.0135
beta_max = ((num_colonnes + shift_j * 4) / 2 * pixel) / f2
betas = np.arange(-beta_max, beta_max, 2 * beta_max / num_colonnes / 1) 

def execu(input):
    m = input[0]
    decal_y = input[1]
    matrice_ordre = np.zeros((num_ligne + shift_i * 10, num_colonnes + shift_j * 10))
    lambdas = zip([((np.sin(alpha) + np.sin(alpha + beta))*np.cos(gamma) )/ (G * m) for beta in betas], betas)
    print m

    for _lambda, beta in lambdas:
        
        # Calcul de l'indice de refraction
        #PBL25Y
        n = (1+1.31960626/(1-1.01863415e-02/(_lambda*10**(3))**2)+1.23752633E-01/ \
            (1-4.83593508E-02/(_lambda*10**(3))**2)+2.10055351E-01/(1-2.73272029E+01/(_lambda*10**(3))**2))**.5
        
        '''formule tiree de "Spectrometrie optique" de Patrick Bouchareine 
         page R 6 310 11  . La formule est avec 8 au denominateur, mais defini pour beta/2 donc avec le carre y a un facteur 4. '''
        n1 = n + ((beta*2)**2*(n**2-1)) / (8*n) 
        # PRISME 1
        i1 = math.radians(28.17) # 28.17
        r1 = np.arcsin(np.sin(i1) / n1)
        r2 = A - r1
        i2 = np.arcsin(n1 * np.sin(r2))

        # PRISME 2
        i1_ = chi - i2
        r1_ = np.arcsin(np.sin(i1_) / n1)
        r2_ = A - r1_
        i2_ = np.arcsin(n1 * np.sin(r2_))

        # PRISME 3
        i1__ = chi - i2_
        r1__ = np.arcsin(np.sin(i1__) / n1)
        r2__ = A - r1__
        i2__ = np.arcsin(n1 * np.sin(r2__))


        x_reel = beta * f2 * (1/pixel) + image.shape[1]/2 # position reelle en x("indice decimal")
        y_reel = (i2__ * f2) * (1/pixel) - 12000 + decal_y # position reelle en y ("indice decimale")
        # result.append([_lambda, x_reel])
        
        if y_reel < (num_ligne + 26):
            _CCD_factor = meth.CCD_factor(_lambda)
            _blaze = meth.blaze(alpha, beta, m)
            _spectre = 1
            if type_image == "fp0":
                _spectre = meth.fabry_perot(_lambda) * meth.planck(_lambda, 3454.35)
            elif type_image == "fla":
                _spectre = meth.planck(_lambda, 3454.35)
            elif type_image == "st0":
                _spectre = meth.get_spectre_nearest(_lambda)
            # sortie.append(_spectre)
            # if _spectre > 0.3 and up == False:
            #     up = True
            #     cpt_pique += 1
            # elif _spectre < 0.3 and up == True:
            #     up = False
            decal_x = 14 * np.tan(np.radians(angle)) # 6.5 avant
            # _spectre_sta = meth.get_spectre_nearest(_lambda)
            ajouter_PSF(matrice_ordre, x_reel + decal_x, y_reel-14,  100* _spectre)#* _CCD_factor * _blaze)
            ajouter_PSF(matrice_ordre, x_reel , y_reel,  100* _spectre)# * _CCD_factor * _blaze) # utiliser la valeur precise plutot que l'entier idem pour y
            ajouter_PSF(matrice_ordre, x_reel - decal_x, y_reel+14,  100* _spectre)# * _CCD_factor * _blaze)
    
    return scipy.sparse.bsr_matrix(matrice_ordre)

def try_execu(m):
    try:
        return execu(m)
    except:
        print "Problem with %d" % m[0]

def main():
    for decal_y in [0]: # pour les boucles de decalage
        image2 = np.zeros((num_ligne + shift_i * 10, num_colonnes + shift_j * 10))
        
        # Multi-processing
        p = Pool(4)
        matrices_ordre = p.map(try_execu, zip(range(22, 62), [decal_y] * len(range(22, 62))))
        p.terminate()
        p.join()
        # matrices_ordre = []
        # for m in [42]:
        #     matrices_ordre.append(execu((m, 0)))

        print "Calculs termines"
        for mat in matrices_ordre:
            image2 += mat.toarray()
        print "Somme terminee"
        # plt.plot(nombre_piques[:, 0], nombre_piques[:, 1])
        # plt.show()
        # Ajout du dark et calcul du bruit
        image = image2[shift_i: shift_i + num_ligne, shift_j:shift_j + num_colonnes] \
             + 300 * np.ones((num_ligne, num_colonnes)) 
        print "maximum = %f, minimum = %f" % (image.max(), image.min())
        del image2
        for i, ligne in enumerate(image):
            for j, case in enumerate(ligne):
               image[i, j] += np.random.normal(loc=0.0, scale=np.sqrt(np.sqrt(image[i, j])))
        
        meth.ecriture_FITS(image, "01", type_image)
        del matrices_ordre
        gc.collect()
        tps = time.clock() - tps1
        print "Termine ! Le temps d'execution vaut : %f s" % tps
if __name__ == '__main__':
    main()
