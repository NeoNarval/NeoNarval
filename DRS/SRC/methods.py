#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate as interp
import pyfits
from datetime import datetime, timedelta
import settings


def desplaza_perfil(prof, pos):
    """
    Decale un vecteur.

    Permet de décaler un vecteur prof par un flottant position
    """

    pF = np.fft.rfft(prof)
    n = pF.size
    k = -2 * np.pi * np.linspace(0, n + 1, n) * pos / (2. * n)
    phase = np.cos(k) + np.sin(k) * 1j
    pF[0:n] = phase * pF[0:n]
    prof = np.real(np.fft.irfft(pF))
    return prof


def planck(wav, T):
    """
    Calcule l'intensité de Planck normalisée
    wav : longueur d'onde en mm
    T : température du corps noir en K
    """
    wav = wav * 10 ** -3
    h = 6.626e-34
    c = 3.0e+8
    k = 1.38e-23
    a = 2.0 * h * c ** 2
    b = h * c / (wav * k * T)
    intensity = a / ((wav ** 5) * (np.exp(b) - 1.0)) / 2004623760881.600098
    return intensity


def blaze(alpha, beta, m):
    """
    Calcule le blaze
    alpha : angle alpha en radians
    beta : angle beta en radians
    m : ordre (en valeur absolue, premier pour Narval = 22)
    """
    facteur = np.cos(alpha) - np.sin(alpha) * (1 / np.tan((2 * alpha + beta) / 2))
    intensity = np.sinc(m * np.cos(alpha) * facteur) ** 2
    return intensity + 0.01     # +0.01 pour eviter d'avoir un blaze nul à certains endroits


# Calcul de la courbe de sensibilité de la CCD
f = open('../DATA/NARVAL_QE.txt', 'r')
f.readline()
wvl = []
QE = []
for line in f:
    wvl.append(float(line.split()[0]))
    QE.append(float(line.split()[1]))
f.close()
wvl = np.asarray(wvl)
QE = np.asarray(QE)

QEf = np.zeros(len(np.arange(300, 1100)))
for i, w in enumerate(np.arange(300, 1000)):
    idx = (np.abs(wvl - w)).argmin()
    idx0 = np.amax([idx - 4, 0])
    idx1 = np.amin([idx + 4, len(wvl) - 1])
    q = interp.UnivariateSpline(wvl[idx0:idx1], QE[idx0:idx1])
    QEf[i] = q(w)

# On allonge artificiellement la courbe pour les longueurs d'onde extrêmes
pente = (QEf[660] - QEf[699]) / 40
for i in np.arange(700, 800):
    nQEf = QEf[699] - pente * (i - 699)
    if nQEf > 0:
        QEf[i] = nQEf
    else:
        QEf[i] = 0


def CCD_factor(wav):
    """Calcule la sensibilité du CC (wav en mm)"""
    wav *= 1e6
    maximum = 90.0083891879
    wav -= 300
    return QEf[int(round(wav))] / maximum


if settings.type_image == "st0":
    # Chargement du fichier de spectre d'étoile
    try:
        f_spectre = open("../DATA/spectre_kappaceti.spctr")
        long_file = sum(1 for line in f_spectre)
        wvl = []
        sp = []
        matrice_spectre = np.zeros((long_file, 2))
        f_spectre.close()
        f_spectre = open("../DATA/spectre_kappaceti.spctr")
        for i, line in enumerate(f_spectre):
            matrice_spectre[i] = [float(line.split()[0]), float(line.split()[1])]
        f_spectre.close()

    except IOError:
        print "Fichier de spectre introuvable"
        matrice_spectre = np.ones((1100 - 300, 2))

    indices = []
    matrice_spectre = matrice_spectre[np.argsort(matrice_spectre[:, 0])]

    # Suppression des longueurs d'ondes dans les recouvrements
    for i in xrange(len(matrice_spectre[:, 0]) - 1):
        if matrice_spectre[i, 0] == matrice_spectre[i + 1, 0]:
            indices.append(i)
    matrice_spectre = np.delete(matrice_spectre, indices, axis=0)
    # Interpolation du spectre
    for i, val in enumerate(matrice_spectre[:, 1]):
        if val < 0:
            matrice_spectre[i, 1] = 0.001
    interp_spectre = interp.interp1d(matrice_spectre[:, 0], matrice_spectre[:, 1], kind='cubic')


def get_spectre_nearest(wav):
    """Retourne le spectre stellaire au plus proche (wav en mm)"""
    wav *= 1e6
    tableau = matrice_spectre[:, 0]
    indice = np.searchsorted(tableau, wav)
    if indice > 0 and (indice == len(tableau) or np.abs(wav - tableau[indice - 1]) < np.abs(wav - tableau[indice])):
        return matrice_spectre[indice - 1, 1]
    else:
        return matrice_spectre[indice, 1]


def get_spectre(wav):
    """Retourne le spectre stellaire interpolé (wav en mm)

    Normalement il n'y a pas besoin d'interpolations, nearest suffit
    """

    wav *= 1e6
    try:
        valeur = interp_spectre(wav)
        if valeur > 0:
            return interp_spectre(wav)
        else:
            raise ValueError
    except ValueError:
        print "erreur !"
        return get_spectre_nearest(wav * 1e-6)


def fabry_perot(wav):
    """Retourne l'intensité du Fabry-Pérot calculé (wav en mm)"""
    wav *= 1e-3
    R = 0.80
    delta = (2351495 * 1e-9) / wav * 2 * np.pi
    return ((1 - R) ** 2) / (1 - 2 * np.cos(delta) * R + R ** 2)


def ecriture_FITS(matrice, type_algo, type_image):
    imageFITS = pyfits.PrimaryHDU(matrice.T)
    imageFITSlist = pyfits.HDUList([imageFITS])
    header = imageFITSlist[0].header
    temps = datetime.utcnow()
    header["DATE"] = temps.strftime("%Y-%m-%dT%H:%M:%S.%f")
    header["TIMEEND"] = (temps + timedelta(minutes=1)).strftime("%H:%M:%S.%f")
    name = "../Brut/simu" + type_algo + temps.strftime("_%Y%m%d_%H%M%S_") + type_image + ".fts"
    imageFITSlist.writeto(name)
    imageFITSlist.close()
    print "ecriture fits terminee"
