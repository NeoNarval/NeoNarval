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
import scipy.optimize
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
from pylab import *
from mpl_toolkits.mplot3d import Axes3D




def gaussian2D(grid,cenX,cenY,amp,widX,widY):
    """
    This function returns a 2D gaussian using the given parameters. 
    Inputs :
    - grid = [X,Y] where X and Y are 1D vectors
    - cenX and cenY are floats
    - amp is a float
    - widX and widY are floats
    """
    x = grid[0]
    y = grid[1]
    return(np.matrix([ [amp*np.exp(-(x[i]-cenX)**2/(2*widX**2))*np.exp(-(y[j]-cenY)**2/(2*widY**2)) for j in range(len(y))] for i in range(len(x))]))
    


def gfit2D(data):
    
    data = np.matrix(data)
    data = data - np.min(data)
    data = data / float(np.sum(data))
    
    X = np.arange(np.shape(data)[0]) 
    Y = np.arange(np.shape(data)[1])
    grid = [X,Y]

    fig = figure(100)
    ax = Axes3D(fig)
    X, Y = np.meshgrid(X, Y)
    Z = data
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')
    show()

    sumraw = np.zeros_like(data[0])
    for i in range(np.shape(data)[0]):
        sumraw += data[i]
    naive_cenY = np.sum([(k+1)*sumraw[0,k] for k in range(np.shape(data)[1])]) -1
    
    sumcol = np.zeros_like(data[:,0])
    for i in range(np.shape(data)[1]):
        sumcol += data[:,i]
    naive_cenX = np.sum([(k+1)*sumcol[k,0] for k in range(np.shape(data)[0])]) -1
        
    naive_ampl = np.max(data)
    naive_widX = np.std(sumcol)
    naive_widY = np.std(sumraw)
    
    try :
        gaussian2D_model = lmfit.Model(gaussian2D)
        params = gaussian2D_model.make_params(cenX=naive_cenX,cenY=naive_cenY,amp=naive_ampl,widX=naive_widX,widY=naive_widY)
        result = gaussian2D_model.fit(data,params,grid=grid)
        best_gaussian2D_fit = result.best_fit
        best_gaussian2D_fit = np.matrix(best_gaussian2D_fit)
        fig = figure(101)
        ax = Axes3D(fig)
        X = np.arange(np.shape(data)[0])
        Y = np.arange(np.shape(data)[1])
        X, Y = np.meshgrid(X, Y)
        Z = best_gaussian2D_fit
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')
        show()        
        return(best_gaussian2D_fit)
    
    except :
        
        print("Fit failed")
        fig = figure(102)
        ax = Axes3D(fig)
        X = np.arange(np.shape(data)[0])
        Y = np.arange(np.shape(data)[1])
        X, Y = np.meshgrid(X, Y)
        Z = gaussian2D(grid,naive_cenX,naive_cenY,naive_ampl,naive_widX,naive_widY)
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='hot')
        show()
        return(gaussian2D(grid,naive_cenX,naive_cenY,naive_ampl,naive_widX,naive_widY))
        
    
M = np.matrix([[0,0,1,1,1,1,0,0,0],[0,0,0,2,3,2,1,0,0],[0,1,1,2,3,4,1,0,0],[0,1,1,2,5,2,2,1,0],[0,1,2,2,4,3,1,0,0],[0,1,2,3,4,2,1,0,0],[0,0,0,1,1,2,1,0,0],[0,0,1,2,3,1,1,1,0],[0,0,1,1,2,0,0,0,0]])
gfit2D(M)