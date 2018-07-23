#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python's modules imports
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit
import math as maths
import numpy.ma as ma
from scipy import interpolate

# Definition of global values
global order_len # length of an order
order_len = 4612
global order_total # number of orders
order_total = 36

# Opening some files to get data :
f = open("ThAr_calibered_original_intensitites.pkl",'r')
Sr = pickle.load(f)
f.close()

f = open("ThAr_calibered_lambdas.pkl",'r')
Lr = pickle.load(f)
f.close()

HDUlist = pyfits.open("spectrum_Narval_20161201_154102_th1.fits")
data = HDUlist[1].data
Sth1 = data['intensity_lane1'][0:4612*36]

HDUlist = pyfits.open("Narval_20161003_160900_th2.fits")
data = HDUlist[1].data
S1 = data['intensity_lane1'][0:4612*36]

HDUlist = pyfits.open("Narval_20161004_053934_th2.fits")
data = HDUlist[1].data
S2 = data['intensity_lane1'][0:4612*36]

"""
The following function is used to test chelli algorithm. It takes a rerefrence spectrum and shift it from the given delta.
Inputs :
- Sr : reference spectrum
- delta : shift to compute
Outputs :
- S : shifted spectrum
"""

def shifter(Sr, shift):
    
    fSr = np.fft.rfft(Sr)
    iList = np.arange(len(fSr))
    k = -2j*np.pi*iList*1.0/order_len/order_total*shift
    fS = np.exp(k)*fSr
    S = np.fft.irfft(fS)
    
    # plt.plot(Sr,color='black')
    # plt.plot(S,color='red')
    # plt.show()
    # 
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
# # Then, let's grab some more accuracy (1/100 of a pixel)
    # jlist = [i*1.0/100.0 for i in range(-11,12)]
    # for j in jlist :
    #     Srj = shifter(Sr,shift+j)
    #     # Computing the current least square indicator for the current k shift    
    #     Qj = np.sum( [ (Srj[i]-S[i])**2 for i in range(len(Srj)) ] )
    #     # Updating only if the lq indicator is better ad recording the best found shift until here...
    #     if  Qj < Q :
    #         Q = Qj
    #         shift = shift+j 
    #         
        # Update of the condition
        condition = (abs(Qfinal-Qinit) > 1e-20) and cpt < cpt_max

        cpt += 1

    if (cpt == cpt_max):             # Case of non convergence
        return (-1, 0, 0)
        print(" /!\ Chelli has not converged! ")
    return (v, Qfinal, sigma_noise)  # Case of convergence


""" 
Since chelli_shift seems to compute better results on small arays and relatively bad ones on big arays, we are going to write a function splitting the input into several parts, computing the shift on each part and then computing the average shift using the shifts computed no each part. We will use it in our algorithms. The Inputs and Outputs we be the same as chelli_shift.
"""
def chelli_shift_v2(Sr,S,step,shift0):
    
    # We will suppose that the lengths of Sr and S are approximatively the same. We will only take a look at the length of Sr to determine if we need (or not) to split the spectrum and compute chelli_shift on smaller parts of it or directly on the global spectrum.
    if len(Sr) < 200 :
        return(chelli_shift(Sr,S,1,shift0))
        
    else :
        
        N = len(Sr)/200        
        shifts = []
        
        for k in range(N) :
            
            Srk = Sr[200*k:200*(k+1)]
            Sk = S[200*k:200*(k+1)]
            result = chelli_shift(Srk,Sk,1,shift0)
            shifts.append(result[0])
            
        # We have to compute the end of the spectrum too...
        Srk = Sr[len(Sr)-200 : len(Sr)]
        Sk = S[len(Sr)-200 : len(Sr)]
        result = chelli_shift(Srk,Sk,1,0)
        shifts.append(result[0])
        
        shift = np.mean(shifts)
        return(shift)


"""
The following function will find the shift between two spectrum but only with an error of 1 pixel : it is a way to find a big shift. Since chelli shift is only working on less than one pixel shifts we need to first find the shift with a precision of one pixel, then apply chelli shift. This function is made for this objective.
Input :
- Sr : reference spectrum.
- S : shifted spectrum.
Output :
- shift : the shift in pixel with a precision of one pixel.
"""

def find_big_shift(Sr,S):
    
    shift = 0
    
    # Initiating the least square indicator with an initial value
    Q = np.sum( [ (Sr[i]-S[i])**2 for i in range(len(Sr)) ] )
    
    # We are going to improve the accuracy step by step, dividing it by 10 in each loop... We start between -110 and 110, then -11 to 11, the -1.1 to 1.1 so we don't miss anything.
    for precision in range(6):
        
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
The following function will do the job of shift_finder, but applied to our case with spectrum divided into several orders.
Inputs :
- Sr : the reference spectrum.
- S : the shifted spectrum (which we want to identify accurately).
Ouputs :
- shift_results : the shift found with high accuracy in indices, in Angstroms and its equivalent in m/s
"""
S = shifter(Sr,1e-6)
    
def find_shift(Sr,S,Lr):
    
    plt.plot(Sr,color='black')
    plt.plot(S,color='red')
    plt.show()
    
    # Creation of the list which will contain the results of the computation.
    shift_results = []
    
    # First of all, since chelli shift is only working on less than one pixel shifts we need to first find the shift with a precision of one pixel, then apply chelli shift. Actually, find_big_shift is more accurate than 1 pixel but it is bonus ofr chelli.
    
    global_approx_shift = find_big_shift(Sr,S)
    print(" Approximate shift for the global spectrum : "+str(global_approx_shift)+" pixels")
    
    
    # The Chelli will be applied to each order individually using the global shift found previously, for more precision. Note that Chelli seem to work better on a small spectrum so we also ave interest to split the work and compute each order individually.
    
    for order in range(len(Sr)/order_len):
        
        
        print("========== Order "+str(order)+" ==========")
        # Selection of the order's intensities for Sr, S and of the order's wavelengths
        order_Sr = Sr[order*order_len:(order+1)*order_len]
        order_S = S[order*order_len:(order+1)*order_len]
        order_lambdas = Lr[order*order_len:(order+1)*order_len]
        
        # Intializing order_shift
        order_shift = global_approx_shift
        
        # Then we shift Sr using the order_shift to help Chelli by giving Chelli a spectrum which has a small shift to compute.
        shifted_Sr = shifter(Sr,order_shift)
        order_shifted_Sr = shifted_Sr[order*order_len:(order+1)*order_len]
        
        # Initalizing the shift for Chelli (with an arbitrary value)
        local_shift = 0.0
        
        # Then we let Chelli do his job, each iteration using the result of the previous one to improve it.
        for i in range(5):
            result = chelli_shift_v2(order_shifted_Sr,order_S,1,local_shift)
            local_shift = result
        
        # Once the local accurate shift is found, let's just add it to the global one and the job is done!
        order_shift = order_shift + local_shift
        print("Shift : " + str(order_shift) + " pixels")
        
        # Finally we compute the radial velocity which goes with the wavelengths shift and we add the result to the output list...
        order_lambdas_shift = order_shift*( Lr[order*(order_len)+int(order_len/2) + 1] - Lr[order*(order_len)+int(order_len/2)] )
        order_rvel = order_lambdas_shift*3*10**8/Lr[order*(order_len)+int(order_len/2)]
        
        print("Radial velocity corresponding to the detected shift in wavelengths : "+str(order_rvel)+" m/s")
        shift_results.append( [ order,order_shift,order_lambdas_shift,order_rvel ] )
    
    return(shift_results)
    

"""
We can also fidnd a big shift with the help of a cross correlation function. This function will implement the finding of a shift thanks to the cross correlation.
Inputs :
- a : array or list
- b : array or list (a with a shift which we want to comute)
Output :
- ind : the shift between a and b
"""

def correlator(a,b):
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
    
    
"""
Same function as find_shift but using the correlator function to compute the approximate shift between Sr and S.
- Sr : the reference spectrum.
- S : the shifted spectrum (which we want to identify accurately).
Ouputs :
- shift_results : the shift found with high accuracy in indices, in Angstroms and its equivalent in m/s
"""

def find_shift2(Sr,S,Lr) :

    plt.plot(Sr,color='black')
    plt.plot(S,color='red')
    plt.show()
    
    # Creation of the list which will contain the results of the computation.
    shift_results = []
    
    # First of all, since chelli shift is only working on less than one pixel shifts we need to first find the shift with a precision of one pixel, then apply chelli shift. Actually, find_big_shift is more accurate than 1 pixel but it is bonus ofr chelli.
    
    global_approx_shift = correlator(Sr,S)
    print(" Approximate shift for the global spectrum : "+str(global_approx_shift)+" pixels")
    
    
    # The Chelli will be applied to each order individually using the global shift found previously, for more precision. Note that Chelli seem to work better on a small spectrum so we also ave interest to split the work and compute each order individually.
    
    for order in range(len(Sr)/order_len):
        
        
        print("========== Order "+str(order)+" ==========")
        # Selection of the order's intensities for Sr, S and of the order's wavelengths
        order_Sr = Sr[order*order_len:(order+1)*order_len]
        order_S = S[order*order_len:(order+1)*order_len]
        order_lambdas = Lr[order*order_len:(order+1)*order_len]
        
        # Intializing order_shift
        order_shift = global_approx_shift
        
        # Then we shift Sr using the order_shift to help Chelli by giving Chelli a spectrum which has a small shift to compute.
        shifted_Sr = shifter(Sr,order_shift)
        order_shifted_Sr = shifted_Sr[order*order_len:(order+1)*order_len]
        
        # Initalizing the shift for Chelli (with an arbitrary value)
        local_shift = 0.0
        
        # Then we let Chelli do his job, each iteration using the result of the previous one to improve it.
        for i in range(5):
            result = chelli_shift_v2(order_shifted_Sr,order_S,1,local_shift)
            local_shift = result
        
        # Once the local accurate shift is found, let's just add it to the global one and the job is done!
        order_shift = order_shift + local_shift
        print("Shift : " + str(order_shift) + " pixels")
        
        # Finally we compute the radial velocity which goes with the wavelengths shift and we add the result to the output list...
        order_lambdas_shift = order_shift*( Lr[order*(order_len)+int(order_len/2) + 1] - Lr[order*(order_len)+int(order_len/2)] )
        order_rvel = order_lambdas_shift*3*10**8/Lr[order*(order_len)+int(order_len/2)]
        
        print("Radial velocity corresponding to the detected shift in wavelengths : "+str(order_rvel)+" m/s")
        shift_results.append( [ order,order_shift,order_lambdas_shift,order_rvel ] )
    
    return(shift_results)
    