import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import pyfits

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
    condition = (abs(Qfinal-Qinit) > 1e-7) and cpt < cpt_max
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

        k = -2j * np.pi * iList * v / (2*step*lenI)
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

        # Update of the condition
        condition = (abs(Qfinal-Qinit) > 1e-7) and cpt < cpt_max

        cpt += 1

    if (cpt == cpt_max):             # Case of non convergence
        return (-1, 0, 0)
    return (v, Qfinal, sigma_noise)  # Case of convergence
    
    
def find_shift(fichier):
	
	# We find the data of the ref thar
	hdulist = pyfits.open("../DATA/th_calibered.fits")
	ref_inti = hdulist[1].data.field(1)				# We store the data of the intensity		
	ref_absc = hdulist[1].data.field(0)				# We store the data	of the pixels
	hdulist.close()
	
	# We find the data of the current thar
	hdulist = pyfits.open(fichier)
	inti = hdulist[1].data.field(1)				# We store the data of the intensity		
	absc = hdulist[1].data.field(0)				# We store the data	of the pixels
	hdulist.close()
	
	len_order=4612
	
	all_shift=[]
	for o in range(len(inti)/len_order):
		cur = inti[o*len_order:(o+1)*len_order]
		ref_cur = ref_inti[o*len_order:(o+1)*len_order]
		
		(shift,Q,S)=chelli_shift(ref_cur,cur,1)
		
		wave_shift=shift*(ref_absc[2306]-ref_absc[2305])
		new_absc=[i-wave_shift for i in ref_absc[o*len_order:(o+1)*len_order]]
		
		plt.plot(new_absc,cur,label="le shift")
		plt.plot(ref_absc[o*len_order:(o+1)*len_order],cur,label="original")
		plt.plot(ref_absc[o*len_order:(o+1)*len_order],ref_cur,label="ref")
		plt.legend()
		plt.show()
		
find_shift("../FILES/Narval_20161201_154102_th1.fts")
