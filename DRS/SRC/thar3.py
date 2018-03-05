import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy.ma as ma
import cPickle as pickle
from scipy import interpolate
import pyfits
import sys
import utils

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
    
	
def find_wave(fichier):
	print "ok"
	# On lit l'atlas solaire sur a
	atlas="../DATA/SolarSpectrum.pkl"
	f=open(atlas,'r')
	a=pickle.load(f)
	f.close()
	
	# We open the fits files 
	hdulist = pyfits.open(fichier)
	inti = hdulist[1].data.field(1)				# We store the data of the intensity		
	absc = hdulist[1].data.field(0)				# We store the data	of the pixels
	hdulist.close()

	# The len of one order
	len_order=4612

	# Here are all the value to help chelli for each order. The first two zeros are for the order 0 and 1, which I was unabled to caliber
	found_v=[0,0,-5,-2,-2,-2,-0.8,-2,0,-0.5,-1,0.5,0.5,0.8,0.8,0.8,0.8,0.8,0.8,2,4,4,5,5,6,7,7,7.5,7.5,8,8,12,12,15,15,16]
	for o in range(len(inti)/len_order):
		# On lit le spectre de l'ordre considere
		b=inti[o*len_order:(o+1)*len_order]
		
		# Je vais faire l'hypothese que le blaze est pil dans le pixel centrale de la CCD
		# alors la il y p0 pixels a droite et autant a gauche
		p0=int(len(b)/2)
		
		#Les characteristiques du spectro.  Je change le blaze de 63.4 a 63.495
		alpha_Blaze=63.495*np.pi/180.
		gamma=0.6*np.pi/180.
		G=79. # grooves/mm
		F=388.  #mm focal length
		p=13.5e-3 # 12 micron pixel in mm
		order=21+o
		#La longueur d'onde de blaze est...
		lamda0=1e7*2*np.sin(alpha_Blaze)*np.cos(gamma)/(G*order) #in Angstroms
		#mais je ne m'en sert pas
		
		#Si le blaze tombe sur le pixel centrale, ou beta=alpha_Blaze, a gauche je sors avec l'angle
		beta=alpha_Blaze-(p0*p)/F
		# qui correspond a la longueur d'onde
		l0=1e7*(np.sin(alpha_Blaze)+np.sin(beta))*np.cos(gamma)/(G*order) #in Angstroms
		#Je calcule la meme chose pour la droite de la CCD
		beta=alpha_Blaze+(p0*p)/F
		l1=1e7*(np.sin(alpha_Blaze)+np.sin(beta))*np.cos(gamma)/(G*order) #in Angstroms
		# Ce calcul est redondant car....voir les lignes qui suivent.
		
		
		# Je cree aussi un array w2 que, pixel a pixel, me donne la longueur d'onde theorique
		pixel=np.arange(-p0,p0,1)
		beta=alpha_Blaze+(pixel*p)/F
		w2=1e7*(np.sin(alpha_Blaze)+np.sin(beta))*np.cos(gamma)/(G*order) #in Angstroms
		# l0 et l1 auraient pu etre definis comme les extremes de w2.....
	
		# Je vais chercher dans l'atlas le morceaux de spectre choisi
		donde=ma.masked_inside(a["Wavelength"],l0,l1)
		w=np.asarray(a["Wavelength"])[donde.mask]
		atlas=np.asarray(a["Intensity"])[donde.mask]
	
		# We interpolate the atlas in order to create a new atlas with the matching length
		p=interpolate.InterpolatedUnivariateSpline(w,atlas)
		glob_ref=[p(i) for i in w2]
		
		# We prepare the coordinates for the polynomial fit
		pix=[]
		wave=[]
		#plt.plot(b)
		#plt.plot(glob_ref)
		#plt.show()
		
		# The value to help the Chelli routine
		help_v=found_v[o]
		for i in range(100,len(b)-100,10):
	
			# We take a window of 200 pixels centered on the pixel we are currently studying
			left=i-100
			right=i+100
			
			# We create the two images for Chelli
			current=b[left:right]
			ref=np.asarray(glob_ref[left:right])
			
			# We use chelli to find the local shift(in pixel) between the atlas an our image
			(shift,Q,S)=chelli_shift(ref/max(ref),current/max(current),1,help_v)
	
			# For the 8th order only
			if o==8:
				if abs(help_v-shift)>=10:
					shift=help_v 
			
			# And we add the coordinate
			pix.append(i)
			wave.append(w2[i]-(shift)*(w2[i]-w2[i-1]))		# We add the shift in Angstroms to the first value of the wavelength
			
			# We change the value to help Chelli
			help_v=shift
	
		# We interpolate 
		z=np.polyfit(pix,wave,4)
		p2=np.poly1d(z)
		print (l0,l1)
		
		# And create our new array which give us four ach pixel the corresponding wavelength
		result=[p2(i) for i in range(len(b))]
		
		for i in range(len(result)):
			absc[o*len_order+i]=result[i]
	
	#plt.plot(absc,inti)
	#plt.show()
	
	return (absc,inti)

# The file of the normalized lunar spectrum			
fichier = "../FILES/table2.fits"

# We store the new wavelength and the intensity of the lunar spectrum	
(absc,inti)=find_wave(fichier)	

# We create the new fits file
col1 = pyfits.Column(name='wavelength_lane1', format='E', array=absc)
col2 = pyfits.Column(name='intensity_lane1', format='E', array=inti)
# col3 = pyfits.Column(name='wavelength_lane2', format='E', array=tableau[2])
# col4 = pyfits.Column(name='intensity_lane2', format='E', array=tableau[3])

cols = pyfits.ColDefs([col1, col2]) # , col3, col4])
tbhdu2 = pyfits.BinTableHDU.from_columns(cols)

tbhdu1 = pyfits.PrimaryHDU(header=pyfits.open(fichier)[0].header)
#tbhdu1.header["FP_FILE"] = cfg.get_param("Fabry Perot fts file").split("/").pop()
#tbhdu1.header["STAR_FILE"] = cfg.get_param("Test fts file").split("/").pop()

final_HDU = pyfits.HDUList([tbhdu1, tbhdu2])

final_HDU.writeto('../DATA/luna_calibered.fits', clobber=True) # changer le nom du fichier pour mettre le nouveau

print "Reussi !"
