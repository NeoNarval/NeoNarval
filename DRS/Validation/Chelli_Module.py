#! /usr/bin/env python

import numpy as np


def Chelly(fr,f,v0=0.1):
	IFou=np.fft.rfft(fr)*np.conjugate(np.fft.rfft(f))

		#################################
		# inicio del Marquardt-Levenberg
		#

		

	NW=len(IFou)
	

	v=v0	
	delta=1.
	Bucle=True
	iter=0
	
	k=-2.*np.pi*np.arange(NW)*v/(2.*NW)

	CFou=IFou*(np.cos(k)+1j*np.sin(k))
	
	sigma=0.5*(np.sum(fr)*np.sum(f)+(np.sum(fr)+np.sum(f))*np.abs(CFou))
	
	Chi2=np.sum(np.imag(CFou)**2/sigma**2)
	
	while(Bucle):
	#v given in pixels
		iter+=1
	#Exponent of phase of the shift 
		
		Num=np.sum(np.arange(NW)*np.real(CFou)*np.imag(CFou)/sigma**2)
		Den=np.sum(np.arange(NW)**2*np.real(CFou)**2/sigma**2)


		deltav=delta*(NW/(2*np.pi))*Num/Den
		
		k=-2.*np.pi*np.arange(NW)*(v+deltav)/(2.*NW)

		CFou=IFou*(np.cos(k)+1j*np.sin(k))
		sigma=0.5*(np.sum(fr)*np.sum(f)+(np.sum(fr)+np.sum(f))*np.abs(CFou))
		
		NewChi2=np.sum(np.imag(CFou)**2/sigma**2)
		
		
		if (NewChi2 > Chi2):
			delta=delta/10.
			#stop
		else:	
			#print('New Chi2:{0}'.format(v+deltav))
			#if ((Chi2-NewChi2))<1e-20:
			if (abs(deltav))<1e-6:
				Bucle=False
				v=v+deltav	
				#print('Converged {0} {1}'.format(deltav,NewChi2))
			#print('New shift {0} pixels Chi2 {1}'.format(v+deltav,NewChi2))
			delta=1.
			Chi2=NewChi2
			v=v+deltav
		if iter > 500:
			print("Too Long")
			print(deltav)
			#stop
			Bucle=False
	return v
	
def xshift(kfft,phi):
	tam=np.shape(kfft)[0]
	#Nyquist frequency.
	nq=tam  
	#This assumes that np.rfft has been used, so that only one half of the array is returned
	
	#Exponent of phase of the shift except for sqrt(-1).
	k=-2.*np.pi*np.arange(nq)*phi/(2*tam)
                
	#Set complex phase.
	phase=np.cos(k)+1j*np.sin(k) 
    
	#Apply phase shift from 0 to Nyquist frequency.
	kfft1=phase*kfft[0:nq]
	
	return kfft1

def shift_profile(perf,dx):
	kfft0 = np.fft.rfft(perf)
	kfft1 = xshift(kfft0,dx)
	xxxx1 = np.fft.irfft(kfft1)
	return xxxx1
	
