#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from math import *
import cPickle as pickle
import os
import sys
import utils

# ================================================ Chelli ========================================================================
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

	if (cpt == cpt_max):			 # Case of non convergence
		return (-1, 0, 0)
	return (v, Qfinal, sigma_noise)  # Case of convergence

# ================================================================================================================================   

# ================================================ Main ==========================================================================
"""
	This section contains the main function of the algorithm
"""

def Ref_localisation(fichier,nb_voie=2):
	"""
				Function to find the reference function for each order. This reference function will be use to locate each order
					fichier : the flat_field document
					nb_voie : number of lane for each order ( 2 or 3)
				Return an array of the left and right limit for the ref function for each order
	
	"""
	
	# We open the data of the flat field and store them in image_data
	hdulist = pyfits.open(fichier)
	image_data = hdulist[0].data
	(lenX,lenY) = image_data.shape
	hdulist.close()
	
	# We add an off set to avoid negative values
	image_data=image_data+500
	
	
	# Our list of parameter
	if (nb_voie == 2 or 3) :
		search_window=35			# a window in which we will search the minimum
		cross_correl_limit = 700		# the limit in pixel after which we use the cross correlation to find the order for the central row
	else :
		search_window=60
		cross_correl_limit = 1250		# the limit in pixel after which we use the cross correlation to find the order for the central row
	nb_order_right=35				# the number of order we want after the brighter order
	nb_order_left=6  				# the number of order we want before the brighter order
	central_row_index=lenX/2		# the index of the central row
	
	# now we search our reference function in the central row for the brighter order
	central_row=image_data[central_row_index,:]													# the data in the central row
	max_loc=np.argmax(central_row)																# the pixel of the global max for the central_row
	left_min_central=np.argmin(image_data[central_row_index,max_loc-search_window : max_loc])	# the left border of the ref function. We search at the left of the global maximum within the search_window
	right_min_central=np.argmin(image_data[central_row_index,max_loc:max_loc+search_window])	# the right border of the ref function. We search at the right of the global maximum within the serach_window
	
	# Now we have our ref function
	ref_function_central=image_data[central_row_index,max_loc-search_window+left_min_central:max_loc+right_min_central]
	
	# We correlate this function with the entire row
	correl=np.correlate(central_row,ref_function_central,mode="same")
	
	# we look at the right of the brighter order
	# up to now, we work with the data from the camera.

	envelope_right=np.zeros((1,nb_order_right))			# we create the envelope in which the ref function will be stored
	envelope_right=envelope_right.astype(int)			# we make the as int to avoid any problem because they will be used as indexes in image_data
	left_min=max_loc-search_window+left_min_central		# we create left_min and right_min which are up to now the left and right of the ref_function_central
	right_min=max_loc+right_min_central
	search_window_right=search_window					# and a new search_window
	j=0
	# now we can search for the other ref function	 
	# for this part we work on the camera data

	while (right_min<cross_correl_limit and j<nb_order_right-1):
		left_min=right_min																							# we assume that the right limit of the function n-1 is the left limit of the function n
		right_min=left_min+5+np.argmin(image_data[central_row_index,left_min+5:left_min+5+search_window_right])		# we look at the next minimum. the '+5' is a safety margin so we are sure the next minimum will be at the right of the order
		envelope_right[0,j]=left_min																				# we store left and right in the envelope
		envelope_right[0,j+1]=right_min																				
		search_window_right=search_window_right+1																	# at the right of the brighter order, the space between tow order increase. Therefore we have to increase our serach_window
		j=j+1
	#up to now we use the cross_correlation
	
	# But before we create a boolean to see if we have reach the limit of the order requested right after the brighter one
	order_limit = False
	if j == nb_order_right :
		order_limit = True
	j=j-1
	while (j<len(envelope_right[0,:])-1 and order_limit==False):
		left_min=envelope_right[0,j]																							# once again, we assume the right limit of the function n-1 is the left limit of the function n
		max_loc_correl=left_min+3+np.argmax(correl[left_min+3:left_min+3+search_window_right])									# we look at the maximum of the correl function within our search_window.There again, the '+3' is a safety margin
		right_min=max_loc_correl+np.argmin(image_data[central_row_index,max_loc_correl:max_loc_correl+search_window_right])		# the max of the correl function give us the position of the order. So at the right of this position we search for the next minimum
		envelope_right[0,j+1]=right_min
		j=j+1
	
	#now we look at the left of the brighter order
	#for Narval, we can't use the cross_correl. The space between two order is to tight and so the paeks in the cross correlation gives us the position of the order but also of the inter_order.
	#Therefore is complicated to identified which on is which
	# the routine is the same but we go backwards now, starting at the brighter order once again
	envelope_left=np.zeros((1,nb_order_left))
	envelope_left=envelope_left.astype(int)
	right_min=max_loc+right_min_central
	left_min=max_loc-search_window+left_min_central
	j=1

	while j<=nb_order_left-1:
		right_min=left_min
		#print central_row_index
		#print str(right_min - 1 - search_window)
		left_min=right_min-1-search_window+np.argmin([image_data[central_row_index,right_min-1-search_window:right_min-1]])
		envelope_left[0,nb_order_left-j]=right_min
		envelope_left[0,nb_order_left-j-1]=left_min
		search_window=search_window-1
		j=j+1

	# to conclude, we make the two envelope as one
	envelope=np.zeros((1,nb_order_left+nb_order_right))
	for j in range (nb_order_left):
		envelope[0,j]=envelope_left[0,j]
	for j in range (nb_order_right):
		envelope[0,j+nb_order_left]=envelope_right[0,j]

	envelope=envelope.astype(int)
	print('envelope',envelope)
	plt.figure()
	plt.plot(image_data[central_row_index,:])
	for i in envelope[0]:
		plt.axvline(i)
	plt.show()
	
	return envelope   
	
def Order_loacalisation(fichier,nb_voie=2):
	"""
					Function that gives the envelope of all the orders  bright enough on the CCD.
						fichier : the file of the flat field
						nb_voie : the number of lanes for each order ( 2 or 3 )
					return a matrix of all the border for all the order, and an array with the width for each order
	
	"""
	
	# We open the data of the flat field and store them in image_data
	hdulist = pyfits.open(fichier)
	image_data = hdulist[0].data
	(lenX,lenY) = image_data.shape
	hdulist.close()
	
	# We add an off set to avoid negative values
	image_data=image_data+500
	
	
	# our list of parameters							
	central_row_index=lenX/2							# the index of the central row
	step=1												# the step for chelli function
	epaisseur=[]										# the width of each order
	a=0													# all those paramaters are for the last orders
	polys=[]
	b=0
	c=0
	d=0
	old_a=0
	old_b=0
	old_c=0
	old_d=0
	old_a2=0
	old_b2=0
	old_c2=0
	old_d2=0

	# we use the previous function to get the ref function for all the order
	ref_envelope=Ref_localisation(fichier,nb_voie)
	#print(ref_envelope)
	envelope_global=np.zeros((lenX,len(ref_envelope[0,:])))		# the matrix in whih the order limits will be stored
	pas=1														# the step for chelli function
		
	for j in range (len(ref_envelope[0,:])-1):

		#print("OK")
		# we create our ref function for the current order and we normalize it
		ref_function_central=image_data[central_row_index,ref_envelope[0,j]:ref_envelope[0,j+1]]/np.max(image_data[central_row_index,ref_envelope[0,j]:ref_envelope[0,j+1]])
		
		# We add the length of the order where it belongs
		ref_length=ref_envelope[0,j+1]-ref_envelope[0,j]
		epaisseur.append(ref_length)
		
		left_central=ref_envelope[0,j]			
		right_central=ref_envelope[0,j+1]
		old_v=0									# we create two velocity both given by chelli function.
		new_v=0
		envelope=np.zeros((lenX,2))				# the matrix in which will be stored the data for th current order 
		down_left=left_central					# we use the left and right limit of the ref function for the beginning
		down_right=right_central
		#print('ref done')
		
		# now we find the border of the order for the down part of the ccd
		border_down=lenX
		for i in range(central_row_index,lenX):
			# we find the function of the current row. We know the shift is little so we assume that the border ( left and right ) between the row n-1 and n is quite the same.
			# Therefore,for the row n, we look between the border of n-1
			
			# we find the current function and we normalize i
			current_function = image_data[i,down_left:down_right]/np.max(image_data[i,down_left:down_right])
			
			# if the length between our current array and the reference one isn't the same, we have reach the end of the CCD. Therefore our border down is the current one
			if len(current_function)!=len(ref_function_central):
				border_down=i
				break
			
			#print(ref_function_central,current_function)
			(new_v,Q,S)=chelli_shift(ref_function_central,current_function,step,0.20)							# chelli gives us the exact shift
	
			envelope[i,0]=down_left+new_v					# with chelli's velocity, we now have the border for this order for the current row
			envelope[i,1]=down_right+new_v
			down_left=down_left+int(old_v+new_v)			# now, we add to our search indexes the value of the old_v and the new_v. We make them as int because they will be used as indexes
			down_right=down_right+int(old_v+new_v)
			old_v=new_v
			
		#print('down done')
		
		# now we do the exact same for the up part of the ccd
		up_left=left_central
		up_right=right_central
		border_up=0
		for i in range(central_row_index,-1,-1):

			current_function = image_data[i,up_left:up_right]/np.max(image_data[i,up_left:up_right])
			
			if len(current_function)!=len(ref_function_central):
				border_up=i
				break
				
			(new_v,Q,S)=chelli_shift(ref_function_central,current_function,step,0.20)
			
			envelope[i,0]=up_left+new_v
			envelope[i,1]=up_right+new_v
			
			up_left=up_left+int(old_v+new_v)
			up_right=up_right+int(old_v+new_v)
			old_v=new_v

		#print('up done')
		if j<len(ref_envelope[0,:])-6:	
			
			# and now we fit the true envelope
			x=np.arange(border_up,border_down)

			y_left=envelope[border_up:border_down,0]
			z_left=np.polyfit(x,y_left,4)
			p_left=np.poly1d(z_left)
			
			y_right=envelope[border_up:border_down,1]
			z_right=np.polyfit(x,y_right,4)
			p_right=np.poly1d(z_right)
			
			a=p_left[4]
			b=p_left[3]
			c=p_left[2]
			d=p_left[1]		
			absc=np.arange(lenX)
			
			polys.append(p_left)
			# and we store the fitted envelope
			for i in range(border_up,border_down):
				envelope[i,0]=p_left(absc[i])
				envelope[i,1]=p_right(absc[i])
		
		else :
			# for the last order, due to the noise, the chelli routine fails a little to find the envelope on the edge of the camera.
			# So we use the previous envelope to find the new one
			#print('fin')
			a=old_a+(old_a-old_a2)
			b=old_b+(old_b-old_b2)
			c=old_c+(old_c-old_c2)
			d=old_d+(old_d-old_d2)
			absc=np.arange(lenX)
			polys.append(np.poly1d([a,b,c,d,ref_envelope[0,j]]))
			
			for i in range(border_up,border_down):
				envelope[i,0]=a*(i**4-central_row_index**4)+b*(i**3-central_row_index**3)+c*(i**2-central_row_index**2)+d*(i-central_row_index)+ref_envelope[0,j]
				envelope[i,1]=a*(i**4-central_row_index**4)+b*(i**3-central_row_index**3)+c*(i**2-central_row_index**2)+d*(i-central_row_index)+ref_envelope[0,j+1]
				if envelope[i,0]>=lenY:
					envelope[i,0]=0
				if envelope[i,1]>=lenY:
					envelope[i,1]=0
					envelope[i,0]=0

		envelope_global[:,j]=envelope[:,0]			# and now we store the border of the order in the global matrix
		envelope_global[:,j+1]=envelope[:,1]
		

		old_a2=old_a
		old_b2=old_b
		old_c2=old_c
		old_d2=old_d
		
		old_a=a
		old_b=b
		old_c=c
		old_d=d
		print('Order {0} of 40'.format(j+1))
	plt.figure()
	plt.imshow(image_data.T, aspect ='auto')
	plt.figure()
	for j in range (len(ref_envelope[0,:])-1):
		plt.plot(envelope_global[:,j])
	plt.show()
	return (envelope_global,epaisseur,polys)

def lanes2(fichier,envelope,epaisseur_order):
	"""
						Function that give for each order the limit of the 2 lanes within
							fichier :  the data
							envelope :  the envelope of the order previously find
							epaisseur_order :  the width for the each order previsously find
						Return a matrix of all the border for each lanes and the width for each lanes
	"""
	
	# We open the data of the flat field and store them in image_data
	hdulist = pyfits.open(fichier)
	image_data = hdulist[0].data
	(lenX,lenY) = image_data.shape
	hdulist.close()
	
	image_data=image_data+500				# We add an off-set
	
	ref_envelope=Ref_localisation(fichier)	# We will need the precise pixels of the reference
	
	# We initialize the matrix to store all the borders
	matrix=np.zeros((lenX,2*len(envelope[1,:])))
	
	epaisseur=[]							# Our array to store the width of all the orders
	
	for j in range(len(envelope[1,:])-1):
		
		# We store the previous border in the new matrix
		matrix[:,2*j]=envelope[:,j]
		matrix[:,2*(j+1)]=envelope[:,j+1]
		
		# We identify the reference function for the current order
		ref_function_central=image_data[lenX/2,ref_envelope[0,j]:ref_envelope[0,j+1]]
		
		# We derive the max and min value in the current order
		max_value=max(image_data[lenX/2,ref_envelope[0,j]:ref_envelope[0,j+1]])
		min_value=min(image_data[lenX/2,ref_envelope[0,j]:ref_envelope[0,j+1]])
		
		# We initialize a safety value  which will be usefull to precisely find where the order is within the previous border
		
		safety = (max_value-min_value)/3
		k=0
		l=0
		#print j
		
		# We find where the order is precisely
		while (image_data[lenX/2,ref_envelope[0,j]+k]<=min_value+safety and k<=ref_envelope[0,j+1]-ref_envelope[0,j]):
			k=k+1
		

		while image_data[lenX/2,ref_envelope[0,j+1]-l]<=min_value+safety and ref_envelope[0,j+1]-l>=ref_envelope[0,j]+k:
			l=l+1
		
		# We store the left and right index previously found
		left_index=ref_envelope[0,j]+k
		right_index=ref_envelope[0,j+1]-l
		print('Order {0}, Left: {1} //  Right: {2} // Thickness: {3}'.format(j+1,left_index,right_index, right_index-left_index))

		# We found the minimum value in the new window which is the the inter_lane
		min_loc_central=left_index+np.argmin(image_data[lenX/2,left_index:right_index])
					
		# We add the width of the 2 lanes in the array
		epaisseur.append(min_loc_central-ref_envelope[0,j])
		epaisseur.append(ref_envelope[0,j+1]-min_loc_central)
		
		# And we find and store the minimum for each lanes for each order, starting from the central lane
		for i in range(lenX/2,lenX):
			if envelope[i,j]==0:
				matrix[i,2*j+1]=0
			else :
				shift=envelope[i,j]-envelope[lenX/2,j]
				matrix[i,2*j+1]=min_loc_central+shift
				if matrix[i,2*j+1]>=lenY:
					matrix[i,2*j+1]=0
		for i in range(lenX/2,-1,-1):
			if envelope[i,j]==0:
				matrix[i,2*j+1]=0
			else :
				shift=envelope[i,j]-envelope[lenX/2,j]
				matrix[i,2*j+1]=min_loc_central+shift
				if matrix[i,2*j+1]>=lenY:
					matrix[i,2*j+1]=0
	
	return (matrix,epaisseur)

def lanes3(fichier,envelope,epaisseur_order):
	"""
						Function that give for each order the limit of the 2 lanes within
							fichier :  the data
							envelope :  the envelope of the order previously find
							epaisseur_order :  the width for the each order previsously find
						Return a matrix of all the border for each lanes and the width for each lanes
	"""
	
	# We open the data of the flat field and store them in image_data
	hdulist = pyfits.open(fichier)
	image_data = hdulist[0].data
	(lenX,lenY) = image_data.shape
	hdulist.close()
	
	# We add an off set to avoid negative values
	image_data=image_data+500
	
	# We initialize the matrix to store all the borders
	matrix=np.zeros((lenX,3*len(envelope[1,:])))
	
	epaisseur=[]
	
	for j in range(len(envelope[1,:])-1):
		
		# We store the previous border in the new matrix
		matrix[:,3*j]=envelope[:,j]
		matrix[:,3*(j+1)]=envelope[:,j+1]
		
		# We identify the reference function for the current order
		ref_function_central=image_data[lenX/2,envelope[lenX/2,j]:envelope[lenX/2,j+1]]
		
		# We derive the max and min value in the current order
		max_value=max(image_data[lenX/2,envelope[lenX/2,j]:envelope[lenX/2,j+1]])
		min_value=min(image_data[lenX/2,envelope[lenX/2,j]:envelope[lenX/2,j+1]])
		
		# We initialize a safety value  which will be usefull to precisely find where the order is within the previous border
		
		safety = (max_value-min_value)/3
		k=0
		l=0
		#print j
		
		# We find where the order is precisely
		while (image_data[lenX/2,envelope[lenX/2,j]+k]<=min_value+safety and k<=envelope[lenX/2,j+1]-envelope[lenX/2,j]):
			k=k+1
		

		while image_data[lenX/2,envelope[lenX/2,j+1]-l]<=min_value+safety and envelope[lenX/2,j+1]-l>=envelope[lenX/2,j]+k:
			l=l+1
		
		# We store the left and right index previously found
		left_index=envelope[lenX/2,j]+k
		right_index=envelope[lenX/2,j+1]-l
		#print(j,left_index,right_index)
		print('Order {0}, Left: {1} //  Right: {2}'.format(j+1,left_index,right_index))
		# We found the minimum value in the new window which is the the inter_lane
		current_min_loc=left_index+np.argmin(image_data[lenX/2,left_index:right_index])
		
		if abs(current_min_loc-left_index)<abs(current_min_loc-right_index):
			left_min_loc=current_min_loc
			matrix[lenX/2,3*j+1]=left_min_loc
			right_min_loc=left_min_loc+3+np.argmin(image_data[lenX/2,left_min_loc+3:right_index])
		else :
			right_min_loc=current_min_loc
			matrix[lenX/2,3*j+2]=right_min_loc
			left_min_loc=left_index+np.argmin(image_data[lenX/2,left_index:right_min_loc-3])
			
		epaisseur.append(left_min_loc-envelope[lenX/2,j])
		epaisseur.append(right_min_loc-left_min_loc)
		epaisseur.append(envelope[lenX/2,j+1]-right_min_loc)
		
		# And we find and store the minimum for each lanes for each order, starting from the central lane
		for i in range(lenX/2,lenX):
			if envelope[i,j]==0:
				matrix[i,2*j+1]=0
			else :
				shift=envelope[i,j]-envelope[lenX/2,j]
				matrix[i,3*j+1]=left_min_loc+shift
				matrix[i,3*j+2]=right_min_loc+shift
		for i in range(lenX/2,-1,-1):
			if envelope[i,j]==0:
				matrix[i,2*j+1]=0
			else :
				shift=envelope[i,j]-envelope[lenX/2,j]
				matrix[i,3*j+1]=left_min_loc+shift
				matrix[i,3*j+2]=right_min_loc+shift
				
	return (matrix,epaisseur)

# ================================================================================================================================

			

# ============================================ Run Subroutine ====================================================================
"""
	This section contains subroutine for the run algorithm
"""	

def Run_envelope(CCD_file,envelope_file,epaisseur_file,polys_file, matrix_file, epaisseur_lane_file,nb_voie=2):
	"""
						Routine to calculate the envelope and the width of each order.
							CCD_file : The fits with the data of the flat
							envelope_file : the name of the pickle file to create, which will contains the matrix of all the envelope
							epaisseur_file : the name of the pickle file to create, which will contains the array of the width of each order
						The envelope is stored in the pickle file, created by the routine
	"""

	(envelope,epaisseur,polys)=Order_loacalisation(CCD_file,nb_voie)

	pickle.dump(envelope,open(envelope_file,"wb"))
	pickle.dump(epaisseur,open(epaisseur_file,"wb"))
	pickle.dump(polys,open(polys_file,"wb"))
	if nb_voie == 1:
		pickle.dump(envelope, open(matrix_file, 'wb'))
		pickle.dump(epaisseur, open(epaisseur_lane_file, 'wb'))

def Open_envelope(epaisseur_file,envelope_file):
	"""
						Routine to open a pickle file to get the envelope and the width of each order
							epaisseur_file : the pickle file to open, which contains the array of the width for each order
							envelope_file : the pickle file to ope, which contains the envelope of each order
						return : 
							the envelope (which is a matrix)
							the width (which is an array)
	"""
	pkl_file=open(envelope_file,"rb")
	envelope=pickle.load(pkl_file)
	
	pkl_file2=open(epaisseur_file,"rb")
	epaisseur=pickle.load(pkl_file2)
	
	return(envelope,epaisseur)
	
def Run_lane(CCD_file,matrix_file,epaisseur_lane_file,epaisseur_file,envelope_file,nb_voie=2):
	"""
						Routine to calculate the envelope and the width for each lane
							CCD_file : The fits with the data of the flat
							matrix_file : the name of the pickle file to create, which will contains the matrix of all the envelope for each lane
							epaisseur_lane_file : the name of the pickle file to create, which will contains the array of the width of each lane
							epaisseur_file : the name of the pickle file which contains the width of all the orders
							envelope_file : the name of the pickle file which contains the envelope of all the orders
						The matrix is stored in a pickle file created by the routine
	"""
	(envelope,epaisseur)=Open_envelope(epaisseur_file,envelope_file)
	if nb_voie==2:
		(matrix,epaisseur_lane)=lanes2(CCD_file,envelope,epaisseur)
	if nb_voie==3:
		(matrix,epaisseur_lane)=lanes3(CCD_file,envelope,epaisseur)
	
	pickle.dump(matrix,open(matrix_file,"wb"))
	pickle.dump(epaisseur_lane,open(epaisseur_lane_file,"wb"))
	
def Open_lane(epaisseur_lane_file,matrix_file):
	"""
						Routine to open a pickle file to get the envelope and the width of each lane
							epaisseur_lane_file : the pickle file to open, which contains the array of the width for each lane
							matrix_file : the pickle file to open, which contains the envelope of each lane
						return : 
							the envelope (which is a matrix)
							the width (which is an array)
	"""
	pkl_file=open(matrix_file,"rb")
	matrix=pickle.load(pkl_file)
	
	pkl_file2=open(epaisseur_lane_file,"rb")
	epaisseur_lane=pickle.load(pkl_file2)
	
	return(matrix,epaisseur_lane)

def find_10flats(dossier):
	"""
						Routine to find the 10f.fts files in a directory and store them in a list
							dossier : the path of the directory in which we search the 10f files
	"""
	docs=os.listdir(dossier)
	flats=[]
	for names in docs :
		if names.endswith("10f.fts"):
			print("10flats found")
			flats.append(dossier+"/"+names)
	flats.sort()
	return flats

def find_path(path_file,param):
	"""
					A method to find the path of the directory we want in a txt file
						path_file : the path of the file which contains all the paths for all the directory
						param : the name of the directory we want the path
					return : 
						path_name : the path of the directory
	"""
	
	my_file=open(path_file,"r")
	within=my_file.read()
	line=within.split("\n")						
	for j in range(len(line)):
		path_name=line[j].split(" : ")			# We assume that in the file : "Name_path : path"
		if path_name[0]==param:
			return path_name[1]
	
	return ("not found")
	
# ================================================================================================================================

# ====================================================== Run =====================================================================
def Run(nb_voie=2):
	j=0
	os.chdir(r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18")
	path_file=r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18"
	directory_path=r"C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\13mar18\ccd_order_data"	
	
	flats10=find_10flats(path_file)

				
	print(len(flats10))

	
	while (len(flats10)!=0):
		j=j+1
		
		fichier=flats10[0]

		envelope_file=directory_path+"/order_position_"+str(j)+".p"			
		epaisseur_file=directory_path+"/order_thickness_"+str(j)+".p"
		polys_file=directory_path+"/polys_envelope_"+str(j)+".p"
		matrix_file=directory_path+"/lanes_position_"+str(j)+".p"
		epaisseur_lane_file=directory_path+"/lanes_thickness_"+str(j)+".p"
		CCD_file=fichier
		
		Run_envelope(CCD_file,envelope_file,epaisseur_file,polys_file, matrix_file, epaisseur_lane_file, nb_voie)
		if nb_voie > 1: # Au cas où il n'y a qu'une voie, on ne calcule que les envelopes des ordres
			Run_lane(CCD_file,matrix_file,epaisseur_lane_file,epaisseur_file,envelope_file,nb_voie)
					
		
		flats10.remove(fichier)
	print("Geometry info from envelope written in {0}".format(directory_path))


Run(nb_voie=2)


	
