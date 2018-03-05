#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from datetime import datetime
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


def find_flats(path):
	"""
					A method to find all the flat in the file and store them in a list
						path : the path of the directory where we are searching the flat
	"""
	docs=os.listdir(path)
	flats=[]
	for names in docs :
		if names.endswith("fla.fts"):			# We assume that the original flats file are "XXXXXXf.fits"
			flats.append(path+"/"+names)		# We store the entire path in one array
	flats.sort()
	return flats

def find_10bias(path):
	
	"""
					A method to find all the 10bias in the file and store them in a list
						path : the path of the directory where we are searching the 10b files
	"""
	docs=os.listdir(path)
	bias=[]
	for names in docs :
		if names.endswith("10b.fts"):			# We assume that the name file are : "Narval_YYYYMMDD_HHMMSS_10b.fts"
			bias.append(path+"/"+names)			# We store the entire path in one array
	bias.sort()
	return bias


def find_10bias_date(bias,ref_time):
	
	"""
					A method to compare the date of a file with all the 10b.fts file in order to find the matching bias file
						bias : an array with all the 10b.fts files
						ref_time : the time of the reference file
					return :
						image_bias : the image data of the matching 10b files				
	"""
	
	# parameters
	condition=0										# condition == 1 when we find the matching 10b file
	j=0												# To be sure we stay in the array 
	
	while(condition==0 and j<len(bias)):
		hdulist=pyfits.open(bias[j])
		hdir_current=hdulist[0].header
		date_current=hdir_current["DATE"]			# We find the date of the current 10b file
		hdulist.close()
		
		# We convert the date from the reference file in a datetime object
		current_time=datetime.strptime(date_current,'%Y-%m-%dT%H:%M:%S.%f')
		
		# We need the total time between the two file
		diff = abs(current_time-ref_time)
		
		if diff.total_seconds()<14440 :				# If the total time is less than 4 hours, we consider they are from the same night
			hdulist=pyfits.open(bias[j])
			image_bias = hdulist[0].data.astype(np.float32)			# So we store the image of the matching 10b file
			hdulist.close()
			condition=1								# And the condition to end the loop
		else :
			j=j+1

	
	if j==len(bias):
		print ("No 10b")
	else :
		return (image_bias)
				
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

def combine_flats(dossier):
	"""
					A method which combine all the flats from the same night
						dossier : the directory in which the routine search the files
					For each night, it create a 10f.fts file which is the mean of all the flats from the same night. 
					However, for each flat, we substract the corresponding 10b.fts file.
	"""
	# First, we store all the flat and the black files
	path_file=os.path.dirname(os.path.dirname(os.path.abspath("Search_bias.py")))+"/DATA/paths.txt"
	path_brut=find_path(path_file,dossier)
	flats=find_flats(path_brut)
	# print(flats)
	path_files10=find_path(path_file,"FILES")
	bias=find_10bias(path_files10)
	# print(bias)
	# Our list of parametres
	count=0						# count the number of flats that are use for one file

	dump=[]	
	
	while (len(flats)!=0):
		use_flat=[]
		# We open the first flat which is our reference
		ref_name=flats[0]
		hdulist=pyfits.open(flats[0])
		image_ref = hdulist[0].data.astype(np.float32)					 # We store the image of the flat
		(lenX,lenY) = image_ref.shape								     # We find the dimension of the image
		hdir_ref=hdulist[0].header
		date_ref=hdir_ref["DATE"]										 # We will use the date of the reference many times
		end_current=hdir_ref["TIMEEND"]
		hdulist.close()
		
		time_end=end_current
		np.asarray(image_ref)
		# We remove the reference from the list
		flats.remove(ref_name)
		
		add_name=ref_name.split("/")
		use_flat.append(add_name[len(add_name)-1])
		# We prepare the final image
		image_combine=image_ref

		
		# We convert the date from the reference file in a datetime object
		ref=datetime.strptime(date_ref,'%Y-%m-%dT%H:%M:%S.%f')
		# We put the Day and Hour information for the file we are going to create
		name_fichier=ref.strftime('%Y%m%d_%H%M%S')
		count+=1
		
		# We find the bias file
		image_bias=find_10bias_date(bias,ref)
		
		# And we substract the bias file to the image
		image_combine-=image_bias
		
		j=0
		while (j<len(flats)):
                        # print (len(flats))
			# We find the flat from the same night as the reference flat
			current_name=flats[j]
			hdulist=pyfits.open(flats[j])
			image= hdulist[0].data.astype(np.float32)				# We store the image of the flat
			hdir1=hdulist[0].header
			date_current=hdir1["DATE"]								# And the date
			end_current=hdir1["TIMEEND"]							# And the end time of exposure
			hdulist.close()
			
			np.asarray(image)
			# We convert the date of the current flat
			current=datetime.strptime(date_current,'%Y-%m-%dT%H:%M:%S.%f')
			
			duree =abs(current-ref)                                # We need the total time between the two file
			

	
			if duree.total_seconds() < 14400 :						# If the time between the two files is less than 4 hours, we consider they are from the same night
							
				(shift,Q,s)=chelli_shift(image_ref[lenX/2,:],image[lenX/2,:],1,0.05)
				print("Decalage flat n {} chelli : {} ".format(count, shift))
				# We check if the flat is similar to the reference image, thanks to Chelli's algorithm
				if shift<0.5:
					# If the shift is small enough, we add the current image and substract the bias image
					image_combine=image_combine+(image-image_bias)
					count+=1										# We increase the number of flat that are used to create the new image
					flats.remove(current_name)						# We delete the current flat from the original array
					add_name=current_name.split("/")
					use_flat.append(add_name[len(add_name)-1])		# We add the name of the used flat in the list
					# print ("count =",count)
					time_end=end_current
				else :
					# Or we store the flat in the dump list
					dump.append(flats[j])
					flats.remove(current_name)
			else :
				print("nouvelle nuit")
				j=j+1
		
		# Then we store the image in the new document
		hdir_ref["TIMEEND"]=time_end							# We change the timeend in the header by the timeend of the last flat we used
		hdu=pyfits.PrimaryHDU((image_combine/count))			# We store the mean image of all the flat we used in the new file
		hdulist=pyfits.HDUList([hdu])
		hdulist[0].header=hdir_ref								# We use the reference header to create the new one
		hdulist[0].header.add_history(str(use_flat))
                path = path_files10 + "/Narval_" + name_fichier + "_10f.fts"
                hdulist.writeto(path, clobber=True)		# We create the name of our new file : "Narval_YYYYMMDD_HHMMSS_10f.fts"
		hdulist.close()
                cfg_Amatrix = utils.CfgFile("../DATA/Amatrix_DATA.txt")
                cfg_Amatrix.modif_param("Flat fts file", path)
		count=0
		print("{0} written".format(path))
		
	return dump

dump=combine_flats("Brut")	
