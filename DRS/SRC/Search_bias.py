import pyfits
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from datetime import datetime


def find_bias(path):
	"""
					A method to find all the bias in the file and store them in a list
						path : the path of the directory where we are searching the bias file
	"""
	
	docs=os.listdir(path)
	bias=[]
	for names in docs :
		if names.endswith("bia.fts"):				# We assume that the file are : "XXXXXXb.fits"
			bias.append(path+"/"+names)				# We store the entire path in the array
	bias.sort()
	return bias

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

def combine_bias(dossier):
	"""
					A method to create the 10b files.
						dossier : the name of the directory in which we have to find the original flats
	"""
	
	path_file=os.path.dirname(os.path.dirname(os.path.abspath("Search_bias.py")))+"/DATA/paths.txt"		# The path of the path_file
	path_brut=find_path(path_file,dossier)					# The path of the directory which contains all the original file
	bias=find_bias(path_brut)								# An array which contains all the bias file
	count=0													# count the number of flats that are use for one file
		

	
	while (len(bias)!=0):
		use_bias=[]									# An array to store all the name of the bias file which we use to create the current 10b file
		# We open the first bias which is our reference
		ref_name=bias[0]
		hdulist=pyfits.open(bias[0])
		image_ref = hdulist[0].data					# We store the data
		(lenX,lenY) = image_ref.shape				# We find the dimension
		hdir_ref=hdulist[0].header
		date_ref=hdir_ref["DATE"]					# We store the date of the header
		end_current=hdir_ref["TIMEEND"]
		hdulist.close()
		
		time_end=end_current
		np.asarray(image_ref)
		# We remove the reference frome the list
		bias.remove(ref_name)
		
		# We store the name of the file in the use_bias array
		add_name=ref_name.split("/")
		use_bias.append(add_name[len(add_name)-1])
		
		# We prepare the final image
		image_combine=image_ref
		image_null=np.zeros((lenX,lenY))
		
		# We convert the date from the reference file in a datetime object
		ref=datetime.strptime(date_ref,'%Y-%m-%dT%H:%M:%S.%f')
		# We put the Day and Hour information for the file we are going to create
		name_fichier=ref.strftime('%Y%m%d_%H%M%S')
		
		# The number of bias file we use
		count+=1
	
		j=0
		
		while (j<len(bias)):
			
			# We find the bias from the same night as the reference bias
			current_name=bias[j]
			hdulist=pyfits.open(bias[j])
			image= hdulist[0].data.astype(float)				# We store the data
			hdir1=hdulist[0].header
			date_current=hdir1["DATE"]							# We store the date of the current flat file
			end_current=hdir1["TIMEEND"]
			hdulist.close()
			
			np.asarray(image)
			# We convert the date of the current bais
			current=datetime.strptime(date_current,'%Y-%m-%dT%H:%M:%S.%f')
			
			duree =abs(current-ref)								# We need the total time between the two files

			# we check if we are not looking into the reference once again
	
			if duree.total_seconds() < 3600 :					# If the time between two files is less than 1 hour, we consider they are from the same night

				image_null=image
				image_combine=image_combine+image_null
				count+=1										# We increase the number of bias that we use to create the current 10b file
				flats.remove(current_name)						# We delete the current bias from the orignial array
				time_end=end_current
				add_name=bias[j].split("/")						# We add the name of the current bias in the use_bias array
				use_bias.append(add_name[len(add_name)-1])

			else :
				
				j=j+1
		
		# Then we store the image in the new document
		hdir_ref["TIMEEND"]=time_end							# We change the timeend in the header by the timeend of the last bias we used
		hdu=pyfits.PrimaryHDU(image_combine/count)				# We store the mean image of all the bias we used in the nwe file
		hdulist=pyfits.HDUList([hdu])
		hdulist[0].header=hdir_ref								# We use the reference header to create the new one 
		hdulist[0].header.add_history(str(use_bias))
		end_path=find_path(path_file,"FILES")
		hdulist.writeto(end_path+"/Narval_"+name_fichier+"_10b.fts",clobber=True)		# We create the name of our new file : "Narval_YYYYMMDD_HHMMSS_10b.fts"
		print("Writing {0}".format(end_path+"/Narval_"+name_fichier+"_10b.fts"))
		hdulist.close()
		count=0
		
combine_bias("Brut")
	

