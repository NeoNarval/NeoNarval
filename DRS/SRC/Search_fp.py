import pyfits
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from datetime import datetime
import utils

def find_fp(path):
	"""
					A method to find all the flat in the file and store them in a list
	"""
	docs=os.listdir(path)
	fps=[]
	for names in docs :
		if names.endswith("fp0.fts"):			# We assume the file : "XXXXXXa.fits"
			fps.append(path+"/"+names)			# We store the path of the file
	fps.sort()
	return fps
	
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
	
def fps_unbiased(dossier):
	"""
						A method to substract the 10b file to the thar file
							dossier : the name of the directory in which we search the thars file
	"""
	# The path of the path_file
	path_file=os.path.dirname(os.path.dirname(os.path.abspath("Search_bias.py")))+"/DATA/paths.txt"
	
	# The path of the directory wich contains all the original file
	path_brut=find_path(path_file,"Brut")
	
	# The array which contains all the fabry-perrot file
	fps=find_fp(path_brut)
	
	# We find all the 10b file
	path_files10=find_path(path_file,"FILES")				# The path of the directory which contains all the 10b files
	bias=find_10bias(path_files10)							# An array of all the 10b files
	
	while(len(fps)!=0):
		ref_name=fps[0]
		hdulist=pyfits.open(ref_name)
		image_ref = hdulist[0].data.astype(np.float32)		# We store the data of the current f-p
		hdir_ref=hdulist[0].header							# We store the data of the header
		date_ref=hdir_ref["DATE"]
		hdulist.close()
		
		# We remove the f-p file from the list
		fps.remove(ref_name)
		
		# We prepare the name of the new file
		add_name=ref_name.split("/")
		
		# We convert the date from the reference file in a datetime object
		ref=datetime.strptime(date_ref,'%Y-%m-%dT%H:%M:%S.%f')
		# We put the Day and Hour information for the file we are going to create
		name_fichier=ref.strftime('%Y%m%d_%H%M%S')
		
		# We find the bias file
		image_bias=find_10bias_date(bias,ref)
		
		# And we substract the bias file to the image
		image_combine=image_ref-image_bias
		
		# Then we store the image in the new document
		hdu=pyfits.PrimaryHDU((image_combine))
		hdulist=pyfits.HDUList([hdu])
		hdulist[0].header=hdir_ref												# We keep the header
		path = path_files10 + "/Narval_" + name_fichier + "_fp1.fts"
		hdulist.writeto(path,clobber=True)		# the name of the new file : "Narval_YYYYMMDD_HHMMSS_fap.fts"
		hdulist.close()
		print("FP file {0} written".format(path))

		# MODIFICATION PAR ANTHONY 19/06/2017
		cfg1 = utils.CfgFile("../DATA/find_fp_slits_DATA.txt")
		cfg2 = utils.CfgFile("../DATA/Amatrix_DATA.txt")
		cfg1.modif_param("Fabry Perot fts file", path)
		cfg2.modif_param("Fabry Perot fts file", path)

	
fps_unbiased("Brut")
