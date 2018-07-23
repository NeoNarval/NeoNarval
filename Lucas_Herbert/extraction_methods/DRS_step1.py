#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
from Chelli_Module import *


workdir='./temp'
CCDsize=4612
ref=2320  # Central column



dir='13mar18/'


# Get FLAT done
files=sorted([f for f in os.listdir(dir) if f.endswith('fla.fts')])
#f='/Users/arturo/NeoNARVAL/Espadons/1772472f.fits'

img=0.
for f in files:
	a=pyfits.open(dir+f)
	img=img+a[0].data
	a.close()

flat=img/len(files)


#Writing flat file
a=pyfits.open(dir+files[0])
header=a[0].header
a.close()

hdu = pyfits.PrimaryHDU(flat)
hdu.header=header
nombre2=files[0].replace('fla','f10')
hdu.writeto('./temp/'+nombre2,clobber=True)      

#Finding orders

perfref=flat[ref,:]
#For Espadons use
#interorden=[72, 103.0, 134.0, 164.0, 196.0, 228.0, 258.0, 294.0, 327.0, 360.0, 394.0, 428.0, 470.0, 507.0, 546.0, 585.0, 632.0, 668.0, 712.0, 759.0, 810.0, 857.0, 904.0, 958.0, 1006.0, 1065.0, 1119.0, 1183.0, 1236.0, 1300.0, 1359.0, 1421.0, 1492.0, 1560.0, 1637.0, 1715.0]
#For Narval use
interorden=[50, 76.0, 107.0, 138.0, 169.0, 200.0, 231.0, 263.0, 295.0, 329.0, 363.0, 396.0, 431.0, 468.0, 505.0, 543.0, 583.0, 623.0, 664.0, 708.0, 752.0, 798.0, 844.0, 892.0, 943.0, 993.0, 1046.0, 1101.0, 1157.0, 1215.0, 1274.0, 1337.0, 1401.0, 1465.0, 1535.0, 1602.0, 1676.0, 1751.0, 1834.0,1910.0,2000.]
# 

interpol=map(int,[12.0, 16.0, 16.0, 16.0, 15.0, 15.0, 16.0, 16.0, 16.0, 16.0, 16.0, 17.0, 18.0, 18.0, 
19.0, 20.0, 20.0, 20.0, 21.0, 20.0, 23.0, 21.0, 24.0, 25.0, 25.0, 27.0, 27.0, 28.0, 29.0, 30.0, 
31.0, 31.0, 31.0, 34.0, 33.0, 38.0, 37.0, 38.0,33,39])

# plt.plot(perfref)
# for i,l in enumerate(interpol):
# 	plt.axvline(interorden[i],linestyle='-')
# 	plt.axvline(interorden[i]+l,linestyle='--')
# 	plt.text(interorden[i]+l,0,l)
# plt.show()



f=open('./REF/oder_reference.pkl','w')
pickle.dump({'Profil':perfref,'Orders':interorden,'Beams':interpol,'info':'Reference profile across orders and the position of the interorder and interbeam'},f)
f.close()

f=open('./REF/oder_reference.pkl','r')
a=pickle.load(f) 
f.close()
perfref0=a['Profil']
interorden=a['Orders']

v=Chelly(perfref0,perfref,v0=.1)     ##### Needs a gross first step
interorden=interorden+np.floor(v)



cuantos=len(interorden)

mapa=np.zeros((cuantos,CCDsize))
plt.ion()
for orden in range(0,cuantos-2):
	print("Orden {0}".format(orden))
	ini=int(interorden[orden])
	fin=int(interorden[orden+1])
	ordenref=perfref[ini:fin]
	iniD=int(interorden[orden])
	finD=int(interorden[orden+1])
	iniI=int(interorden[orden])
	finI=int(interorden[orden+1])
	corrD=0
	corrI=0
	refintensity=np.amax(ordenref)
	vIold=0
	vDold=0.
	for delta in range(2000):
		plt.clf()
		vI=0
		vD=0
		ordenI=img[ref-delta,iniI:finI]
		ordenD=img[ref+delta,iniD:finD]
		if np.amax(ordenI) > 0.05*refintensity:
			ordenI=ordenI*refintensity/np.amax(ordenI)
			vI=Chelly(ordenref,ordenI,v0=vIold)
		if np.amax(ordenD) > 0.05*refintensity:
			ordenD=ordenD*refintensity/np.amax(ordenD)
			vD=Chelly(ordenref,ordenD,v0=vDold)
			
		# plt.plot(shift_profile(ordenref,0))
# 		plt.plot(shift_profile(ordenI,-vI))
# 		plt.plot(shift_profile(ordenD,-vD))
# 		plt.draw()
# 		print(vI,vD)
		
		if vI> 1:
			iniI+=1
			finI+=1
			corrI+=1
			vI-=1
		if vD> 1:
			iniD+=1
			finD+=1
			corrD+=1
			vD-=1
		# if np.abs(vI-vIold)>1:
# 			stop
# 		if np.abs(vD-vDold)>1:
# 			stop
		vIold=vI
		vDold=vD
		# print(vI,vD)
		# plt.plot(shift_profile(ordenref,0))
# 		plt.plot(shift_profile(ordenI,-vI))
# 		plt.plot(shift_profile(ordenD,-vD))
# 		plt.draw()
		mapa[orden,ref-delta]=corrI+vI
		mapa[orden,ref+delta]=corrD+vD
	plt.plot(mapa[orden,:])
	x=np.arange(ref-2000,ref+2000)
	c,var=np.polyfit(x,mapa[orden,ref-2000:ref+2000],3,cov=True)
	plt.plot(x,c[0]*x**3+c[1]*x**2+c[2]*x+c[3])
	x=np.arange(CCDsize)
	mapa[orden,:]=c[0]*x**3+c[1]*x**2+c[2]*x+c[3]
	plt.title("Orden {0}".format(orden))
	plt.draw()
	
	fout=open("./REF/Order_Curvature.pkl","w")
	pickle.dump({'map':mapa,'info':"Shift computed from Chelli's algorithm using reference profile at pixel 2320"},fout)
	fout.close()




