# -*- coding: utf-8 -*-

"""

"""

import os
import gc
import utils
import pyfits
import tables
import cPickle
import numpy as np
from scipy import sparse, linalg,signal
from scipy.sparse import linalg as slinalg
import matplotlib.pyplot as plt
from csv import reader
import pywt

def rint(x): return np.round(x).astype(int)
def fint(x): return np.floor(x).astype(int)
def cint(x): return np.ceil(x).astype(int)
def sign(x): return int(np.sign(x))

def desplaza_perfil(prof,pos):
	pF=np.fft.rfft(prof)
	n=pF.size
	#nq=n/2
	k=-2*np.pi*np.linspace(0,n+1,n)*pos/(2.*n)
	phase=np.cos(k)+np.sin(k)*1j
	pF[0:n+1]=phase*pF[0:n+1]
	#pF[n-nq:n-1]=np.conj(pF[nq:1:-1])
	prof=np.real(np.fft.irfft(pF,n=len(prof)))
	return prof
def desplaza_x(img,cx):

	tam=np.shape(img)
	x=np.arange(tam[0])
	if cx!=0:
		for i in range(tam[1]):
			img[:,i]=np.interp(x,x+cx,img[:,i])
	
	return img

def fill_A_fp(left_lane, thickness, ref_data, tab_centre, art, ini_index, end_index):
    pix = 1. / art

    (lenX, _) = ref_data.shape
    nb_col = rint(thickness)

    row, col, data = [], [], []
    search_window = 5 * pix

    for [x_pos, _] in tab_centre[0:10]:
        if (ini_index < x_pos < end_index):
            min_x = rint(max(ini_index, x_pos-search_window))
            max_x = rint(min(end_index, x_pos+search_window+1))

            for x in range(min_x, max_x):
                left = left_lane[x]
                for y in range(rint(left), rint(left+thickness)):
                    row.append((x-ini_index)*nb_col + y-rint(left))
                    col.append(rint((x_pos-ini_index) * pix))

                    data.append(float(ref_data[x, y]))

    A_shape = (lenX*nb_col, rint(lenX*pix))
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float')

    return A


def fill_A_matrix(left_lane, thickness, tab_centre, ref_data, ini_index, end_index):
    art = 1
    pix = 1. / art
    rpix = rint(pix)
    (lenX, _) = ref_data.shape
    nb_col = rint(thickness)

    A_fp = fill_A_fp(left_lane, thickness, ref_data, tab_centre, art, ini_index, end_index)
    print(A_fp.todense().shape)	
    print(np.shape(np.squeeze(A_fp.getrow(0).todense()))) 
    A=np.zeros((2000,4612))
    for i in np.arange(0,2000):
	A[i,:]=  np.squeeze(A_fp.getrow(i).todense())

    plt.imshow(A,aspect='auto')
    plt.show()
    
    u, s, v = slinalg.svds(A_fp, k=len(tab_centre) + 10)

    del(A_fp)
    gc.collect()

    C1 = np.matrix(v.transpose()) * np.matrix(linalg.pinv(np.diag(s)))
    U = np.matrix(u.transpose())

    A_fp_inv = np.dot(C1, U)

    row, col, data = [], [], []
    search_window = 5 * pix

    for i in range(0, len(tab_centre) - 1):
        x0, x1 = tab_centre[i][0], tab_centre[i+1][0]

        if (ini_index < x0 < end_index and ini_index < x1 < end_index):
            ax0, ax1 = pix * x0, pix * x1
            arx0, arx1 = rint(ax0), rint(ax1)

            for j in range(arx0+1, arx1):
                g = j
                j = rint(j*art)
                min_x = rint(max(ini_index, g - search_window))
                max_x = rint(min(end_index * rpix, g + search_window+1))

                for x in range(min_x, max_x):
                    for y in range(thickness):

                        row.append(g - rint(pix * ini_index))
                        col.append(rint((x-ini_index)*art)*nb_col + y)

                        pos_0 = rint((x-g)*art + x0 - (1+art) * ini_index)*nb_col + y
                        pos_1 = rint((x-g)*art + x1 - (1+art) * ini_index)*nb_col + y

                        p0 = A_fp_inv[arx0, min(rint(lenX*thickness)-1, rint(pos_0))]
                        p1 = A_fp_inv[arx1, min(rint(lenX*thickness)-1, rint(pos_1))]

                        coeff = float(g-arx0) / (rint(x1-x0) * rpix)
                        p = (p0 + coeff * float((p1-p0))) * art

                        data.append(p)

    A_shape = (rint(lenX*pix), rint(lenX*thickness))
    A = sparse.coo_matrix((data, (row, col)), shape=A_shape, dtype='float').tolil()

    for [x_pos, _] in tab_centre:
        if (ini_index < x_pos < end_index):
            min_x = rint(max(0, x_pos-search_window))
            max_x = rint(min(lenX, x_pos+search_window+1))

            for x in range(min_x, max_x):
                for y in range(thickness):
                    A[rint(x_pos * pix), x*nb_col +
                      y] = A_fp_inv[rint(x_pos * pix), x*nb_col + y]

    A_full = A.tocsr()

    del(A)
    del(A_fp_inv)
    gc.collect()

    return A_full


def cut_border_edge(left_lane):
    lenX = len(left_lane)
    init = 0
    end = lenX - 1
    test = False

    for i in range(lenX):
        if not left_lane[i]:
            if not test and i >= init:
                init = min(i+1, end)
            elif test and i >= init:
                end = max(init, i-1)
                break
        else:
            test = True

    return (init, end)


def find_path(path_file, key):
    """
        A method to find the path of the directory we want in a txt file
            path_file   : the path of the file which contains all the paths for all the directory
            key         : the name of the directory we want the path
        return : 
            value       : the value corresponding to the key
    """

    file = open(path_file, "r")
    text = file.read().split("\n")
    for j in range(len(text)):
        # We assume that in the file : "Key : value"
        line = text[j].split(" : ")
        if line[0] == key:
            return line[1]

    return ("not found")



    """
        Returns an array containing the FP centroids of all considered orders.
        it also gives a colourized of the image_data to have a quick overview
        of the result.
    """
    # path_file = "Amatrix_DATA.txt"

cfgFile = utils.CfgFile("./Amatrix_DATA.txt") #ALA: use my own config file

thickness_data_file = cfgFile.get_param("Lane thickness file")
envelope_data_file = cfgFile.get_param("Lane envelope file")
fp_position_file = cfgFile.get_param("FP slits position file")
fp_file = cfgFile.get_param("Fabry Perot fts file")
order = int(float(cfgFile.get_param("Order n°")))
lane = int(float(cfgFile.get_param("Lane n°")))

    # thickness_data_file = find_path(path_file, "Lane thickness file")
    # envelope_data_file = find_path(path_file, "Lane envelope file")
    # fp_position_file = find_path(path_file, "FP slits position file")
    # fp_file = find_path(path_file, "Fabry Perot fts file")
    # order = int(float(find_path(path_file, "Order n°")))
    # lane = int(float(find_path(path_file, "Lane n°")))

pos_global = cPickle.load(open(fp_position_file, 'rb'))
envelope_data = cPickle.load(open(envelope_data_file, 'rb'))
thickness_data = cPickle.load(open(thickness_data_file, 'rb'))

image_file = pyfits.open(fp_file)
fp_data = image_file[0].data.astype(np.float32)
image_file.close()

Ordering=2*order + lane - 1 #Setup 2 voies
Ordering=order	#Setup 1 voie e.g. Lune	


left = envelope_data[:,Ordering]
tab_pos = pos_global[Ordering]
#thickness = rint(thickness_data[2*order + lane - 1])
thickness = rint(thickness_data[Ordering]) #Pour la Lune

#print(np.shape(tab_pos))	
#plt.plot(tab_pos[0:30],'o')
#plt.show()	
    		
data_file = cfgFile.get_param("Test fts file")
(ini, end) = (200,4400)#cut_border_edge(left)
df=pyfits.open(data_file)
data = df[0].data.astype(np.float32)
df.close()
img=np.zeros(((end-ini+1),thickness))
for i in range(ini,end):
	img[i-ini,:]=fp_data[i,rint(left[i]):rint(left[i]+thickness)]-10


PSFsource='Gaussian'
Zoo=1.
r=3
aZ=1./Zoo
bZ=(1-aZ)*ini


A_matrix=np.zeros(((end-ini+1)*thickness,rint((end-bZ)*Zoo-ini+1)))

PSF=np.zeros((2*r+1,thickness))

wavelet=pywt.Wavelet('db4')


if PSFsource=='Zeemax':
	# IMPORT DE LA PSF
	file='/home/main/Documents/DRS/DATA/IMA_45Langst5019_S2.csv'
	f=open(file,'r')
	a=reader(f,delimiter=';')
	PSF_th=np.zeros((320,320))
	i=0
	for row in a:
	    PSF_th[i,:]=np.asarray(map(float,row))
	    i=i+1
	f.close()

	for i in xrange(len(PSF)):
	        PSF_th[i, :] = desplaza_perfil(PSF_th[i, :], np.tan(np.radians(10.)) * ((len(PSF_th)/2) -i))
elif PSFsource=='ThAr':
	psff=pyfits.open('/home/main/Documents/Arturo/Source_files/Narval_20170609_175308_th1.fts')
	Thp = psff[0].data.astype(np.float32)	
	psff.close()
	i=207+ini
	PSF_th=Thp[i-r:i+r+1,rint(left[i]):rint(left[i]+thickness)]-10
	
elif PSFsource=='Gaussian':
	x=np.arange(-10*r-5,(10*r)+5)
	x0=(np.arange(0,thickness)-thickness/2.)*np.sin(-10*np.pi/180.)	
	
#plt.ion()
for wvl in np.arange(ini,(end-bZ)*Zoo):

	wvlPix=aZ*wvl+bZ 
	residuo=wvlPix-rint(wvlPix)

	cual=np.argmin(abs(wvlPix-np.asarray(tab_pos)[:,0]))
	
	try:
		aqui=rint(np.asarray(tab_pos)[cual,0])-ini
		norma=img[aqui,:]
		norma=norma/norma.max()
	except:
		pass


	if PSFsource=='Gaussian':
		for i in range(thickness):

			perfh=np.exp(-(x-x0[i]-10*residuo)**2/(2*(2.**2)))
			PSF[:,i]=np.sum(perfh.reshape((2*r+1,10)),1)
		PSF=PSF/np.amax(PSF)
		B=np.copy(img)*0.
			
		fi=min(rint(wvlPix-ini)+r+1,len(B))
		inic=max(0,rint(wvlPix-ini)-r)
		#B[inic:fi,:]=PSF.dot(np.diag(norma))[0:fi-inic,:]
		B[inic:fi,:]=PSF[0:fi-inic,:]
	elif PSFsource=='ThAr':
		B=np.copy(img)*0.
		PSF=np.zeros((2*r+1,thickness))
		for i in range(thickness):
			PSF[:,i]=desplaza_perfil(PSF_th[:,i],residuo)
		
		
		PSF=PSF/PSF.max()
		fi=min(rint(wvlPix-ini)+r+1,len(B))
		inic=max(0,rint(wvlPix-ini)-r)
		B[inic:fi,:]=PSF[0:fi-inic,:]

	elif PSFsource=='FP':

		
		B=np.copy(img)*0.
		FP0=np.zeros((2*r+1,thickness))
		FP1=np.zeros((2*r+1,thickness))
		if rint(np.asarray(tab_pos)[cual,0]) > wvlPix:
			#Nearest FP is at larger wavelengths
			aqui=rint(np.asarray(tab_pos)[cual,0])
			FP1=fp_data[aqui-r:aqui+r+1,rint(left[aqui]):rint(left[aqui]+thickness)]-10
			deltaW=aqui-rint(np.asarray(tab_pos)[cual-1,0])
			aqui=rint(np.asarray(tab_pos)[cual-1,0])
			FP0=fp_data[aqui-r:aqui+r+1,rint(left[aqui]):rint(left[aqui]+thickness)]-10
			wv0=np.asarray(tab_pos)[cual-1,0]
		else:
			aqui=rint(np.asarray(tab_pos)[cual,0])
			FP0=fp_data[aqui-r:aqui+r+1,rint(left[aqui]):rint(left[aqui]+thickness)]-10
			deltaW=-aqui+rint(np.asarray(tab_pos)[cual+1,0])
			aqui=rint(np.asarray(tab_pos)[cual+1,0])
			FP1=fp_data[aqui-r:aqui+r+1,rint(left[aqui]):rint(left[aqui]+thickness)]-10
			wv0=np.asarray(tab_pos)[cual,0]


		FP=FP0+(wvlPix-wv0)*(FP1-FP0)/deltaW
		
		#plt.clf()
		#plt.subplot(3,1,1)
		#plt.imshow(FP0)
		#plt.title('aqui={0} wvlPix={1}'.format( rint(np.asarray(tab_pos)[cual,0]),wvlPix))
		#plt.subplot(3,1,2)
		#plt.imshow(FP1)
		#plt.subplot(3,1,3)
		#plt.imshow(FP)
		#plt.show()
		#stop
		mask=np.zeros(2*r+1)
		mask[r]=1.
		for i in range(thickness):
			FP[:,i]=desplaza_perfil(FP[:,i],residuo)*mask
			
		fi=min(rint(wvlPix-ini)+r+1,len(B))
		inic=max(0,rint(wvlPix-ini)-r)
		B[inic:fi,:]=FP[0:fi-inic,:]
		if (rint(wvl-ini)==150):
			stop
	elif PSFsource=='Zeemax': ##Only valid for 2x3 setup
		PSF0=np.zeros((320,22))
		for i in range(320):
			perf=desplaza_perfil(PSF_th[:,i],14*residuo)
			PSF0[i,:]=np.sum(perf[7:315].reshape((22,14)),1)
		PSF1=np.zeros((22,22))
		for i in range(22):
			perf=PSF0[7:315,i].reshape((22,14))
			PSF1[i,:]=np.sum(perf,1)
		PSF=PSF1.T
		#####Not finished
	A_matrix[:,rint(wvl-ini)]=np.hstack(B) 
	#if (150<wvl-ini<250):
	#	print(wvlPix,inic,fi,rint(wvl-ini))


print('{0} Arturo/pixels'.format(1./aZ))
print('{0} km/sec per bin'.format(2.6*aZ))		
#plt.subplot(1,2,2)
#plt.imshow(A_matrix,aspect='auto')
#plt.show()






img=np.zeros((end-ini+1,thickness))
for i in range(ini,end):
	f=1.
	if i-ini <10:
		f=(i-ini)/10.
	if end-i <10:
		f=(end-i)/10.
	
	img[i-ini,:]=(data[i,rint(left[i]):rint(left[i]+thickness)]-10)

C=np.dot(A_matrix.T,A_matrix)
C=C+np.identity(len(C))*0.0001
G=np.linalg.cholesky(C)
y=np.linalg.solve(G,np.dot(A_matrix.T,np.hstack(img)))
x=np.linalg.solve(G.T,y)
#x=signal.wiener(x,5)

#nI=np.reshape(np.dot(A_matrix,x),img.shape)

#sky=signal.medfilt(img-nI,51)
#img=img-sky
#y=np.linalg.solve(G,np.dot(A_matrix.T,np.hstack(img))[4:len(C)-4])
#x=np.linalg.solve(G.T,y)

plt.subplot(2,1,1)
plt.imshow(img.T,aspect='auto')

data_file = cfgFile.get_param("Flat file")
#(ini, end) = (2000,2499)#cut_border_edge(left)
df=pyfits.open(data_file)
flat = df[0].data.astype(np.float32)
df.close()
imgf=np.zeros((end-ini+1,thickness))
for i in range(ini,end):
	imgf[i-ini,:]=flat[i,rint(left[i]):rint(left[i]+thickness)]-10
y=np.linalg.solve(G,np.dot(A_matrix.T,np.hstack(imgf)))
fsp=np.linalg.solve(G.T,y)
fsp=fsp/np.mean(fsp)


#sigma=np.sqrt(np.clip(np.abs(np.dot(A_matrix,x)),1,65000))

#ny=np.hstack(img)/sigma
#A2=A_matrix.copy()
#for i in range(50,451):
#	A2[i,:]=A2[i,:]/sigma[i]
#C=np.dot(A_matrix.T,A2)
#C=C+np.identity(len(C))*0.0001
#G=np.linalg.cholesky(C)
#y=np.linalg.solve(G,np.dot(A_matrix.T,ny))
#x2=np.linalg.solve(G.T,y)

file_Lune='/home/main/Documents/DRS/DATA/lune_calibered.fits'
l=pyfits.open(file_Lune)
a=l[1].data
l.close()
wvl=a['wavelength_lane1'][order*4612:(order+1)*4612-1]

plt.subplot(2,1,2)
l=np.arange(len(C))/Zoo
plt.plot(wvl[ini:end+1],x/x.max())
xf=x/fsp
plt.plot(wvl[ini:end+1],xf/x.max())
#coefs=pywt.wavedec(x[20:580]/x.max(),wavelet)
c=0
d=0
#for i in range(len(coefs)):
#	d=d+np.count_nonzero(coefs[i])
#	cuales=np.where(abs(coefs[i])<500/x.max())
#	coefs[i][cuales[0]]=0.
#	c=c+np.count_nonzero(coefs[i])
#plt.plot(l[20:580],pywt.waverec(coefs,wavelet))
#plt.xlim((0,500))
plt.ylim((0,5))
plt.show()

stop
#A_matrix = fill_A_matrix(left, thickness, tab_pos, fp_data, ini, end)

#name_A_matrix = './Amatrix_' + 'OR' + str(order) + '_LA' + str(lane) + '.npz'
#sparse.save_npz(name_A_matrix, A_matrix)

#cfgFile.modif_param('A matrix', name_A_matrix)


