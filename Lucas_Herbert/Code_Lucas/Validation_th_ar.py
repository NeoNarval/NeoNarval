"""
Documentation and blabla
IN : fichier calibered.fits (fichier dépouillé avec une série de valeurs de lambdas et leurs intensités, toutes les voies sont mêlées)

OUT : un tableau avec une liste longueurs d'ondes de raies à comparer à un atlas
"""

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit 
from scipy.optimize import leastsq

##
# Opens a cPickle file and read it



def open_pkl(path):
    
    file = open(path,'rb')
    data = cPickle.load(file)
    file.close()
    return(data)
    
# dealing with the .pkl file of the slits

path2 = '/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/ThAr_Line_Cal.pkl'
lambdas2 , intensities2 = open_pkl(path2)['Lines'], open_pkl(path2)['FWHM']




# Open the .fits file which we want to reduce in a table

file_path = '/home/stagiaire/Documents/Données utiles/th_calibered.fits'

data_file = pyfits.open(file_path)

# We search for our data in the .fits file : the wawelenghts and intensities of the thorium argon spectrum, it is in the first extension of the .fits file

data = data_file[1].data

# Now we have the right table, let's builds the two lists of values which interest us

lambdas = data['wavelength_lane1']
intensities = data['intensity_lane1']




# We can plot the spectrum :
#plt.plot(lambdas, intensities)
#plt.show()


# Now let's work on this spectrum : let's try to reduce it in a table of slits to compare with an Athlas

# Since some  wavelengths appear several times in the spectrum (1 time for each order where the wavelength belongs), we need to make the difference between the different orders to have only one slit per wavelength


##______________________________________________________________
#Function which returns a list of wavelenghts for each order (without the jumps in the wavelenghts in the full spectrum due to the different orders)
##

def search_orders(lambdas,intensities):
    
    orders_lambdas = []
    orders_intensities = []
    i=0
    
    while ( i <= len(lambdas)-2 ) :
        current_order_lambdas = []
        current_order_intensities = []
        while ( i<= len(lambdas)-2 and lambdas[i] - lambdas[i+1] <= 50 ):
            current_order_lambdas.append(lambdas[i])
            current_order_intensities.append(intensities[i])
            i+=1
        else  :
            i+=1
                
        orders_lambdas.append(current_order_lambdas)
        orders_intensities.append(current_order_intensities)

        
    return(orders_lambdas,orders_intensities)
        
##       


# For each order, we need to convert the spectrum into a table of slits



## 
#In order to find the splits, we need to know the offset, and we are gonna soustract the #offset to the original data to obtain the pure slits data. (à reformuler). To compute the #offset data, we need to find all the minima and then interpolate between to minima. By doing #this we will obtain a list of intensities representing the offset.

def find_offset(lambdas, intensities):
    
    offset = [0*len(intensities)]
    i_minima = []
    l_minima = []
    minima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
        if ( intensities[i-1] >= intensities[i] and intensities[i] <= intensities[i+1] ) :  
            
            minima_indices.append(i)
            i_minima.append(intensities[i])
            l_minima.append(lambdas[i])
        
    # Now that we have our minima, we need to fill the list with interpolated values between the minima in order to reach the lenght of the intensities list...
    
    for k in range(len(l_minima)-1):
            
        delta_lambda  = l_minima[k+1] - l_minima[k]
        coeff_dir = (1/delta_lambda)*(i_minima[k+1]-i_minima[k])
        for i in range(minima_indices[k],minima_indices[k+1]):
            offset.insert(i,  i_minima[k]+coeff_dir*(lambdas[i]-l_minima[k] ) )
        
    for i in range(0,minima_indices[0]-1):
        offset.insert(i,intensities[i])
        
    for i in range(minima_indices[-1],len(intensities)):
         offset.insert(i,intensities[i]) 
        
    return(offset)    


def find_offsetv2(lambdas, intensities):
    
    offset = [0*len(intensities)]
    i_minima = []
    l_minima = []
    minima_indices = []
    
    
    average = (np.sum(intensities))/(len(intensities))
    
    for i in range ( 2 , len(lambdas)-1 ):
        
        if ( intensities[i-1] >= intensities[i] and intensities[i] <= intensities[i+1] and intensities[i] <= 3*average ) :  
            
            minima_indices.append(i)
            i_minima.append(intensities[i])
            l_minima.append(lambdas[i])
    
    # Some minima are too high because of blended slits, we can delete it from the list


    # Now that we have our minima, we need to fill the list with interpolated values between the minima in order to reach the lenght of the intensities list...
    
    for k in range(len(l_minima)-1):
            
        delta_lambda  = l_minima[k+1] - l_minima[k]
        coeff_dir = (1/delta_lambda)*(i_minima[k+1]-i_minima[k])
        for i in range(minima_indices[k],minima_indices[k+1]):
            offset.insert(i,  i_minima[k]+coeff_dir*(lambdas[i]-l_minima[k] ) )
        
    for i in range(0,minima_indices[0]-1):
        offset.insert(i,intensities[i])
        
    for i in range(minima_indices[-1],len(intensities)):
         offset.insert(i,intensities[i]) 
        
    return(offset)    
    
    


            
## 
# We normalize our spectrum using the offsrt we just computed

def normalize_offset(intensities,offset):
    
    spectrum = [intensities[i] - offset[i] for i in range(len(intensities)) ]
    
    return(spectrum)
    
## 
#  Usefull algorithms



def compute_order(n,switch) :
    
    orders = search_orders(lambdas, intensities)
    order_lambdas = orders[0][n]
    order_intensities = orders[1][n]
    
    offsetv1 = find_offset(order_lambdas, order_intensities)
    offsetv2 = find_offsetv2(order_lambdas, order_intensities)
    
    spectrumv1 = normalize_offset(order_intensities, offsetv1)
    spectrumv2 = normalize_offset(order_intensities, offsetv2)
    
    if switch == 1 :
        offset = offsetv1
        spectrum = spectrumv1
    
    if switch == 2 :
        offset = offsetv2
        spectrum = spectrumv2
    
    return(order_lambdas,spectrum)
    

def print_order( n, switch):
    
    orders = search_orders(lambdas, intensities)
    order_lambdas = orders[0][n]
    order_intensities = orders[1][n]
    
    offsetv1 = find_offset(order_lambdas, order_intensities)
    offsetv2 = find_offsetv2(order_lambdas, order_intensities)
    
    spectrumv1 = normalize_offset(order_intensities, offsetv1)
    spectrumv2 = normalize_offset(order_intensities, offsetv2)
    
    if switch == 1 :
        offset = offsetv1
        spectrum = spectrumv1
    
    if switch == 2 :
        offset = offsetv2
        spectrum = spectrumv2
    
    plt.plot(order_lambdas, order_intensities, 'red')
    plt.plot(order_lambdas, offset, 'black')
    plt.plot(order_lambdas, spectrum, 'blue')
    plt.show()




##
#Function which finds the wavelenghts corresponding to the local maxima of the intensities and the associated intensities values in order to find the slits locations 
##

def find_maxima(lambdas,intensities):
    
    intensities_maxima = []
    lambdas_maxima = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] ) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
        
    return(lambdas_maxima,intensities_maxima)   
    
     
        
def find_order_maxima(n,switch):
    
    lambdas = compute_order(n,switch)[0]
    intensities = compute_order(n,switch)[1]
    
    return(find_maxima(lambdas,intensities))

##
# In this function, we search for each list of lambdas and intensities corresponding to each maximum in order to fit a gaussian with this data and find the exact localisation of the slit which corresponds to the maximum we study
    
def find_slits(lambdas,intensities):
    
    slits = []
    
    intensities_maxima = []
    lambdas_maxima = []
    maxima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] ) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
            maxima_indices.append(i)
    # Now we need to find the slit around each maximum
    
    for j in range(len(lambdas_maxima)):
        
        local_slit = []
        
        # left minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index-1] ):
            local_slit.append([lambdas[index],intensities[index]])
            index -= 1
        local_slit.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        # right minimum research
        index = maxima_indices[j]
        while ( intensities[index] > intensities[index+1] ):
            local_slit.append([lambdas[index],intensities[index]])
            index += 1
        local_slit.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        local_slit.sort() # We sort the list according to the lambdas order
        local_slit_lambdas = [ local_slit[i][0] for i in range(len(local_slit)) ]
        local_slit_intensities = [ local_slit[i][1] for i in range(len(local_slit)) ]
        slits.append([local_slit_lambdas,local_slit_intensities])
    
    return(slits)
    
# Same function but for the order n :
    
def find_slits_order(n):
    
    lambdas = compute_order(n,2)[0]
    intensities = compute_order(n,2)[1]
    slits = find_slits(lambdas,intensities)
    
    for i in range(len(slits)):
        
        x = slits[i][0]
        y = slits[i][1]
        #plt.plot(x,y)
    #plt.show()   
    
    return(slits)
    
    
##
# Functions which fits a gaussian model on the given slit to dertermine its center with the higher posible precision given the data we have

def find_slit_center(lambdas,data):
    y = np.copy(data)-np.min(data)
    X = np.arange(len(y))
    centre = float(np.sum(X*y))/np.sum(y) # centre of the gaussian (= avergae)
    ampl = np.max(y) # amplitude of the gaussian
    width = np.sqrt(np.abs(np.sum((X-centre)**2*y)/np.sum(y))) # width of the gaussian
    

    def gaussian(x,amp,cen,wid):
        return(amp*np.exp(-(x-cen)**2/wid))
    
    
    try:
        # we use the lmfit algorithm to improve our fit's precision
        gaussian_model = lmfit.Model(gaussian)
        params = gaussian_model.make_params(cen=centre,amp=ampl,wid=width)
        result = gaussian_model.fit(y,params,x=X)
        [cen,wid,amp] = list(result.best_values.values())
        if 0 < cen < len(y) :
            centre = cen
    except :
        pass
    # We have the center on a normalized scale which start from 0 and goes to len(y), we need to bring it back to the real lambda value of the center of the slit
    
    lambda_min = lambdas[0]
    lambda_max = lambdas[len(lambdas)-1]
    step = (lambda_max - lambda_min)/len(lambdas)
    lambda_centre = lambda_min + centre*step

    return(lambda_centre)
            
    
##
# Function which finds all the slit's centers and return the list of wavelengths corresponding to each slit (that is to say each maximum) for a given order

def slit_center_listing(n):
    
    slits = find_slits_order(n)
    slit_centers = []
    
    for i in range(len(slits)):
        
        slit_center = find_slit_center(slits[i][0],slits[i][1])
        slit_centers.append(slit_center)
    
    return(slit_centers)



##
 # Prints the selected order of the pkl file


def find_arturo_slit_order(n) :
    
    return(search_orders(lambdas2,intensities2)[0][n])
    
def print_arturo_slits_order():
    
    for n in range(36) :
        orders = search_orders(lambdas2,intensities2)
        x , y = orders[0][n],orders[1][n]
        x.sort()
        plt.plot(x,y)
        plt.show()
        
        


##
# This function will compare two lists of numbers in order to find the better possible matching for each value of each list with the values of the other lists and thus determine the correlation rate between those two lists.

def data_matching(x,y,precision):
    
    best_matching_values = []
    
    for i in range(len(x)):
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        best_matching_values.insert(i,best_matching_value)
        print(x[i],best_matching_value)
    while (best_matching_values.__contains__(None)):
        best_matching_values.remove(None)

    return(float(len(best_matching_values))/float(len(x)))


def order_matching(n,precision):
    
    x = slit_center_listing(n)
    y = select_atlas(x[0],x[len(x)-1])
    
    return(data_matching(x,y,precision))


##
#

path = '/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/thar_UVES_MM090311.dat'

file = open(path)    
wavelength_slits = []
for i in range(0,3005):
     
    line = file.readline()
    line_lambda = float(line[12:23])
    wavelength_slits.append(line_lambda)
        
    
for i in range(3006,3223):
       
    line = file.readline()
    line_lambda = float(line[11:23])
    wavelength_slits.append(line_lambda)

file = open('/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/ThAr_Atlas.pkl','w')
ThAr_pickle = pickle.dump(wavelength_slits,file)
file.close()


def select_atlas(lmin,lmax):
    wavelength_slits= cPickle.load(open('/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/ThAr_Atlas.pkl','rb'))
    selected_wavelength_slits = []
    for slit in wavelength_slits :
        if lmin <= slit <= lmax:
            selected_wavelength_slits.append(slit)
    
    
    return(selected_wavelength_slits)
    





    