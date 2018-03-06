"""
Documentation and blabla
IN : fichier calibered.fits (fichier dépouillé avec une série de valeurs de lambdas et leurs intensités, toutes les voies sont mêlées)

OUT : un tableau avec une liste longueurs d'ondes de raies à comparer à un atlas
"""

import pyfits
import numpy as np
import matplotlib.pyplot as plt


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
        while ( i<= len(lambdas)-2 and lambdas[i] < lambdas[i+1] ):
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
#Function which finds the wavelenghts corresponding to the local maxima of the intensities and the associated intensities values in order to find the slits locations 
##

def find_maxima_v1(lambdas,intensities):
    
    intensities_maxima = []
    lambdas_maxima = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their                      intensities, but we are gonna detect some local maxima which are not slits, so we have to try to discriminate those false slits. We are gonna  find the local average value of the offset in intensity (which is not zero, and is not always the same) and then compare the local maximum value with this local offset : if the maximum is a slit, it will be far bigger than the offset, if not we juste delete it from the maxima list.
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] ) :    
        
            # lolcal offset computation
            nb = 10
            if ( i <= nb ) :
                local_offset = (1/nb)*np.sum([intensities[k] for k in range(0,nb)])
            if ( i >= len(intensities) - nb ):
                local_offset = (1/nb)*np.sum([intensities[k] for k in range(len(intensities)-nb, len(intensities)) ] )
            else :
                local_offset = (1/nb)*np.sum( [ intensities[k] for k in range( i-int(round(nb/2)) , i+int(round(nb/2)) )] )      
            
            # additionnal discrimination
            if ( (intensities[i] >= 2*local_offset ) and  ( intensities[i] >= 1,5*(intensities[i-1]+intensities[i+1])/2 ) ) :
                intensities_maxima.append(intensities[i])
                lambdas_maxima.append(lambdas[i])
        
    return(lambdas_maxima,intensities_maxima)    
        

## 
#In order to find the splits, we need to know the offset, and we are gonna soustract the #offset to the original data to obtain the pure slits data. (à reformuler). To compute the #offset data, we need to find all the minima and then interpolate between to minima. By doing #this we will obtain a list of intensities representing the offset.

def find_offset(lambdas, intensities):
    
    offset = [0*len(intensities)]
    i_minima = []
    l_minima = []
    minima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
        if ( intensities[i-1] > intensities[i] and intensities[i] < intensities[i+1] ) :  
            
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
        
        if ( intensities[i-1] > intensities[i] and intensities[i] < intensities[i+1] and intensities[i] <= 1.5*average ) :  
            
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
#  Usefull algorithm

def print_order(lambdas, intensities ,n, switch):
    
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
# Now that we have a normalized spectrum, we can fit a gaussian for each maximum and determine its center in order to finds the precise locations of the slits. The idea is to do the same thing for all maxima using the local wavelength and intensities data.


def find_center(lambdas,intensities):
    
    total = np.sum(intensities)
    
























    