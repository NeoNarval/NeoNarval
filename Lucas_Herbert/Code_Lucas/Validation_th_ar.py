import pyfits
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import cPickle
import pickle
import lmfit 
from scipy.optimize import leastsq
from scipy.optimize import curve_fit

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
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] and intensities[i] >= 0  and intensities[i] >= 0.5) :    
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
    
def find_center_order(n):
    
    lambdas = compute_order(n,2)[0]
    intensities = compute_order(n,2)[1]
    plt.figure(1)
    plt.plot(lambdas,intensities,color='black')
    slits = find_slits(lambdas,intensities)
    
    for i in range(len(slits)):
        
        x = slits[i][0]
        y = slits[i][1]
        plt.figure(1)
        plt.plot(x,y, color='blue')
    plt.show()   
    
    return(slits)
    
    
    
##
# Function which finds all the slit's centers and return the list of wavelengths corresponding to each slit (that is to say each maximum) for a given order

def lambda_center_listing(n):
    
    slits = find_center_order(n)
    lambdas_data = []
    
    for i in range(len(slits)):
        
        lambdas_data.insert(i,fitting_slit(slits[i][0],slits[i][1]))
        
    return(lambdas_data)

def computed_order_matching(n,precision):
    
    lambdas_data = lambda_center_listing(n)
    lambdas_centers = []
    lambdas_widths = []
    for i in range(len(lambdas_data)):
        lambdas_centers.insert(i,lambdas_data[i][0])
        lambdas_widths.insert(i,lambdas_data[i][1])
    atlas_lambdas = select_atlas(min(lambdas_centers),max(lambdas_centers))
    data = data_matching(lambdas_centers,atlas_lambdas,precision)
    gaps = []
    lambdas_gaps = []
    matching_data = data[0]
    for i in range(len(matching_data)):
        gaps.insert(i,matching_data[i][2])
        lambdas_gaps.insert(i,matching_data[i][0])
    plt.figure(2)
    plt.bar(lambdas_gaps,gaps, color='brown')
    plt.title("Error = f(Lambda)")
    plt.show() 
    average_error = float(np.sum(gaps)/len(matching_data))
    average_width = float(np.sum(lambdas_widths)/len(lambdas_widths))
    print("Average width",average_width)
    return(data[1],average_error)
    
def global_computed_matching(precision):
    
    for n in range(36):
        print(computed_order_matching(n,precision))
    
##
 # Prints the selected order of the pkl file


def find_arturo_slit_order(n) :
    path2 = '/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/ThAr_Line_Cal.pkl'
    lambdas2 , intensities2 = open_pkl(path2)['Lines'], open_pkl(path2)['FWHM']
    return(search_orders(lambdas2,intensities2)[0][n])
    
def find_arturo_data_order(n):
     return(search_orders(lambdas2,intensities2)[0][n],search_orders(lambdas2,intensities2)[1][n])
     
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
    
    matching_data = []
    
    for i in range(len(x)):
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        data = [x[i],best_matching_value,gap_min]            
        if (best_matching_value != None ):
            matching_data.insert(i,data)
    
    matching_rate = float(len(matching_data))/float(len(x))
    print("Matched, Total",len(matching_data),len(x))
    return(matching_data,matching_rate)

# Does the same but applied to the given order, comparing our finding of the slits and the atlas

def order_matching(n,precision):
    
    x = slit_center_listing(n)
    y = select_atlas(x[0],x[len(x)-1])
    
    return(data_matching(x,y,precision),data_matching(y,x,precision))
    
def global_matching(precision):
    
    matching1 = []
    matching2 = []
    
    for n in range(36):
            matching1.append(order_matching(n,precision)[0])
            matching2.append(order_matching(n,precision)[1])
    global_matching1 = np.sum(matching1)/float(len(matching1))  
    global_matching2 = np.sum(matching2)/float(len(matching2))  
    return(global_matching1,global_matching2)    
        
    

# Same idea on Arturo's slits selection and the atlas given an precise order

def arturo_order_matching(n,precision):
    
    x = find_arturo_slit_order(n)
    y = select_atlas(x[0],x[len(x)-1])
    
    return(data_matching(x,y,precision),data_matching(y,x,precision))

def arturo_global_matching(precision):
    
    matching1 = []
    matching2 = []
    for n in range(36):
        matching1.append(arturo_order_matching(n,precision)[0])
        matching2.append(arturo_order_matching(n,precision)[1])
    global_matching1 = np.sum(matching1)/float(len(matching1))  
    global_matching2 = np.sum(matching2)/float(len(matching2))  
    return(global_matching1,global_matching2)    

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
    




def print_atlas():

    path = '/home/stagiaire/depot_git/NeoNarval/Lucas_Herbert/Documents/thar_UVES_MM090311.dat'
    
    file = open(path)    
    wavelength_slits = []
    intensities = []
    for i in range(0,3005):
        
        line = file.readline()
        line_lambda = float(line[12:23])
        if line[25] == "-" :
            intensities.append(float(line[25:29])*-1)
        else :
            intensities.append(float(line[25:29]))
            
        wavelength_slits.append(line_lambda)
    plt.plot(wavelength_slits,intensities)
    
    plt.show()
    
    

##
# Localizes the slits found by Arturo using the same method as the one used before.
# Firstly, we search for the maxima corresponding to Arturo's slits and then we do the same job as before with each maximum


def localize_arturo_slits(n):
    
    slits = []
    lambdas = compute_order(n,2)[0]
    intensities = compute_order(n,2)[1]
    plt.figure(1)
    plt.plot(lambdas,intensities,'blue', color = 'black')
    intensities_maxima = []
    lambdas_maxima = []
    maxima_indices = []
    
    for i in range ( 2 , len(lambdas)-1 ):
        
    # We use the mathematical definition of a maximum to find the maxima and their intensities 
        
        if ( intensities[i-1] < intensities[i] and intensities[i] > intensities[i+1] and intensities[i] >= 0 ) :    
            intensities_maxima.append(intensities[i])
            lambdas_maxima.append(lambdas[i])
            maxima_indices.append(i)

            
    # Now we need to find the slit around each maximum which has also been found by Arturo, so we keep only the maxima found by Arturo. Thus, we have to delete the other maxima and their indices from their lists if they don't match a maximum found by arturo
    
    arturo_slits = find_arturo_slit_order(n)
    selected_maxima_indices = []

    for i in range(len(arturo_slits)) :
        
        plt.axvline(arturo_slits[i])
        
        gap_min = 0.5
        matching_maximum = None
        index = 0
        for j in range(len(lambdas_maxima)) :
            
            if ( abs(arturo_slits[i] - lambdas_maxima[j]) <= gap_min ):
                gap_min =  abs(arturo_slits[i] - lambdas_maxima[j])
                matching_maximum = lambdas_maxima[j]
                matching_maximum_index = maxima_indices[j]
        if not matching_maximum == None :
            selected_maxima_indices.append(matching_maximum_index)
            
            
    if ( len(arturo_slits) != len(selected_maxima_indices) ) :
        print("PROBLEM : Some of Arturo's lambdas are unmatched!")
    

    for j in range(len(selected_maxima_indices)):
        
        local_slit = []
        
        # left minimum research
        index = selected_maxima_indices[j]
        while ( intensities[index] > intensities[index-1] ):
            local_slit.append([lambdas[index],intensities[index]])
            index -= 1
        local_slit.append([lambdas[index],intensities[index]]) # don't forget the last point
        # right minimum research
        index = selected_maxima_indices[j]
        while ( intensities[index] > intensities[index+1] ):
            local_slit.append([lambdas[index],intensities[index]])
            index += 1
        local_slit.append([lambdas[index],intensities[index]]) # don't forget the last point
        
        local_slit.sort() # We sort the list according to the lambdas order
        local_slit_lambdas = [ local_slit[i][0] for i in range(len(local_slit)) ]
        local_slit_intensities = [ local_slit[i][1] for i in range(len(local_slit)) ]
        slits.append([local_slit_lambdas,local_slit_intensities])
    
    for i in range(len(slits)):
        x = slits[i][0]
        y = slits[i][1]
        plt.plot(x,y,color = 'blue')
    
    plt.show() 
        
    return(slits)

    
##
# Function fitting a gaussian  to the given data :

def fitting_slit(lambdas,data):
    
    y = np.copy(data)
    X = lambdas
    def gaussian(x,cen,amp,wid):
        return(amp*np.exp(-(x-cen)**2/(2*wid**2)))
    
    # printing the initial gaussian (before the optimization of the fit)
    naive_center = float(np.sum(lambdas*y))/np.sum(y)
    naive_width = np.sqrt(abs((np.sum((lambdas-naive_center)**2*y)/np.sum(y))))
    naive_ampl = np.max(y)
    naive_gaussian = [ gaussian(x,naive_center,naive_ampl,naive_width) for x in lambdas ]
    # plt.plot(lambdas,naive_gaussian,color='green')
    # plt.axvline(naive_center,color='green')
    # plt.show()
    lambda_centre = naive_center
    lambda_width = naive_width
    try :
        # we use the lmfit algorithm to improve our fit's precision
        gaussian_model = lmfit.Model(gaussian)
        params = gaussian_model.make_params(cen=naive_center,amp=naive_ampl,wid=naive_width)
        result = gaussian_model.fit(y,params,x=X) 
        # printing the best gaussian fit
        best_gaussian_fit = result.best_fit
        best_cen = result.best_values['cen']
        best_wid = result.best_values['wid']
        best_amp = result.best_values['amp']
        best_params_gaussian = [ gaussian(x,best_cen,best_amp,best_wid) for x in lambdas ]
        plt.figure(1)
        plt.plot(lambdas, best_params_gaussian, color='red')
        plt.axvline(best_cen, color = 'red')
        plt.show()
        #computed_centre = float(np.sum(X*best_gaussian_fit))/np.sum(best_gaussian_fit) 
        #plt.plot(lambdas, best_gaussian_fit, 'b--' ,color='purple')
        #plt.axvline(computed_centre, color='purple')
        lambda_centre = best_cen
        lambda_width = best_wid
        # we need the report data to improve our understanding of the results
        report = result.fit_report()
        
        
    except : 
        report = ""
        pass
        
    return(lambda_centre,lambda_width)
    
##
# This function fits a gaussian to each slit found previously thanks to Arturo's list

def localize_arturo_slits_centers(n) :
    
    slits = localize_arturo_slits(n)
    lambdas_data = []
    
    for i in range(len(slits)):
        
        slit_infos = []
        data = fitting_slit(slits[i][0],slits[i][1])
        lambdas_data.append(data)
        
    return(lambdas_data)   

def arturo_fitted_order_matching(n,precision):
    
    lambdas_values = []
    widths_values = []
    data = localize_arturo_slits_centers(n)

        
    for j in range(len(data)):
        lambdas_values.insert(j,data[j][0])
        widths_values.insert(j,data[j][1])
        
    average_width = float(np.sum(widths_values)/len(widths_values))
    
    print("Average width",average_width)
    return( data_matching(lambdas_values,select_atlas(min(lambdas_values),max(lambdas_values)),precision) )

def plot_order_matching(n,precision):
    
    matching = arturo_fitted_order_matching(n,precision)
    matching_data = matching[0]
    matched_lambdas = []
    gaps = []
    for i in range(len(matching_data)):
        matched_lambdas.insert(i,matching_data[i][0])
        gaps.insert(i,matching_data[i][2])
    plt.figure(2)
    plt.bar(matched_lambdas,gaps, color='brown')
    plt.title("Error = f(Lambda)")
    plt.show()   
    average_gap = float(np.sum(gaps)/len(gaps))
    print(matching[1],average_gap)
    return(matching[1],average_gap)

def plot_all_matching(precision):
    
    for n in range(36):
        plot_order_matching(n,precision)

    #     slit_infos.append(data)     # adding the center of each slit to the info list
    #     
    #     result = data[1]
    #     report = result.fit_report()
    #     #print(report)
    #     chi_square_str = report[189:198]
    #     chi_square = float(chi_square_str)
    #     slit_infos.append(chi_square)    # adding the chi square to the info list
    #     
    #     step = data[2]
    #     error_str = report[392:416]
    #     def local_reader(str):
    #         nbr_str = ''
    #         i = 0
    #         while str[i] != '-' :
    #             i+=1
    #         debut = i
    #         i+=1
    #         while str[i] != '(' :
    #             nbr_str+= str[i]    
    #             i+=1
    #         nbr = float(nbr_str)
    #         return(nbr)
    #     error = local_reader(error_str)
    #     slit_infos.append(error)
    #     
    #     
    #     slits_infos.append(slit_infos)
    #     
    #     
    # return(slits_infos)













