import pyfits
import matplotlib.pyplot as plt

"""
read .fits file encoding a spectrum : 1 column for the wavelength another for the intensity

    : path : string with the relative path toward the file
"""
path =r"C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\DRS\DATA\th_calibered.fits"

def order_gen(lambd, intensity):
    L = len(lambd)
    k = 1
    order_l = []
    order_i = []
    current_order_l = [lambd[0]]
    current_order_i = [intensity[0]]
    for k in range (L):
        if lambd[k-1]<lambd[k] and k<L-2:
            current_order_l.append(lambd[k])
            current_order_i.append(intensity[k])
        else :
            order_l.append(list(current_order_l))
            order_i.append(list(current_order_i))
            current_order_l = [lambd[k]]
            current_order_i = [intensity[k]]
            
    return order_l, order_i
    
def plot_all(lambd, intensity):
    plt.plot(lambd, intensity)
    plt.xlabel(r'wavelength ($\AA$)')
    plt.ylabel('intensity')
    plt.show()
    
def plot_orders(lambd, intensity, n):
    order_l, order_i = order_gen(lambd, intensity)
    plt.plot(order_l[n], order_i[n])
    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Intensity')
    plt.title('order n°{0}'.format(n))
    plt.show()
        
        
def fits_reader(path):
    l = pyfits.open(path)
    a = l[1].data
    l.close()
    lambd = a['wavelength_lane1']
    intensity = a['intensity_lane1']
    return lambd, intensity
 
lambd, intensity = fits_reader(path)
order_l, order_i = order_gen(lambd,intensity)
print(len(order_l))
#print(order_i[0])
#plot_all(lambd, intensity)
plot_orders(lambd, intensity, 35)
