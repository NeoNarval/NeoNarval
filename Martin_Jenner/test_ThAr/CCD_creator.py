import numpy as np
import collections
import matplotlib.pyplot as plt
import pyfits
import cPickle
import numpy.random as rd
import matplotlib.pyplot as plt
import os

def rint(x) : return int(round(x))

def randomizer(nb_slits):
    lenCCD = end_index-ini_index
    CCD = 0.0001*np.ones((4000,4000))
    slits = []
    while len(slits)<nb_slits:
        index = rd.randint(0,lenCCD)
        flag = False
        for s in slits:
            if s==index:
                flag = True
        if flag == False:
            slits.append(index)
    slits.sort()
    print(slits)
    for s in slits:
        for y in range(thickness):
            CCD[s+ini_index][rint(y+left_lane[s+ini_index])] = 1
            
    return CCD
  
    
def gaussian(x,amp,wid):
    return(amp*np.exp(-(x)**2/(2*wid**2)))
    
    
def gaussianizer(nb_slits):
    window = 6
    lenCCD = end_index-ini_index
    CCD = 30000*np.ones((4000,4000))
    slits = []
    while len(slits)<nb_slits:
        index = rd.randint(0,lenCCD)
        flag = False
        for s in slits:
            if s==index:
                flag = True
        if flag == False:
            slits.append(index)

    slits.sort()
    print(slits)
    for s in slits:
        # wid = rd.uniform(1,5)
        wid = 2
        # amp = 1
        
        #--------sans repartion aleatoire des amplitudes-------
        # amp = rd.randint(0,175000/thickness)
        # for y in range(thickness):
        #     CCD[s+ini_index][rint(y+left_lane[s+ini_index])]+=gaussian(0,amp,wid)
        #     for w in range(1,window):           
        #         CCD[s+ini_index+w][rint(y+left_lane[s+ini_index])]+=gaussian(w,amp,wid)
        #         CCD[s+ini_index-w][rint(y+left_lane[s+ini_index])]+=gaussian(-w,amp,wid)
        #-----------------------------------------------------
        
        #--------------------avec-------------------------
        # amp = rd.randint(0,175000)
        # 
        # sum = gaussian(0,amp,wid)
        # for y in range(thickness-1):
        #     value = rd.uniform(0,sum)
        #     CCD[s+ini_index][rint(y+left_lane[s+ini_index])]+=value
        #     sum-=value
        # y = thickness-1
        # CCD[s+ini_index][rint(y+left_lane[s+ini_index])]+=sum
        # for w in range(1,window):
        #     sum_plus = gaussian(w,amp,wid)
        #     sum_moins = gaussian(-w,amp,wid)
        #     for y in range(thickness-1):
        #         value_plus = rd.uniform(0,sum_plus)
        #         value_moins = rd.uniform(0,sum_moins)
        #         CCD[s+ini_index+w][rint(y+left_lane[s+ini_index])]+=value_plus
        #         CCD[s+ini_index-w][rint(y+left_lane[s+ini_index])]+=value_moins
        #         sum_plus-=value_plus
        #         sum_moins-=value_moins
        #     y = thickness-1
        #     CCD[s+ini_index+w][rint(y+left_lane[s+ini_index])]+=sum_plus
        #     CCD[s+ini_index-w][rint(y+left_lane[s+ini_index])]+=sum_moins
        #-------------------------------------------------- 
        
        # #--------avec répartition gaussienne de l'amplitude--------
        # amp = rd.randint(0,175000/thickness)
        # sig = 4
        # for y in range(thickness):
        #     coeff = gaussian(y-(thickness//2),1,sig)
        #     CCD[s+ini_index][rint(y+left_lane[s+ini_index])]+=gaussian(0,amp,wid)*coeff
        #     for w in range(1,window):           
        #         CCD[s+ini_index+w][rint(y+left_lane[s+ini_index])]+=gaussian(w,amp,wid)*coeff
        #         CCD[s+ini_index-w][rint(y+left_lane[s+ini_index])]+=gaussian(-w,amp,wid)*coeff
        # 
        #-----------------------------------
        val_col = rd.randint(0,175000)
        sig = 10
        sum = gaussian(0,val_col,wid)
        amp = sum/thickness
        err = np.sqrt(amp)
        for y in range(thickness):
            coeff = gaussian(y-(thickness//2),1,sig)
            CCD[s+ini_index][rint(y+left_lane[s+ini_index])]+= err*rd.randn()+amp*coeff
        
        for w in range(1,window):
            sum_plus = gaussian(w,val_col,wid)
            sum_moins = gaussian(-w,val_col,wid)
            amp_plus = sum_plus/thickness
            amp_moins = sum_moins/thickness
            err_plus=np.sqrt(amp_plus)
            err_moins = np.sqrt(amp_moins)
            for y in range(thickness):
                coeff = gaussian(y-(thickness//2),1,sig)
                CCD[s+ini_index+w][rint(y+left_lane[s+ini_index])]+= err_plus*rd.randn() + amp_plus*coeff
                CCD[s+ini_index-w][rint(y+left_lane[s+ini_index])]+= err_moins*rd.randn() +amp_moins*coeff
            
        
            
    return CCD
            
def regulizer(nb_slits):
    lenCCD = end_index-ini_index
    CCD = 0.0001*np.ones((4000,4000))
    delta_slits = rint(lenCCD/nb_slits)
    for i in range(nb_slits):
        s = i*delta_slits
        print(s)
        for y in range(thickness):
           CCD[s + ini_index][rint(y+left_lane[s + ini_index])] = 1
    return CCD

def CCD_creator(nb_slits, test):
    path = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt'
    global ini_index       # first index of the window of ThAr to be processed
    global end_index       # last index of the window of ThAr to be processed
    global envelope_data    # Upper envelope of each lane
    global thickness_data   # Thickness of the lanes
    global lane             # considered lane
    global nbr_lanes        # Number of lanes per order
    global left_lane
    global thickness
    
    # Import of data from data file
    dic = collections.OrderedDict()
    with open(path, 'r') as file:
        content = file.readlines()

    content = [line.strip() for line in content]
    for line in content:
        param = line.split(" : ")
        if len(param) == 2:
            nom_param = param[0]
            value_param = param[1]
            dic[nom_param] = value_param
                 
    envelope_data_file  = dic["Lane envelope file"]
    thickness_data_file = dic["Lane thickness file"]
    ini_index           = int(dic["initial index"])
    end_index           = int(dic["final index"])
    order               = int(dic["order"])
    nbr_lanes           = int(dic["nb lane per order"])
    lane                = int(dic["lane"])
    envelope_data  = cPickle.load(open(envelope_data_file, 'r'))
    thickness_data = cPickle.load(open(thickness_data_file, 'r'))
    thickness_float = thickness_data[nbr_lanes * order + lane - 1]
    thickness       = rint(thickness_float)   # The thickness of the considered lane
    first_pix       = thickness*ini_index
    
    left_lane     = envelope_data[:,nbr_lanes*order+lane-1]
    print(thickness)
    print(first_pix)
    if test == 'rand':
        CCD = randomizer(nb_slits)
    elif test == 'reg':
        CCD = regulizer(nb_slits)
    elif test == 'gauss':
        CCD = gaussianizer(nb_slits)
    print(np.shape(CCD))
    pklpath = r'C:\Users\Martin\Documents\Stage IRAP 2018\NeoNarval\TEMP_\random_CCD_'+str(ini_index)+'-'+str(end_index)
    file = open(pklpath, 'w')    
    cPickle.dump(CCD, file)
    file.close()
    dic["test file"] = pklpath
    with open(path, 'w') as file:
        for value in dic.items():
            line = value[0] + " : " + value[1] + "\n"
            file.write(line)
    plt.matshow(CCD, aspect='auto')
    plt.show()
CCD_creator(10, 'gauss')
