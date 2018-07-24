import pyfits
import cPickle
import numpy as np
from scipy import sparse
import collections
import matplotlib.pyplot as plt
import math as m
import os
os.chdir(r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr')
from find_ThAr_slits import find_ThAr_slits
from create_B_matrix import create_B_matrix

def full_A(test):
    
    for i in range(0,39):

        path = r'C:\Users\Martin\Documents\Stage_IRAP_2018\NeoNarval\NeoNarval\Martin_Jenner\test_ThAr\Bmatrix_data_sheet.txt'
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
            
        dic['order'] = i
        with open(path, 'w') as file:
            for value in dic.items():
                line = str(value[0]) + " : " + str(value[1])+ "\n"
                file.write(line)
        print('================================')       
        print('computing matrix A : order {0}'.format(i))
        find_ThAr_slits(test)
        print('searching slits data : done')
        create_B_matrix(test)
        print('assembling matrix : done')
        print('================================')

full_A('ThAr')