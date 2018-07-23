import matplotlib.pyplot as plt
import numpy as np
import pickle
import numpy.ma as ma
import cPickle as pickle
from scipy import interpolate
import pyfits
import sys
import utils
import os





dirBrut='/home/main/Documents/Arthur/Brut'

dirs=[f for f in os.listdir(dirBrut) if f.endswith('17') or f.endswith('16')]

for d in dirs:
	files=[f for f in os.listdir(dirBrut+'/'+d+'/BRUT') if f.endswith('c.fits')]
	for f in files:
		fin=pyfits.open(dirBrut+'/'+d+'/BRUT/'+f)
		data=fin[0].data


