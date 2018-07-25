#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


file='extraction_methods/Extraction_files/betcrb_001.s'

f=open(file,'r')
f.readline()
Nwvl,a=map(float,f.readline().split())
Wvl=np.zeros(int(Nwvl))
Sp=np.zeros(int(Nwvl))
for i,l in enumerate(f.readlines()):
	Wvl[i],a,Sp[i],=map(float,l.split())
plt.plot(Wvl,Sp)
plt.show()