#!/usr/bin/env python
# -*- coding: utf-8 -*-


import cPickle
import pickle

"""
This function reads a pickle file, returning its content. The given argument is the path of the file in the directory.
Input : 
- path : string, name of the file to read.
Output :
- data : extracted data.
"""

def read_pickle(path):
    
    file = open(path,'r')
    data = cPickle.load(file)
    file.close()
    return(data)