# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 17:04:20 2021

@author: Alison


    This script takes an excel file of UV-Vis data from the Landry Lab UV-Vis
    machine, baseline subtracts to a certain wavelength if given, and returns 
    the spectral graphs.
    
    If the UV-Vis graphs are used for determining concentation of CF or another
    molecule, this script can also create a linear regression of absorbance
    versus dilution factor and return the bulk concencentration.
"""
fname = 

baseline_wavelength =  #nm

calculate_CF_concentration = True

#######################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
 
#######################################   

def extractdata(wdir, fname, fluor_index):
    
    " Pull out data for one csv file"
    
    # extract header title
    base = fname[:-7]
    title = os.path.basename(base)
    title = title.split(' -  ')
    title = title[0]
 
   # Get Data
    os.chdir(wdir)
    
    # how many fluorophore files?
    N_fluor = len(fluor_index)
    
    # extract the first emission data so we know what size to make the array
    base = fname[:-7]
    df = pd.read_csv(fname)
    data = df.to_numpy()
    H, N = np.shape(data)
    
    # put fluorophore data into empty emission data array
    data = np.empty((H,N,N_fluor))
    for x in range(0,N_fluor):
        i = int(fluor_index[x])
        temp = pd.read_csv(base + fluor_name[i] + ".csv")
        temp = temp.to_numpy()
        data[:,:,x] = temp
    
    # Emission data sits in 3rd and later columns
    data = data[:,2:,:]
    H,N,N_fluor = np.shape(data)

    return H, N, N_fluor, data, title