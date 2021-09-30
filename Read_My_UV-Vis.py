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
fname = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P69_Concentrated_CF_in_HEPES_Stock\2021-09-20 UV Vis\2021-09-20_CF_in_HEPES_dilutions.xlsx"

baseline_wavelength = 800 #nm

calculate_CF_concentration = True
Dilution_Factors_SheetName = 'Sheet2'

#######################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
 
#######################################  


def makeplot(xdata, ydata, headers, title='', xlims=[300,600], ylims=[-0.1,1]):
    
    for i in np.arange(0,len(ydata[0])):
        ax.plot(xdata, ydata[:,i], label=headers[i])
        
    ax.legend(bbox_to_anchor=(1.1, 1.05))
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorbance')
    ax.set_title(title)
    
    plt.xlim(xlims)
    plt.ylim(ylims)

#######################################   
#%%
# Convert to DataFrame, extract headers
absorbance = pd.read_excel(fname)
headers = list(absorbance.columns)

# Convert to numpy array, extract wavelength data, reassign variable 'absorbance' to be just absorbance data
absorbance = absorbance.to_numpy()
wavelengths = absorbance[:,0]

# remove Wavelength column from both headers and absorbance data
headers = headers[1:]
absorbance = absorbance[:,1:]

# Plot raw data with legend
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance, headers, xlims=[200,600], title = 'UV-Vis')


 ####################################### 
#%%
# Create baseline subtracted data and plot the same graph

baseline_index = np.where(wavelengths == baseline_wavelength)[0][0]

absorbance_BC = absorbance - absorbance[baseline_index]


# Plot baselined data with legend
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance_BC, headers, xlims=[200,600], title = 'UV-Vis Baselined')


#######################################   
#%%
# Calculate CF Concentration from linear regression data

if calculate_CF_concentration == True:
    
    " Beer's Law: A = e*L*C "
    
    CF_ext_coefficient = 74000	# M^-1 cm^-1
    pathlength = 1 # cm
    
    wv_450 = np.where(wavelengths == 450)[0][0]
    wv_500 = np.where(wavelengths == 500)[0][0]
    
    # find max absorbances for baselined sample
    abs_max = np.max(absorbance_BC[wv_450:wv_500], axis=0)
    
    # pull out dilution factors from excel shet
    DF = pd.read_excel(fname, sheet_name=Dilution_Factors_SheetName)
    DF = DF.to_numpy()
    DF = DF[0][1:]
    
    # scatter plot of maximum absorbance of baselined sample to dilution factor
    fig, ax = plt.subplots()
    ax.scatter(DF,abs_max)
    
    # linear regression
    slope = np.polyfit(DF.astype(float), abs_max.astype(float), 1)[0]
    
    Xvals = np.linspace(0,np.max(DF))
    Yvals = Xvals * slope
    
    ax.set_xlabel('Dilution Factor')
    ax.set_ylabel('Absorbance')
    ax.set_title('Calculating CF Concentration')
    ax.plot(Xvals, Yvals, label= "y = {}x".format(np.round(slope, decimals=0)))
    fig.legend()
    
    # At dilution factor 1, absorbance is expected to be:
    abs_at_DF_1 = slope
    # Then, A = e*l*c, or, c = A/(e*l)
    conc_M = abs_at_DF_1 / (CF_ext_coefficient * pathlength) # Molar
    conc_mM = conc_M * 10**3 # mM
    
    plt.annotate("Stock Concentration = {} mM".format(np.round(conc_mM, decimals = 1)), [np.min(DF), 0.8*np.max(abs_max)])




























