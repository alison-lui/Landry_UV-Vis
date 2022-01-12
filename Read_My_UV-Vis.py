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

fname = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P123\UV-Vis 12-21\2021-12-21 Free CF Samples.xlsx"

baseline_wavelength = [550, 600] #nm

calculate_Stock_CF_concentration = False
calculate_All_CF_concentrations = True
Dilution_Factors_SheetName = 'Sheet2' # must be in same file as fname

#######################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from tabulate import tabulate
import os
import sys
 
#######################################  


def makeplot(xdata, ydata, headers, title='', xlims=[300,600], ylims=[-0.1,1]):
    
    for i in np.arange(0,len(ydata[0])):
        ax.plot(xdata, ydata[:,i], label=headers[i])
        
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorbance')
    ax.set_title(title)
    
    plt.xlim(xlims)
    plt.ylim(ylims)
    
    plt.tight_layout()

#######################################   
#%%
" FORMAT DATA "

# Extract directory name
dirname = os.path.dirname(fname)

# Convert to DataFrame, extract headers
absorbance = pd.read_excel(fname)
headers = list(absorbance.columns)

# Convert to numpy array, extract wavelength data, reassign variable 'absorbance' to be just absorbance data
absorbance = absorbance.to_numpy()
wavelengths = absorbance[:,0]

# remove Wavelength column from both headers and absorbance data
headers = headers[1:]
absorbance = absorbance[:,1:]

#%%

# Plot raw data with legend
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance, headers, title = 'UV-Vis')

saverawimage = os.path.join(dirname, "Raw Data.jpg")
plt.savefig(saverawimage)


 ####################################### 
#%%
" BASELINE DATA AND PLOT "
# Create baseline subtracted data and plot the same graph

if len(baseline_wavelength) == 1:
    baseline_index = np.where(wavelengths == baseline_wavelength)[0][0]
    absorbance_BC = absorbance - absorbance[baseline_index]
    
elif len(baseline_wavelength) == 2:
    baseline_index = []
    for i in baseline_wavelength:
        baseline_index = np.append(baseline_index, np.where(wavelengths == i)[0][0])
    baseline_index = baseline_index.astype(int)
    baseline_index = np.sort(baseline_index)
    abs_correction = np.average(absorbance[baseline_index[0]:baseline_index[1]], axis=0)
    absorbance_BC = absorbance - abs_correction

# Plot baselined data with legend
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance_BC, headers, title = 'UV-Vis Baselined')

# Save Plot
saveimage1 = os.path.join(dirname, "Baselined Data.jpg")
plt.savefig(saveimage1)

# Plot baselined data with zoomed in ylims
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance_BC, headers, title = 'UV-Vis Baselined', ylims=[None,None])

# Save Plot
saveimage2 = os.path.join(dirname, "Baselined Data Zoom-In.jpg")
plt.savefig(saveimage2)


#######################################   
#%%
# Calculate CF Concentration from linear regression data

" Beer's Law: A = e*L*C "

CF_ext_coefficient = 74000	# M^-1 cm^-1
pathlength = 1 # cm

wv_450 = np.where(wavelengths == 450)[0][0]
wv_550 = np.where(wavelengths == 550)[0][0]
wvrange = np.sort([wv_450, wv_550])

# find max absorbances for baselined sample
peak_abs_range = absorbance_BC[wvrange[0]:wvrange[1]]
abs_max = np.max(peak_abs_range, axis=0)

if calculate_Stock_CF_concentration == True:
    
    # pull out dilution factors from excel shet
    DF = pd.read_excel(fname, sheet_name=Dilution_Factors_SheetName)
    DF = DF.to_numpy()
    DF = DF[0][1:]
    
    if len(np.unique(DF)) == 1:
        
        # if there's only one data point, average the max absorbances 
        # and compute concentration
        
        abs_max_avg = np.average(abs_max)
        conc_avg = abs_max_avg  / (CF_ext_coefficient * pathlength) # M
        
        conc_avg_mM = conc_avg * 10 ** 3 # mM
        conc_avg_uM = conc_avg * 10 ** 6 # nM
        conc_avg_nM = conc_avg * 10 ** 9 # nM
        
        fig, ax = plt.subplots()
        ax.scatter(DF, abs_max)
        ax.set_xlabel('Dilution Factor')
        ax.set_ylabel('Absorbance')
        ax.set_title('Calculating CF Concentration')
        ax.set_xlim(np.unique(DF)*[0.8,1.2])
        ax.set_ylim([0, np.max(abs_max)*1.2])
        plt.annotate("Stock Concentration = {} uM".format(np.round(conc_avg_uM, decimals = 5)), [np.unique(DF)*0.9, 0.8*np.max(abs_max)])


    else:    
        
        # convert absorbance maxes to concentrations
        conc_max_M = [x / (CF_ext_coefficient * pathlength) for x in abs_max] # Molar
        conc_max_mM = [x * 10**3 for x in conc_max_M]
        conc_max_uM = [x * 10**6 for x in conc_max_M]
       
        
        # scatter plot of maximum absorbance of baselined sample to dilution factor
        fig, ax = plt.subplots()
        ax.scatter(DF,abs_max)
        
        # linear regression
        slope = np.polyfit(DF.astype(float), abs_max.astype(float), 1)[0]
        intercept = np.polyfit(DF.astype(float), abs_max.astype(float), 1)[1]
        
        Xvals = np.linspace(0,np.max(DF))
        Yvals = Xvals * slope
        YvalsInt = Xvals * slope + intercept
        
        ax.set_xlabel('Dilution Factor')
        ax.set_ylabel('Absorbance')
        ax.set_title('Calculating CF Concentration')
        ax.plot(Xvals, Yvals, label= "y = {}x".format(np.round(slope, decimals=0)))
        ax.plot(Xvals, YvalsInt, label= "y = " + str(np.round(slope, decimals=0)) + "x + " + str(np.round(intercept, decimals=3)))
        
        # new fit with zero intercept
        x = np.vstack([DF, np.ones(len(DF))]).T
        a, _, _, _ = np.linalg.lstsq(x, abs_max)
        
        plt.plot(x, y, 'bo')
        plt.plot(x, a*x, 'r-')
        plt.show()
        
        
        
        ax.legend()
        
        # At dilution factor 1, absorbance is expected to be:
        abs_at_DF_1 = slope
        # Then, A = e*l*c, or, c = A/(e*l)
        conc_M = abs_at_DF_1 / (CF_ext_coefficient * pathlength) # Molar
        conc_mM = conc_M * 10**3 # mM
        
        plt.annotate("Stock Concentration = {} mM".format(np.round(conc_mM, decimals = 1)), [np.min(DF), 0.8*np.max(abs_max)])

#%% Calculate individual CF concentrations and print out table

if calculate_All_CF_concentrations == True:
    
    #convert max absorbance (after baseline) to concentration
    conc_M  = [x / (CF_ext_coefficient * pathlength) for x in abs_max] # Molar
    conc_mM = [x * 10**3 for x in conc_M] # mM
    conc_uM = [x * 10**6 for x in conc_M] # uM
    conc_nM = [x * 10**9 for x in conc_M] # uM
    
    # convert to dataframe
    results = pd.DataFrame([conc_M,conc_mM,conc_uM, conc_nM],
                           columns=headers, 
                           index=['M conc.', 'mM conc.', 'uM conc.', 'nM conc.'])
    results = results.T
    print(results)     
    
    # save table of concentration data in excel in same folder
    savefile = os.path.join(dirname,"results.xlsx")    
    results.to_excel(savefile)  
    
    
    
    



















