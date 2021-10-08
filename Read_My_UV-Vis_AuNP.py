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

fname = r"G:\My Drive\Research\Landry Lab Summer Research 2021\AL Data\B2P84\UV-Vis 2021-10-08\2021-10-08_AuNP_supernatant.xlsx"

baseline_wavelength = [750, 800] #nm

calculate_AuNP_concentration = True
AuNP_size_in_nm              = 10 # 5 or 10
Dilution_Factors_SheetName = 'Sheet2' # must be in same file as fname

#######################################

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
 
#######################################  

# AuNP OD 1 constants (from nanoComposix spec sheet)
OD1_to_nM = {'5': 80.9800664451827, '10': 8.421052631578947}


#######################################  

fig,ax = plt.subplots()

def makeplot(xdata, ydata, headers, title='', AX=ax, xlims=[400,1000], ylims=[-0.1,1]):
    
    if len(np.shape(ydata)) == 1: # There is only one column to plot 
        AX.plot(xdata, ydata, label=headers[0])
    else:        
        for i in np.arange(0,np.shape(ydata)[1]):
            AX.plot(xdata, ydata[:,i], label=headers[i])
        
    AX.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    AX.set_xlabel('Wavelength (nm)')
    AX.set_ylabel('Absorbance')
    AX.set_title(title)
    
    plt.xlim(xlims)
    plt.ylim(ylims)

#######################################   
#%%
" FORMAT DATA "

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
" RAW DATA PLOT"

# Plot raw data with legend
fig, ax = plt.subplots()
makeplot(wavelengths, absorbance, headers, title = 'UV-Vis')


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
fig2, ax2 = plt.subplots()
makeplot(wavelengths, absorbance_BC, headers, title = 'UV-Vis Baselined', AX=ax2)


#######################################   
#%%
# Calculate AuNP Concentration from absorbance

if calculate_AuNP_concentration == True:
    
    " Beer's Law: A = e*L*C "
    
    OD_conversion = OD1_to_nM[str(AuNP_size_in_nm)]
    
    # look for max between 450 and 600
    
    wv_450 = np.where(wavelengths == 450)[0][0]
    wv_600 = np.where(wavelengths == 600)[0][0]
    
    # find max absorbances for baselined sample
    abs_max = np.max(absorbance_BC[wv_600:wv_450], axis=0)
    # find indices where these maxes are found
    max_index = []
    wv_for_max_abs = []
    for i in np.arange(0,len(abs_max)):
        max_index = int(np.where(absorbance_BC[:,i] == abs_max[i])[0][0])
        wv_for_max_abs = np.append(wv_for_max_abs, wavelengths[max_index])
    
    # pull out dilution factors from excel shet
    DF = pd.read_excel(fname, sheet_name=Dilution_Factors_SheetName)
    DF = DF.to_numpy()
    DF = DF[0][1:]
    
    if len(np.unique(DF)) == 1:
        
        # if there's only one data point, average the max absorbances 
        # and compute concentration
        
        abs_max_avg = np.average(abs_max)
        conc_avg = abs_max_avg * OD_conversion  # nM
                
        fig, ax = plt.subplots()
        ax.scatter(DF, abs_max)
        ax.set_xlabel('Dilution Factor')
        ax.set_ylabel('Absorbance')
        ax.set_title('Calculating AuNP Concentration')
        ax.set_xlim(np.unique(DF)*[0.8,1.2])
        ax.set_ylim([0, np.max(abs_max)*1.2])
        plt.annotate("Stock Concentration = {} nM".format(np.round(conc_avg, decimals = 5)), [np.unique(DF)*0.9, 0.8*np.max(abs_max)])


    # Also, plot the baseline subtracted chart again with arrows to concentrations 
        
    # convert max absorbances to concentrations
    conc = abs_max * OD_conversion #nM
    
    
    # add annotation to baselined plot
    for i in np.arange(0,len(headers)):
        ax2.annotate('AuNP Conc = {} nM'.format(np.round(conc[i], decimals = 5)), xy=(wv_for_max_abs[i], abs_max[i]),  xycoords='data',
            xytext=(wv_for_max_abs[i]+400, abs_max[i]),
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='top',
            )
    

#%%

# just for fun, plot 4 subplot



fig, ax4 = plt.subplots(round(len(headers)/2),2)
ax4 = ax4.flatten()
fig.set_size_inches(8, 6)

H,W = np.shape(absorbance_BC)
for i in np.arange(0,W):
    
    ax4[i].plot(wavelengths, absorbance_BC[:,i])
    ax4[i].set_title(headers[i])
    
    ax4[i].annotate('AuNP Conc = {} nM'.format(np.round(conc[i], decimals = 5)), xy=(wv_for_max_abs[i], abs_max[i]),  xycoords='data',
            xytext=(wv_for_max_abs[i]+400, abs_max[i]),
            arrowprops=dict(facecolor='black', shrink=0.05),
            horizontalalignment='right', verticalalignment='top',
            )

plt.tight_layout()
























