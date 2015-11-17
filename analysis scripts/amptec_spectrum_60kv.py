# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 17:48:08 2015

@author: ga56pan
"""

import numpy as np
import scipy.constants as cons
import re
import os
from scipy import ndimage
import csv
import matplotlib.pyplot as plt
import ddfSetupProcessing.HarmonicAnalysis.lsqprocessing as dpcprocess


plots = False
corrig = True
kvp60 = False
kvp80 = True
#plt.close('all')
#============================================================================================================
#                    Readin data
#============================================================================================================
""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
data_folder = "spectraCDTE/"
if kvp60:
    file_id ="spectra_60kvp.csv"
if kvp80:    
    file_id ="spectra_80kvp.csv"
ff_file=open(data+data_folder+file_id)


ff_csv=csv.reader(ff_file)


ff_data=[]


for row in ff_csv:
    ff_data.append(row[0])      

#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath_p = '/users/Baier/analysis_pictures/spectra_amptec/'
filepath_d = '/users/Baier/analysis_files/'      
#============================================================================================================
#                    fuction definitions 
#============================================================================================================ 
cal_ch_a=float(ff_data[7][0:4]) ## channel for offset value of calibration
cal_en_a=float(ff_data[7][5:7]) ## energy for offset value of calibration
B = float(ff_data[10][10:19])    ## slope of the calibration 
A=cal_en_a-B*cal_ch_a
def energy_calibration(x,A,B):
    return(A+B*x)
#============================================================================================================
#                    analyis
#============================================================================================================
"""
correction factors linearabs coeff
"""
#if kvp80:
abscoeff_cdte_file = open(filepath_p+'lin_abs_coeff_cdte_80kvp.csv')
abscoeff_cdte = csv.reader(abscoeff_cdte_file)
total_ion_mu_cdte =  []
for row in abscoeff_cdte:
    total_ion_mu_cdte.append(row[1])
total_ion_mu_cdte.pop(0)
abscoeff_be_file = open(filepath_p+'lin_abs_coeff_be_80kvp.csv')
abscoeff_be = csv.reader(abscoeff_be_file)
total_ion_mu_be =  []
for row in abscoeff_be:
    total_ion_mu_be.append(row[1])
total_ion_mu_be.pop(0)
abscoeff_csi_file = open(filepath_p+'lin_abs_coeff_csi_80kvp.csv')
abscoeff_csi = csv.reader(abscoeff_csi_file)
total_ion_mu_csi =  []
for row in abscoeff_csi:
    total_ion_mu_csi.append(row[1])
total_ion_mu_csi.pop(0)
#if kvp60:
#    abscoeff_cdte_file = open(filepath_p+'lin_abs_coeff_cdte_60kvp.csv')
#    abscoeff_cdte = csv.reader(abscoeff_cdte_file)
#    total_ion_mu_cdte =  []
#    for row in abscoeff_cdte:
#        total_ion_mu_cdte.append(row[1])
#    total_ion_mu_cdte.pop(0)
#    abscoeff_be_file = open(filepath_p+'lin_abs_coeff_be_60kvp.csv')
#    abscoeff_be = csv.reader(abscoeff_be_file)
#    total_ion_mu_be =  []
#    for row in abscoeff_be:
#        total_ion_mu_be.append(row[1])
#    total_ion_mu_be.pop(0)
#    abscoeff_csi_file = open(filepath_p+'lin_abs_coeff_csi_60kvp.csv')
#    abscoeff_csi = csv.reader(abscoeff_csi_file)
#    total_ion_mu_csi =  []
#    for row in abscoeff_csi:
#        total_ion_mu_csi.append(row[1])
#    total_ion_mu_csi.pop(0)
###
#corresponding material thickness: important!!!! change unit to cm
cdte =  0.1
be   = 0.01
csi = 0.06
"""
"""
p_two=10e-6
d=0.925-0.04

ch_corr=np.empty((8206-14,1))

ff=np.empty((6,8206-14,1))
energies_raw=np.empty((8206-14,1))

for j in np.arange(0,8206-14,1): ### go here through the csv line per line and save the data 
    ff_line=map(int,re.findall(r'\d+',ff_data[j+14]))
    for i in np.arange(0,6,1):    
        ff[i,j,0]=float(ff_line[2*i+1])
        energies_raw[j] = energy_calibration(j,A,B)
        if j==0:
            ch_corr[j]=100000000
        else:    
            ch_corr[j]=(d*cons.h*cons.c)/(1000.*energies_raw[j]*cons.e*p_two)
            
smooth_ff_raw=np.empty((6,8206-14,1))
       
for i in np.arange(0,6,1):  
    smooth_ff_raw[i]=ndimage.uniform_filter(ff[i],size=50)
#coeff_ff = dpcprocess.lsq_fit(smooth_ff, nb_periods = 1, order = 1)
#ff = ff[:,59:6430,:]


if kvp80:         
    energies = energies_raw[59:6430]
    smooth_ff  =smooth_ff_raw[:,59:6430,:]
if kvp60:
    energies = energies_raw[59:4836]
    smooth_ff  =smooth_ff_raw[:,59:4836,:]
if corrig:
    corrig_ff = np.ones((6,np.shape(smooth_ff)[1]))
    corrig_ff_CSI = np.ones((6,np.shape(smooth_ff)[1]))
    for j in range(6):
        for i in range(np.shape(smooth_ff)[1]):
            #corrig_ff_div[j,i] = smooth_ff[j,i,:]*np.exp(-(float(total_ion_mu_cdte[i])*cdte+float(total_ion_mu_be[i])*be))
            corrig_ff_CSI[j,i] = (smooth_ff[j,i,:]/(1-np.exp(-float(total_ion_mu_cdte[i])*cdte)))*np.exp(float(total_ion_mu_be[i])*be)*(1-np.exp(-float(total_ion_mu_csi[i])*csi))
            corrig_ff[j,i] = (smooth_ff[j,i,:]/(1-np.exp(-float(total_ion_mu_cdte[i])*cdte)))*np.exp(float(total_ion_mu_be[i])*be)#*(1-np.exp(-float(total_ion_mu_csi[i])*csi))
#total_ion_mu_csi[59:4836]
#pp = np.ones(np.shape(smooth_ff)[1])            
#qq = np.ones(np.shape(smooth_ff)[1])           
#zz = np.ones(np.shape(smooth_ff)[1])
#for i in range(np.shape(smooth_ff)[1]):
#    zz[i] =(1-np.exp(-float(total_ion_mu_csi[i])*csi))
#    qq[i] =float(total_ion_mu_csi[i])*csi
#    pp[i] =-float(total_ion_mu_csi[i])*csi
if kvp60:
    plt.figure('Source spectra with without gratings')
    plt.plot(energies,smooth_ff[0],lw =.5)
    plt.plot(energies,smooth_ff[1],lw =.5)
    plt.plot(energies,smooth_ff[2],lw =.5)
    plt.plot(energies,smooth_ff[3],lw =.5)
    plt.plot(energies,smooth_ff[4],lw =.5)
    plt.plot(energies,smooth_ff[5],lw =.5)
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(smooth_ff)*1.05,500) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(smooth_ff)*1.05))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('G_0 +G_1 + G_2','G_2','G_1 + G_2','G_1','G_0','Spectrum without gratings') ,loc = 0,fontsize = 10,title = 'Grating combinations @ 60kVp')                
   # plt.savefig(filepath_p+'Spectra_at_60kVp_map.pdf', format = 'pdf',dpi = 300 )    
#plt.figure('g0/g2')
#plt.plot(energies,corrig_ff_CSI[4]/corrig_ff_CSI[1])    
    
    plt.figure('Spectra corrected and uncorrected at 60kvp')
    plt.plot(energies,corrig_ff_CSI[5]*2.,lw = 2, label = 'corrected')
    #plt.plot(energies,corrig_ff[5],lw = 2, label = 'corrected ')
    plt.plot(energies,smooth_ff[5]*2.,lw =2, label ='uncorrected')
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(corrig_ff)*2.*1.05,1000) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(corrig_ff_CSI)*2.*1.03))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend( loc = 0,fontsize = 10,title = 'Spectra @ 60kVp')                
   # plt.savefig(filepath_p+'Spectra_nograt_combi_corrig_and_CsI_at_60kVp_map.pdf', format = 'pdf',dpi = 300 )
    
    
    plt.figure('corrected spectra 60kvp')
    plt.plot(energies,corrig_ff_CSI[0])
    plt.plot(energies,corrig_ff_CSI[1])
    plt.plot(energies,corrig_ff_CSI[2])
    plt.plot(energies,corrig_ff_CSI[3])
    plt.plot(energies,corrig_ff_CSI[4])
    #plt.plot(energies,corrig_ff_CSI[5])
    plt.plot(energies,corrig_ff_CSI[5]*2.)## account half power
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.02)
    plt.yticks( np.arange(0, np.max(corrig_ff_CSI)*2.*1.05,1000) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(corrig_ff_CSI)*2.*1.05))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('G_0 +G_1 + G_2','G_2','G_1 + G_2','G_1','G_0','Spectrum without gratings') ,loc = 0,fontsize = 10,title = 'Grating combinations @ 60kVp')                
   # plt.savefig(filepath_p+'Spectra_corrig_and_CsI_at_60kVp_map.pdf', format = 'pdf',dpi = 300 )
    
    
if kvp80:
    plt.figure('Source spectra with without gratings')
    plt.plot(energies,smooth_ff[0],lw =.5)
    plt.plot(energies,smooth_ff[1],lw =.5)
    plt.plot(energies,smooth_ff[2],lw =.5)
    plt.plot(energies,smooth_ff[3],lw =.5)
    plt.plot(energies,smooth_ff[4],lw =.5)
    plt.plot(energies,smooth_ff[5],lw =.5)
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(smooth_ff)+1,500) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(smooth_ff)*1.01))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('G_0 +G_1 + G_2','G_2','G_1 + G_2','G_1','G_0','Spectrum without gratings') ,loc = 0,fontsize = 10,title = 'Grating combinations @ 80kVp')                
   # plt.savefig(filepath_p+'Spectra_at_80kVp_map.pdf', format = 'pdf',dpi = 300 )
    
    
    plt.figure('Spectra corrected and uncorrected at 80kvp')
    plt.plot(energies,corrig_ff_CSI[5]*2.,lw = 2, label = 'corrected')
    #plt.plot(energies,corrig_ff[5],lw = 2, label = 'corrected')
    plt.plot(energies,smooth_ff[5]*2.,lw =2, label ='uncorrected')
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(corrig_ff)*2.*1.05,1000) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(smooth_ff)*2.*1.01))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend( loc = 0,fontsize = 10,title = 'Spectra @ 80kVp')                
    plt.savefig(filepath_p+'Spectra_norgrat_combi_corrig_and_CsI_at_80kVp_map.pdf', format = 'pdf',dpi = 300 )
    
       
    plt.figure('corrected spectra at 80kvp')
    plt.plot(energies,corrig_ff_CSI[0])
    plt.plot(energies,corrig_ff_CSI[1])
    plt.plot(energies,corrig_ff_CSI[2])
    plt.plot(energies,corrig_ff_CSI[3])
    plt.plot(energies,corrig_ff_CSI[4])
    #plt.plot(energies,corrig_ff_CSI[5])
    plt.plot(energies,corrig_ff_CSI[5]*2.)###account half power
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(corrig_ff_CSI)*2.*1.05,1000) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(corrig_ff_CSI)*2.*1.05))
    plt.ylabel('Intensity $[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('G_0 +G_1 + G_2','G_2','G_1 + G_2','G_1','G_0','Spectrum without gratings') ,loc = 0,fontsize = 10,title = 'Grating combinations @ 80kVp')                
    plt.savefig(filepath_p+'Spectra_corrig_and_CsI_at_80kVp_map.pdf', format = 'pdf',dpi = 300 )




