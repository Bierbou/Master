# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 15:57:01 2015

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


corrig = False
plt.close('all')
#============================================================================================================
#                    Readin data
#============================================================================================================
""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
data_folder = "spectraCDTE/"
file_id_flatfield ='visimap comp.csv'
file_id_sample = 'SIO2_2.79shakecomp.csv'
ff_file=open(data+data_folder+file_id_flatfield)
im_file=open(data+data_folder+file_id_sample)

ff_csv=csv.reader(ff_file)
im_csv=csv.reader(im_file)

ff_data=[]
im_data=[]

for row in ff_csv:
    ff_data.append(row[0])      
for row in im_csv:
    im_data.append(row[0])
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
cal_en_a=float(ff_data[7][5:9]) ## energy for offset value of calibration
B = float(ff_data[10][8:17])    ## slope of the calibration 
A=cal_en_a-B*cal_ch_a
def energy_calibration(x,A,B):
    return(A+B*x)
#============================================================================================================
#                    analyis
#============================================================================================================
"""
correction factors linearabs coeff
"""
 
abscoeff_cdte_file = open(filepath_p+'lin_abs_coeff_cdte_100kvp.csv')
abscoeff_cdte = csv.reader(abscoeff_cdte_file)
total_ion_mu_cdte =  []
for row in abscoeff_cdte:
    total_ion_mu_cdte.append(row[1])
total_ion_mu_cdte.pop(0)
abscoeff_be_file = open(filepath_p+'lin_abs_coeff_be_100kvp.csv')
abscoeff_be = csv.reader(abscoeff_be_file)
total_ion_mu_be =  []
for row in abscoeff_be:
    total_ion_mu_be.append(row[1])
total_ion_mu_be.pop(0)
abscoeff_csi_file = open(filepath_p+'lin_abs_coeff_csi_100kvp.csv')
abscoeff_csi = csv.reader(abscoeff_csi_file)
total_ion_mu_csi =  []
for row in abscoeff_csi:
    total_ion_mu_csi.append(row[1])
total_ion_mu_csi.pop(0)

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

ff=np.empty((8,8206-14,1))
im=np.empty((8,8206-14,1))
energies_raw=np.empty((8206-14,1))

for j in np.arange(0,8206-14,1): ### go here through the csv line per line and save the data 
    ff_line=map(int,re.findall(r'\d+',ff_data[j+14]))
    im_line=map(int,re.findall(r'\d+',im_data[j+14]))
    for i in np.arange(0,8,1):    
        ff[i,j,0]=float(ff_line[2*i+1])
        im[i,j,0]=float(im_line[2*i+1])
        energies_raw[j] = energy_calibration(j,A,B)
        if j==0:
            ch_corr[j]=100000000
        else:    
            ch_corr[j]=(d*cons.h*cons.c)/(1000.*energies_raw[j]*cons.e*p_two)
            
smooth_ff_raw=np.empty((8,8206-14,1))
smooth_im_raw=np.empty((8,8206-14,1))

energies = energies_raw[339:7857]
smooth_ff = smooth_ff_raw[:,339:7857]
smooth_im = smooth_im_raw[:,339:7857]      
for i in np.arange(0,8,1):  
    smooth_ff_raw[i]=ndimage.uniform_filter(ff[i],size=50)
    smooth_im_raw[i]=ndimage.uniform_filter(im[i],size=50)

if corrig:
    corrig_ff_CSI = np.ones((8,np.shape(smooth_ff)[1]))
    corrig_im_CSI = np.ones((8,np.shape(smooth_ff)[1]))
    for j in range(8):
        for i in range(np.shape(smooth_ff)[1]):
            #corrig_ff_div[j,i] = smooth_ff[j,i,:]*np.exp(-(float(total_ion_mu_cdte[i])*cdte+float(total_ion_mu_be[i])*be))
            corrig_ff_CSI[j,i] = (smooth_ff[j,i,:]/(1-np.exp(-float(total_ion_mu_cdte[i])*cdte)))*np.exp(float(total_ion_mu_be[i])*be)*(1-np.exp(-float(total_ion_mu_csi[i])*csi))
            corrig_im_CSI[j,i] = (smooth_im[j,i,:]/(1-np.exp(-float(total_ion_mu_cdte[i])*cdte)))*np.exp(float(total_ion_mu_be[i])*be)*(1-np.exp(-float(total_ion_mu_csi[i])*csi))

    coeff_ff_corrig = dpcprocess.lsq_fit(corrig_ff_CSI, nb_periods = 1, order = 1)        
    coeff_im_corrig = dpcprocess.lsq_fit(corrig_im_CSI, nb_periods = 1, order = 1)
    
    a0_corrig = coeff_ff_corrig[0,0]
    a1_corrig = coeff_ff_corrig[0,1]
    visibility_corrig = a1_corrig/a0_corrig
    visibility_corrig[np.isnan(visibility_corrig)]=0.
    a0_im_corrig = coeff_im_corrig[0,0]
    a1_im_corrig = coeff_im_corrig[0,1]
    visibility_im_corrig = a1_im_corrig/a0_im_corrig
    visibility_im_corrig[np.isnan(visibility_im_corrig)]=0. 

    plt.figure(file_id_flatfield+'corrig')
    plt.plot(energies,visibility_corrig)
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)    
    plt.yticks( np.arange(0, np.max(visibility_corrig)*1.01,.05) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(visibility_corrig)*1.01))
    plt.ylabel('Visibility $[\%]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('Energy dependent visibility',) ,loc = 0,fontsize = 12)
    #plt.savefig(filepath_p+'energy_visi_map.pdf', format = 'pdf',dpi = 300 )
    
    plt.figure('stepping spectra corrig')
    plt.plot(energies,corrig_ff_CSI[0],lw =.5)
    plt.plot(energies,corrig_ff_CSI[1],lw =.5)
    plt.plot(energies,corrig_ff_CSI[2],lw =.5)
    plt.plot(energies,corrig_ff_CSI[3],lw =.5)
    plt.plot(energies,corrig_ff_CSI[4],lw =.5)
    plt.plot(energies,corrig_ff_CSI[5],lw =.5)
    plt.plot(energies,corrig_ff_CSI[6],lw =.5)
    plt.plot(energies,corrig_ff_CSI[7],lw =.5)
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    #plt.yticks( np.arange(0, np.max(corrig_ff_CSI)*1.01,500) )
    plt.yticks(fontsize  = 12)
    #plt.ylim((0,np.max(corrig_ff_CSI)*1.01))
    plt.ylabel('Intensity$[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =12)
    plt.legend(('Step 1','Step 2','Step 3','Step 4','Step 5','Step 5','Step 6','Step 7','Step 8') ,loc = 0,fontsize = 12, title = 'Stepping steps')                
    #plt.savefig(filepath_p+'energy_stepping_map.pdf', format = 'pdf',dpi = 300 )           

else:
    coeff_ff = dpcprocess.lsq_fit(smooth_ff, nb_periods = 1, order = 1)        
    coeff_im = dpcprocess.lsq_fit(smooth_im, nb_periods = 1, order = 1) 

    a0 = coeff_ff[0,0]
    a1 = coeff_ff[0,1]
    visibility = a1/a0
    visibility[np.isnan(visibility)]=0.
    a0_im = coeff_im[0,0]
    a1_im = coeff_im[0,1]
    visibility_im = a1_im/a0_im
    visibility_im[np.isnan(visibility_im)]=0.            
#atc, dpc, dfc = dpcprocess.dpcprocess(smooth_ff, smooth_im, nb_periods = 1, order = 1)
#plt.figure()
#plt.plot(energies,coeff_ff[0,1]/coeff_ff[0,0])        
#plt.figure()
#plt.plot(ch_corr,dfc)            
#
#try:
#    os.makedirs(filepath_p)    
#except:
#    print('Folder already exists')


    plt.figure(file_id_flatfield)
    plt.plot(energies,visibility)
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)    
    plt.yticks( np.arange(0, np.max(visibility)*1.01,.05) )
    plt.yticks(fontsize  = 10)
    plt.ylim((0,np.max(visibility)*1.01))
    plt.ylabel('Visibility $[\%]$')
    plt.xlabel('Energy $[keV]$',fontsize =10)
    plt.legend(('Energy resolved visibility',) ,loc = 0,fontsize = 12)
    #plt.savefig(filepath_p+'energy_visi_map.pdf', format = 'pdf',dpi = 300 )
    
    
    plt.figure('stepping spectra')
    plt.plot(energies,smooth_ff[0],lw =.8,label = 'Step 1',color = ((0,0,1)))
    plt.plot(energies,smooth_ff[1],lw =.8,label = 'Step 2',color = ((1,0,0)))
    plt.plot(energies,smooth_ff[2],lw =.8,label = 'Step 3',color = ((0,0,0)))
    plt.plot(energies,smooth_ff[3],lw =.8,label = 'Step 4',color = ((0,1,1)))
    plt.plot(energies,smooth_ff[4],lw =.8,label = 'Step 5',color = ((1,0,1)))
    plt.plot(energies,smooth_ff[5],lw =.8,label = 'Step 6',color = ((1,0.6,0)))
    plt.plot(energies,smooth_ff[6],lw =.8,label = 'Step 7',color = ((0,1,0)))
    plt.plot(energies,smooth_ff[7],lw =.8,label = 'Step 8',color = ((0.5,0,0.5)))
    plt.xticks( np.arange(0, np.max(energies)+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(energies)*1.01)
    plt.yticks( np.arange(0, np.max(smooth_ff)*1.05,500) )
    plt.yticks(fontsize  = 10)
    plt.ylim((0,np.max(smooth_ff)*1.04))
    plt.ylabel('Intensity$[arb.unit]$')
    plt.xlabel('Energy $[keV]$',fontsize =10)
    plt.legend(loc = 0,fontsize = 12, title = 'Stepping steps')                
    plt.set_cmap('rainbow')
    #plt.savefig(filepath_p+'energy_stepping_map.pdf', format = 'pdf',dpi = 300 )
            