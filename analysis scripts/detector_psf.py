# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 13:44:41 2015

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import scipy.special as scs
import scipy.optimize as opti
import pyCT as ct
import os

sys.path.append('/users/schaff/Python Scripts/')
"""
Script determing the witdth of focal size and shape
with line-spread-function analysis
"""
plt.close("all")
#============================================================================================================
#                        Data input
#============================================================================================================
""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
""" Folder containing the data"""
data_folder = "detector_psf_60kvp"
"""   Vertical edges """
flatfield_v = 78122 # picturenumber of first flat of vertical direction
pic_folder = data+data_folder+"/paxscan/" # data folder 
sub_dirs_v=os.listdir(pic_folder)
data =list() ### the real pictures
flatfields =list()
data_vert = list()
data_hor = list()
dir_num=0
for sub_dir_v in sub_dirs_v:    ##### check the subfolders and sort all pictures in it in the data folder 
    dirs=np.sort(os.listdir(pic_folder+ sub_dir_v)) 
    for dir_name in dirs:                                               
        data.append(pic_folder+ sub_dir_v+'/'+dir_name)
#        else:    
#            flats_v.append(pic_folder+ sub_dir_v+'/'+dir_name)
        dir_num=dir_num+1

flatfields.append(data[0])
flatfields.append(data[2])
flatfields.append(data[4])
flatfields.append(data[6])
flatfields.append(data[8])
flatfields.append(data[10])

data_vert.append(data[1])
data_vert.append(data[3])
data_vert.append(data[5])
data_vert.append(data[7])
data_vert.append(data[9])
data_vert.append(data[11])

data_hor.append(data[12])
data_hor.append(data[13])
data_hor.append(data[14])
data_hor.append(data[15])
data_hor.append(data[16])
data_hor.append(data[17])

#img_vert = e17.io.h5read(data_vert[0])["raw_data"] #[430:630,220:420]
#1/0
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = data_folder+'.txt'
# write the header for txt file
info = 'contains at first column the different Watts, then the detector PSF'

#============================================================================================================
#                       Analysis options
#============================================================================================================
# determine the exact projection angle (set True) or you know it (set False)
proj_angle = True   ##### If you run the script the first time do this option it's very profitable ####
#choose your edge directions
fit_psf_v = True
fit_psf_h = True
# select your detector 
pilatus = False
paxscan = True
ccd = False
# Set True if you want plots
plots = True
# Set True if you want your results saved
saveon = False
#set true analyzing errorfit on data and gausspeak 
edge_fit = True
#============================================================================================================
#                       Function definitions
#============================================================================================================ 
#errorfunction fiting the lineedge
def error_fit_func(x,a,b,c,d): 
    return (a*scs.erf((x-b)/(math.sqrt(2)*c))+d)
#gaussfunction determin the FWHM
def gauss (x,a,b,c,d,):
    return (a*np.exp(-0.5*((x-b)/c)**2) + d)
#============================================================================================================
#                       Analysis parameter
#============================================================================================================

angle_tilt = 0.      
pic_number = 1 # Number of pictures analyzing per power
powersteps = 6
pixel_subdivision = 0.01
#Projector settings
if pilatus:
    detector_pixel_size = 172.
if paxscan:
    detector_pixel_size = 127.

#============================================================================================================
#                    Analysis for vertical edge
#============================================================================================================
if fit_psf_v:
    vertical_angle = angle_tilt - 90.
    P = ct.P.Parallel_AWP()
    P.angles = vertical_angle
    P.is_rescale = True #for subdivision
    
    mean_psf_hor = np.ones((powersteps))
    psf_hor = np.ones((powersteps,pic_number))
    fit_params_v = np.ones((powersteps,3,4))    
    for j in range(powersteps):
        if paxscan:
            #flatfield = e17.io.h5read(data+sam_v_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_v +j*12))["raw_data"]           
            flatfield = e17.io.h5read(flatfields[j])["raw_data"]
        fwhm_vert_pix = np.ones((pic_number))
        for i in np.arange(pic_number):                
            if paxscan:
                #img_vert = e17.io.h5read(data+sam_v_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_v +1+j*12+i))["raw_data"]
                img_vert = e17.io.h5read(data_vert[j])["raw_data"]
                #1/0
            img_vert = (img_vert/(flatfield*1.))[470:670,120:320]
            #1/0
            pic_shape = img_vert.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_vert = P.project(img_vert)[0]/np.sqrt(2)
            proj_norm_vert = np.nan_to_num(proj_vert/norm)
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = (detector_pixel_size * np.sqrt(2)) * pixel_subdivision
            x_val = np.arange(pic_shape*int(1/pixel_subdivision))*1.* projected_pixel_size
            #1/0
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;
            
            cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision)
            
            a = 1. # Amplitude                                                                           #######without max_x_val 
            b = 0.5# offset in x-direction                                                             ####### b = 20000
            c = 0.5#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge   ####### c = 337
            d = 0. # offset in y-direction
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff], proj_norm_vert[cutoff :-cutoff], p0=start_params)
            esf_fit = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and adapt fitparameter 
            
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            if edge_fit:
                plt.figure('edge + fit,Gauss vert_egde')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff],proj_norm_vert[cutoff :-cutoff])
                plt.plot(x_val[cutoff :-cutoff],esf_fit[cutoff :-cutoff])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2] ), best_fitparam[3]))
            fwhm_vert_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
            fit_params_v[j,i] = best_fitparam
        psf_hor[j,0] = fwhm_vert_pix[i]
        
    for x in range (powersteps):
        mean_psf_hor[x] = (psf_hor[0]+ psf_hor[1]+ psf_hor[2]+psf_hor[3]+psf_hor[4]+psf_hor[5])/6

#============================================================================================================
#                   Analysis for horizontal edge
#============================================================================================================
if fit_psf_h:
    horizontal_angle = angle_tilt
    P = ct.P.Parallel_AWP()
    P.angles = horizontal_angle
    P.is_rescale = True #for subdivision    
    
    mean_psf_vert = np.ones((powersteps))
    psf_vert = np.ones((powersteps,pic_number))
    fit_params_h = np.ones((powersteps,3,4))
    for j in range(powersteps):
        if paxscan:
            #flatfield = e17.io.h5read(data+sam_h_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_h +j*12))["raw_data"]           
            flatfield = e17.io.h5read(flatfields[j])["raw_data"]            
        fwhm_hor_pix = np.ones((pic_number))
        for i in np.arange(pic_number):               
            if paxscan:
                #img_hor = e17.io.h5read(data+sam_h_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_h +1+j*12+i))["raw_data"]
                img_hor = e17.io.h5read(data_hor[j])["raw_data"]
                #1/0
            img_hor = (img_hor/(flatfield*1.))[430:630,220:420]
            #1/0
            pic_shape = img_hor.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_hor = P.project(img_hor)[0]/np.sqrt(2)
            proj_norm_hor = np.nan_to_num(proj_hor/norm)
            #1/0
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision
            x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision))
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;
            
            cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision)
            
            a = 1. # Amplitude 
            b = 0.5 # offset in x-direction
            c = 0.5#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge
            d = 0. # offset in y-direction
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets 
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff ], proj_norm_hor[cutoff :-cutoff ], p0=start_params)
            esf_fit = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and amean_psf_vert = np.ones((powersteps))dapt fitparameter 
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            if edge_fit:
                plt.figure('edge + fit, Gauss hor_edge')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff ],proj_norm_hor[cutoff :-cutoff ])
                plt.plot(x_val[cutoff :-cutoff ],esf_fit[cutoff :-cutoff ])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2]), best_fitparam[3]))
    
            fwhm_hor_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])           
            fit_params_h[j,i] = best_fitparam
        psf_vert[j,0] = fwhm_hor_pix[0]
            
    for x in range (powersteps):
        mean_psf_vert[x] = (psf_vert[0]+ psf_vert[1]+ psf_vert[2]+psf_vert[3]+psf_vert[4]+psf_vert[5])/6
#============================================================================================================
#                   Plots and data output
#============================================================================================================
watts = np.asarray([5,10,20,50,100,150])
# bring data in right shape for output
if plots:
    if fit_psf_v:
        x_max_v = (1.1*np.max(watts))
        y_max_v = (1.4*np.max(psf_hor))
        v = [ 0., x_max_v,  200., y_max_v]
        plt.figure('Results HORIZONTAL PSF for 60kvp' )
        plt.axis(v)
        plt.plot(watts, psf_hor, lw=2 )
        plt.plot(watts, mean_psf_hor, lw=2, ls = '--', color = 'k')
        plt.title('Horizontal PSF 60kvp')
        plt.ylabel('PSf [micron]')
        plt.xlabel('Power [Watts]')
        plt.legend(('Detector PSF','Mean Detector PSf') ,loc = 0)
    if fit_psf_h:
        x_max_h = (1.1*np.max(watts))
        y_max_h = (1.4*np.max(psf_vert))
        z = [ 0., x_max_h,  200., y_max_h]
        plt.figure('Results VERTICAL PSF for 60kvp')
        plt.axis(z)
        plt.plot(watts, psf_vert, lw=2)
        plt.plot(watts, mean_psf_vert, lw=2, ls = '--', color = 'k')
        plt.title('Vertical PSF for 60kvp')
        plt.ylabel('PSF [micron]')
        plt.xlabel('Power [Watts]')
        plt.legend(  ('Detector PSF','Mean Detector PSf'),loc = 0)
        
if saveon:
    output_data = np.ones((3,powersteps))
    for y in range(powersteps):
        output_data[0,y] = watts[y]
        output_data[1,y] = psf_hor[y]
        output_data[2,y] = psf_vert[y]
    output_data = np.swapaxes(output_data,0,1)
    try:
        os.makedirs(filepath)    
    except:
        print('Folder already exist')
    try:
       # np.savetxt(filepath+filename,((watts,mean_spot_size_vertical,mean_spot_size_horizontal,mean_spot_size_diagonal)),fmt = '%5.5f',delimiter = ',',header = info)
        np.savetxt(filepath+filename, output_data,fmt = '%10.5f',delimiter = ',',header = info)
    except:
        print('Error, could not save the file (does the script/terminal have write access?)')









