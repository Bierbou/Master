# -*- coding: utf-8 -*-
"""
Created on Thu May 21 13:24:51 2015

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
kv_40 = False
kv_80 = True
""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
""" Folder containing the data"""
"""   Vertical edges """
if kv_40:
#    watts = 60
#    data_folder = "detector_psf_40kvp"
#    flatfield_pic = data+data_folder+"/paxscan/ct/paximage_ct_183796.h5" # data folder
#    vert_edge = data+data_folder+"/paxscan/ct/paximage_ct_183797.h5"
#    hor_edge = data+data_folder+"/paxscan/ct/paximage_ct_183798.h5"
    watts = 50
    data_folder = "detector_psf_60kvp_control/vert"
    flatfield_pic_h = data+data_folder+"/paxscan/ct/paximage_ct_763694.h5" # data folder
    flatfield_pic_v = data+data_folder+"/paxscan/ct/paximage_ct_763696.h5" # data folder
    vert_edge = data+data_folder+"/paxscan/ct/paximage_ct_763697.h5"
    hor_edge = data+data_folder+"/paxscan/ct/paximage_ct_763695.h5"
if kv_80:
#    watts = 100
#    data_folder = "detector_psf_80kvp"
#    flatfield_pic = data+data_folder+"/paxscan/ct/paximage_ct_183799.h5" # data folder
#    vert_edge = data+data_folder+"/paxscan/ct/paximage_ct_183800.h5"
#    hor_edge = data+data_folder+"/paxscan/ct/paximage_ct_183801.h5"
    watts = 100
    data_folder = "detector_psf_60kvp_control/vert"
    flatfield_pic_h = data+data_folder+"/paxscan/ct/paximage_ct_763692.h5" # data folder
    flatfield_pic_v = data+data_folder+"/paxscan/ct/paximage_ct_763698.h5" # data folder
    vert_edge = data+data_folder+"/paxscan/ct/paximage_ct_763699.h5"
    hor_edge = data+data_folder+"/paxscan/ct/paximage_ct_763693.h5"
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath_pic = '/users/Baier/analysis_pictures/'
filepath_data = '/users/Baier/analysis_files/'
#filename = data_folder+'.txt'
if kv_40:
    filename = 'detector_psf_40kvp_control.txt'
if kv_80:
    filename = 'detector_psf_80kvp_control.txt'
# write the header for txt file
info = 'contains at first column the Watt value, then the detector PSF in horizontal direction, then in vertical direction'
#============================================================================================================
#                       Analysis options
#============================================================================================================
proj_angle  =False
#choose your edge directions
fit_psf_v = True
fit_psf_h = True
# select your detector 
paxscan = True
# Set True if you want plots
# Set True if you want your results saved
saveon = True
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
a = 1. # Amplitude                                                                           #######without max_x_val 
b = 0.5# offset in x-direction                                                             ####### b = 20000
c = 0.5#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge   ####### c = 337
d = 0. # offset in y-direction
angle_tilt_v = -0.7 
angle_tilt_h = 1.7   
 # Number of pictures analyzing per power
pixel_subdivision = 0.01
#Projector settings
if paxscan:
    detector_pixel_size = 127.
#============================================================================================================
#                   find best projection angle for vertical edge
#============================================================================================================
if proj_angle:
    start_angle = angle_tilt_v - .5
    end_angle = angle_tilt_v + .5
    proj_number = 1+(end_angle -start_angle)/0.1

 ### take pictures with higher power -> bigger spot -> easier to find the right angle. Choose 80  watts
    flatfield = e17.io.h5read(flatfield_pic_v)["raw_data"]
    img_vert = e17.io.h5read(vert_edge)["raw_data"]
    #1/0
    img_vert = (img_vert/(flatfield*1.))[380:550,155:325]
    pic_shape = img_vert.shape[0]
    fwhm = np.zeros(proj_number)
    angle = np.zeros(proj_number)
    for i in np.arange(proj_number): # loop over all angles in 0.1° steps 
        # projection of data and normation
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = start_angle + 90. # projetion angle
        P.is_rescale = True #for subdivision      
        P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2) # project ones as normation
        proj_vert = P.project(img_vert)[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm = np.nan_to_num(proj_vert/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data 
        max_x_val = max(x_val)   # normation to the max value makes it easier converging the fit all the time 
        x_val = x_val/max_x_val;
        
        cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision) # cause of projecting a square data set 
        start_params = [a,b,c,d] # Startparameter for errorfit                  # left and right zeros are added get rid of them    
        # fitting datasets 
        best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff:-cutoff], proj_norm[cutoff:-cutoff], p0=start_params)
        fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
        start_angle += 0.1
        angle[i] = start_angle
        print fwhm[i]
    plt.figure('fwhm_min_v')
    plt.plot(angle,fwhm)
    best_angle_vert = angle[np.argmin(fwhm)]
else:
    best_angle = angle_tilt_v
#============================================================================================================
#                   find best projection angle for horizontal edge
#============================================================================================================
if proj_angle:
    start_angle = angle_tilt_h - .5
    end_angle = angle_tilt_h + .5
    proj_number = 1+(end_angle -start_angle)/0.1

 ### take pictures with higher power -> bigger spot -> easier to find the right angle. Choose 80  watts
    flatfield = e17.io.h5read(flatfield_pic_h)["raw_data"]
    img_vert = e17.io.h5read(hor_edge)["raw_data"]
    #1/0
    img_vert = (img_vert/(flatfield*1.))[350:520,220:390]
    pic_shape = img_vert.shape[0]
    fwhm = np.zeros(proj_number)
    angle = np.zeros(proj_number)
    for i in np.arange(proj_number): # loop over all angles in 0.1° steps 
        # projection of data and normation
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = start_angle #+ 90. # projetion angle
        P.is_rescale = True #for subdivision      
        P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2) # project ones as normation
        proj_vert = P.project(img_vert)[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm = np.nan_to_num(proj_vert/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data 
        max_x_val = max(x_val)   # normation to the max value makes it easier converging the fit all the time 
        x_val = x_val/max_x_val;
        
        cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision) # cause of projecting a square data set 
        start_params = [a,b,c,d] # Startparameter for errorfit                  # left and right zeros are added get rid of them    
        # fitting datasets 
        best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff:-cutoff], proj_norm[cutoff:-cutoff], p0=start_params)
        fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
        start_angle += 0.1
        angle[i] = start_angle
        print fwhm[i]
    plt.figure('fwhm_min_h')
    plt.plot(angle,fwhm)
    best_angle_hor = angle[np.argmin(fwhm)]
else:
    best_angle = angle_tilt_h
#============================================================================================================
#                    Analysis for vertical edge
#============================================================================================================
if fit_psf_v:
    vertical_angle = angle_tilt_v - 90.
    P = ct.P.Parallel_AWP()
    P.angles = vertical_angle
    P.is_rescale = True #for subdivision 
    if paxscan:           
        flatfield = e17.io.h5read(flatfield_pic_v)["raw_data"]                
        img_vert = e17.io.h5read(vert_edge)["raw_data"]
        #1/0
    #img_vert = (img_vert/(flatfield*1.))[470:670,120:320]### for the ususal masurement 
    img_vert = (img_vert/(flatfield*1.))[380:550,155:325]### for the control measurent
    pic_shape = img_vert.shape[0]
    P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
    # projection of data and normation
    norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
    proj_vert = P.project(img_vert)[0]/np.sqrt(2)
    proj_norm_vert = np.nan_to_num(proj_vert/norm)
    # generation of x-values for plot and normation for fitting
    projected_pixel_size = (detector_pixel_size * np.sqrt(2)) * pixel_subdivision
    x_val_v = np.arange(pic_shape*int(1/pixel_subdivision))*1.* projected_pixel_size
    max_x_val_v = max(x_val_v)
    x_val_v = x_val_v/max_x_val_v;
    
    cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/(pixel_subdivision))
    
    a = 1. # Amplitude                                                                           #######without max_x_val 
    b = 0.5# offset in x-direction                                                             ####### b = 20000
    c = 0.5#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge   ####### c = 337
    d = 0. # offset in y-direction
    start_params = [a,b,c,d] # Startparameter for errorfit
    # fitting datasets
    best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val_v[cutoff :-cutoff], proj_norm_vert[cutoff :-cutoff], p0=start_params)
    esf_fit_vert = error_fit_func(x_val_v, *best_fitparam)
    v = cov_esf
    # convert x-values back to real and adapt fitparameter 
    #1/0
    x_val_v = x_val_v*max_x_val_v
    best_fitparam *= max_x_val_v
    psf_hor = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
    fit_params_v = best_fitparam

#============================================================================================================
#                   Analysis for horizontal edge
#============================================================================================================
if fit_psf_h:
    horizontal_angle = angle_tilt_h - 180.
    P = ct.P.Parallel_AWP()
    P.angles = horizontal_angle
    P.is_rescale = True #for subdivision    

    if paxscan:          
        flatfield = e17.io.h5read(flatfield_pic_h)["raw_data"]            
        img_hor = e17.io.h5read(hor_edge)["raw_data"]
        #1/0
    img_hor = (img_hor/(flatfield*1.))[350:520,220:390]
    pic_shape = img_hor.shape[0]
    P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
    # projection of data and normation
    norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
    proj_hor = P.project(img_hor)[0]/np.sqrt(2)
    proj_norm_hor = np.nan_to_num(proj_hor/norm)
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
    esf_fit_hor = error_fit_func(x_val[cutoff :-cutoff ], *best_fitparam)
    h = cov_esf
    # convert x-values back to real and amean_psf_vert = np.ones((powersteps))dapt fitparameter 
    x_val = x_val*max_x_val
    best_fitparam *= max_x_val
    psf_vert = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])           
    fit_params_h = best_fitparam

#============================================================================================================
#                   Plots and data output
#============================================================================================================
if edge_fit:
    plt.figure('edge + fit,Gauss vert_egde')
    plt.subplot(211)
    plt.plot(x_val[cutoff :-cutoff],proj_norm_vert[cutoff :-cutoff])
    plt.plot(x_val[cutoff :-cutoff],esf_fit_vert[cutoff :-cutoff])
    plt.subplot(212)
    plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2] ), best_fitparam[3]))
#    if kv_40:           
#        plt.savefig(filepath_pic+'Edge_fit and Gauss vert_edge_40kv.tiff', format = 'tiff')
#    if kv_80:           
#        plt.savefig(filepath_pic+'Edge_fit and Gauss vert_edge_80kv.tiff', format = 'tiff')

    plt.figure('edge + fit, Gauss hor_edge')
    plt.subplot(211)
    plt.plot(x_val[cutoff :-cutoff ],proj_norm_hor[cutoff :-cutoff ])
    plt.plot(x_val[cutoff :-cutoff ],esf_fit_hor)
    plt.subplot(212)
    plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2]), best_fitparam[3]))
#    if kv_40:           
#        plt.savefig(filepath_pic+'Edge_fit and Gauss hor_edge_40kv.tiff', format = 'tiff', dpi = 300)
#    if kv_80:           
#        plt.savefig(filepath_pic+'Edge_fit and Gauss hor_edge_80kv.tiff', format = 'tiff', dpi = 300)
        
# bring data in right shape for output
if saveon:
    output_data = np.ones((3))
    output_data[0] = watts
    output_data[1] = psf_hor
    output_data[2] = psf_vert
    output_data = np.swapaxes(output_data,1,0)
    try:
        os.makedirs(filepath_data)    
    except:
        print('Folder already exist')
    try:
       # np.savetxt(filepath+filename,((watts,mean_spot_size_vertical,mean_spot_size_horizontal,mean_spot_size_diagonal)),fmt = '%5.5f',delimiter = ',',header = info)
        np.savetxt(filepath_data+filename, output_data,fmt = '%10.5f',delimiter = ',',header = info)
    except:
        print('Error, could not save the file (does the script/terminal have write access?)')