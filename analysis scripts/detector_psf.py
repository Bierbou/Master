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
#data_folder = "detector_psf_60kvp"
data_folder = "detector_psf_60kvp_control"
"""   Vertical edges """
flatfield_v = 763656 # picturenumber of first flat of vertical direction
#pic_folder = data+data_folder+"/paxscan/" # data folder 
pic_folder = data+data_folder+"/vert/paxscan/" # data folderfor 60kvp control
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

flatfields.append(data[0]),flatfields.append(data[2])##for 60kvp control
flatfields.append(data[4]),flatfields.append(data[6])
flatfields.append(data[8]),flatfields.append(data[10])
flatfields.append(data[12]),flatfields.append(data[14])
flatfields.append(data[16]),flatfields.append(data[18])
flatfields.append(data[20]),flatfields.append(data[22])
flatfields.append(data[24]),flatfields.append(data[26])
flatfields.append(data[28]),flatfields.append(data[30])
flatfields.append(data[32]),flatfields.append(data[34])

data_vert.append(data[1])##for 60kvp control
data_vert.append(data[3])
data_vert.append(data[5])
data_vert.append(data[7])
data_vert.append(data[9])
data_vert.append(data[11])
data_vert.append(data[13])
data_vert.append(data[15])
data_vert.append(data[17])

data_hor.append(data[19])##for 60kvp control
data_hor.append(data[21])
data_hor.append(data[23])
data_hor.append(data[25])
data_hor.append(data[27])
data_hor.append(data[29])
data_hor.append(data[31])
data_hor.append(data[33])
data_hor.append(data[35])

#flatfields.append(data[0])
#flatfields.append(data[2])
#flatfields.append(data[4])
#flatfields.append(data[6])
#flatfields.append(data[8])
#flatfields.append(data[10])
#
#data_vert.append(data[1])
#data_vert.append(data[3])
#data_vert.append(data[5])
#data_vert.append(data[7])
#data_vert.append(data[9])
#data_vert.append(data[11])
#
#data_hor.append(data[12])
#data_hor.append(data[13])
#data_hor.append(data[14])
#data_hor.append(data[15])
#data_hor.append(data[16])
#data_hor.append(data[17])
#


#img_vert = e17.io.h5read(data_vert[0])["raw_data"] #[430:630,220:420]
#1/0
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath_pic = '/users/Baier/analysis_pictures/detector_60kvp_control/'
filepath_data = '/users/Baier/analysis_files/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = data_folder+'nodrop.txt'
filenamemean = data_folder+'nodrop_meanpsfs.txt'
# write the header for txt file
info1 = 'contains at first column the different Watts, then the detector PSF'
info2 = 'contains at first column the mean_vert detectot PSF, then the mean_hor detector PSF'
#============================================================================================================
#                       Analysis options
#============================================================================================================
# find right projection angle
proj_angle_v = False
proj_angle_h = False
#choose your edge directions
fit_psf_v = True
fit_psf_h = True
# Set True if you want plots
plots = True
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
angle_tilt_vert = 1.## for 60 kvp contol
angle_tilt_hor = 1.6  ###for 60kvp control
#angle_tilt_vert = 0.
#angle_tilt_hor = 0. 
pic_number = 1 # Number of pictures analyzing per power
#powersteps = 6
powersteps = 9#for 60 kvp control
pixel_subdivision = 0.01
## startparameter for the errrorfit
a = 1. # Amplitude                 #######without max_x_val 
b = 0.5# offset in x-direction      ####### b = 20000
c = 0.5# width of the edge          ####### c = 337
d = 0. # offset in y-direction

#Projector settings

detector_pixel_size = 127.
#============================================================================================================
#                   find best projection angle vertical edge 
#============================================================================================================
if proj_angle_v:
    start_angle = angle_tilt_vert - .5
    end_angle = angle_tilt_vert + .5
    proj_number = 1+(end_angle -start_angle)/.1

 ### take pictures with higher power -> bigger spot -> easier to find the right angle. Choose 80  watts
    flatfield = e17.io.h5read(flatfields[6])["raw_data"]
    img_vert = e17.io.h5read(data_vert[6])["raw_data"]
    img_vert = (img_vert/(flatfield*1.))[370:560,135:325]
    pic_shape = img_vert.shape[0]
    fwhm = np.zeros(proj_number)
    sigma_v = np.ones(proj_number)
    angle = np.zeros(proj_number)
    for i in np.arange(proj_number): # loop over all angles in 0.1° steps 
        # projection of data and normation
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = start_angle + 90. # projetion angle
        P.is_rescale = True #for subdivision      
        P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2) # project ones as normation
        proj_vert = P.project(img_vert)[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm_v = np.nan_to_num(proj_vert/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data 
        max_x_val = max(x_val)   # normation to the max value makes it easier converging the fit all the time 
        x_val = x_val/max_x_val;
        
        cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision) # cause of projecting a square data set 
        start_params = [a,b,c,d] # Startparameter for errorfit                  # left and right zeros are added get rid of them    
        # fitting datasets 
        best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff:-cutoff], proj_norm_v[cutoff:-cutoff], p0=start_params)
        fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
        sigma_v[i] = best_fitparam[2]
        start_angle += .1
        angle[i] = start_angle
        print fwhm[i]
    plt.figure('fwhm_min_vert')
    plt.plot(angle,fwhm)
    plt.figure('sigma_min_vert')
    plt.plot(angle,sigma_v)
    best_angle = angle[np.argmin(fwhm)]
else:
    best_angle_v = angle_tilt_vert
#============================================================================================================
#                   find best projection angle horizontal edge 
#============================================================================================================
if proj_angle_h:
    start_angle = angle_tilt_hor - .5
    end_angle = angle_tilt_hor + .5
    proj_number = 1+(end_angle -start_angle)/.1

 ### take pictures with higher power -> bigger spot -> easier to find the right angle. Choose 80  watts
    flatfield = e17.io.h5read(flatfields[15])["raw_data"]
    img_hor = e17.io.h5read(data_hor[6])["raw_data"]
    img_hor = (img_hor/(flatfield*1.))[335:525,210:400]
    pic_shape = img_hor.shape[0]
    fwhm = np.zeros(proj_number)
    sigma_h = np.ones(proj_number)
    angle = np.zeros(proj_number)
    for i in np.arange(proj_number): # loop over all angles in 0.1° steps 
        # projection of data and normation
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = start_angle  # projetion angle
        P.is_rescale = True #for subdivision      
        P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2) # project ones as normation
        proj_hor = P.project(img_hor)[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm_h = np.nan_to_num(proj_hor/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data 
        max_x_val = max(x_val)   # normation to the max value makes it easier converging the fit all the time 
        x_val = x_val/max_x_val;
        
        cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision) # cause of projecting a square data set 
        start_params = [a,b,c,d] # Startparameter for errorfit                  # left and right zeros are added get rid of them    
        # fitting datasets 
        best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff:-cutoff], proj_norm_h[cutoff:-cutoff], p0=start_params)
        fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
        sigma_h[i]= best_fitparam[2]
        start_angle += .1
        angle[i] = start_angle
        print fwhm[i]
    plt.figure('fwhm_min_hor')
    plt.plot(angle,fwhm)
    plt.figure('sigma_min_hor')
    plt.plot(angle,sigma_h)
    best_angle = angle[np.argmin(fwhm)]
else:
    best_angle_h = angle_tilt_hor
#============================================================================================================
#                    Analysis for vertical edge
#============================================================================================================
if fit_psf_v:
    vertical_angle = angle_tilt_vert - 90.
    P = ct.P.Parallel_AWP()
    P.angles = vertical_angle
    P.is_rescale = True #for subdivision
    
    mean_psf_hor = np.ones((powersteps))
    psf_hor = np.ones((powersteps,pic_number))
    fit_params_v = np.ones((powersteps,3,4))    
    for j in range(powersteps):
        #flatfield = e17.io.h5read(data+sam_v_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_v +j*12))["raw_data"]           
        flatfield = e17.io.h5read(flatfields[j])["raw_data"]
        fwhm_vert_pix = np.ones((pic_number))
        for i in np.arange(pic_number):                
            #img_vert = e17.io.h5read(data+sam_v_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_v +1+j*12+i))["raw_data"]
            img_vert = e17.io.h5read(data_vert[j])["raw_data"]
            #1/0
            img_vert = (img_vert/(flatfield*1.))[370:560,135:325]### for 60kvp control
            #img_vert = (img_vert/(flatfield*1.))[470:640,125:295]
            pic_shape = img_vert.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_vert = P.project(img_vert)[0]/np.sqrt(2)
            proj_norm_vert = np.nan_to_num(proj_vert/norm)
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = detector_pixel_size  * pixel_subdivision* np.sqrt(2) # ????????????????????????????
            x_val = np.arange(pic_shape*int(1/pixel_subdivision))*1.* projected_pixel_size
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;
            cutoff =  int((pic_shape- pic_shape/np.sqrt(2) )+.5)*int(1/pixel_subdivision)            
            
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff], proj_norm_vert[cutoff :-cutoff], p0=start_params)
            esf_fit_v = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and adapt fitparameter 
            
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            if edge_fit:
                plt.figure('edge + fit,Gauss vert_egde')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff],proj_norm_vert[cutoff :-cutoff])
                plt.plot(x_val[cutoff :-cutoff],esf_fit_v[cutoff :-cutoff])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2] ), best_fitparam[3]))
                
            fwhm_vert_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
            fit_params_v[j,i] = best_fitparam
        psf_hor[j,0] = fwhm_vert_pix[i]
        
    #plt.savefig(filepath_pic+'Edge_fit and Gauss vert_edge.tiff', format = 'tiff', dpi = 300)
    for x in range (powersteps):
        #mean_psf_hor[x] = (psf_hor[0]+ psf_hor[1]+ psf_hor[2]+psf_hor[3]+psf_hor[4])/5.
        #mean_psf_hor[x] = (psf_hor[0]+ psf_hor[1]+ psf_hor[2]+psf_hor[3]+psf_hor[4]+psf_hor[5])/6.
       # mean_psf_hor[x] = (psf_hor[0]+ psf_hor[1]+ psf_hor[2]+psf_hor[3]+psf_hor[4]+psf_hor[5]+psf_hor[6]+psf_hor[7]+psf_hor[8])/9##for 60kvp control
        mean_psf_hor[x] = (psf_hor[0]+ psf_hor[1]+ psf_hor[2]+psf_hor[3]+psf_hor[4]+psf_hor[5]+psf_hor[6]+psf_hor[7])/8##for 60kvp control

#============================================================================================================
#                   Analysis for horizontal edge
#============================================================================================================
if fit_psf_h:
    horizontal_angle = angle_tilt_hor
    P = ct.P.Parallel_AWP()
    P.angles = horizontal_angle
    P.is_rescale = True #for subdivision    
    
    mean_psf_vert = np.ones((powersteps))
    psf_vert = np.ones((powersteps,pic_number))
    fit_params_h = np.ones((powersteps,3,4))
    for j in range(powersteps):
        #flatfield = e17.io.h5read(data+sam_h_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_h +j*12))["raw_data"]           
        #flatfield = e17.io.h5read(flatfields[9+j])["raw_data"]###for 60 kvp control
        flatfield = e17.io.h5read(flatfields[j])["raw_data"]             
        fwhm_hor_pix = np.ones((pic_number))
        for i in np.arange(pic_number):               
            #img_hor = e17.io.h5read(data+sam_h_folder+"/paxscan/ct/paximage_ct_%06d.h5"%(flatfield_h +1+j*12+i))["raw_data"]
            img_hor = e17.io.h5read(data_hor[j])["raw_data"]
            #1/0
            img_hor = (img_hor/(flatfield*1.))[335:525,210:400]# for 60kvp control
            #img_hor = (img_hor/(flatfield*1.))[440:610,240:410]
            pic_shape = img_hor.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_hor = P.project(img_hor)[0]/np.sqrt(2)
            proj_norm_hor = np.nan_to_num(proj_hor/norm)
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = detector_pixel_size  * pixel_subdivision* np.sqrt(2)#??????????????
            x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision))
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;            
            cutoff = int((pic_shape- pic_shape/np.sqrt(2) )+.5)*int(1/pixel_subdivision)

            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets 
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff ], proj_norm_hor[cutoff :-cutoff ], p0=start_params)
            esf_fit_h = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and amean_psf_vert = np.ones((powersteps))dapt fitparameter 
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            
            if edge_fit:
                plt.figure('edge + fit, Gauss hor_edge')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff ],proj_norm_hor[cutoff :-cutoff ])
                plt.plot(x_val[cutoff :-cutoff ],esf_fit_h[cutoff :-cutoff ])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2]), best_fitparam[3]))
            fwhm_hor_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])           
            fit_params_h[j,i] = best_fitparam
        psf_vert[j,0] = fwhm_hor_pix[0]
            
    #plt.savefig(filepath_pic+'Edge_fit and Gauss hor_edge.tiff'[370:560,220:420], format = 'tiff', dpi = 300)
    for x in range (powersteps):
        #mean_psf_vert[x] = (psf_vert[0]+ psf_vert[1]+ psf_vert[2]+psf_vert[3]+psf_vert[4])/5. #test
        #mean_psf_vert[x] = (psf_vert[0]+ psf_vert[1]+ psf_vert[2]+psf_vert[3]+psf_vert[4]+psf_vert[5])/6.
        #mean_psf_vert[x] = (psf_vert[0]+ psf_vert[1]+ psf_vert[2]+psf_vert[3]+psf_vert[4]+psf_vert[5]+psf_vert[6]+psf_vert[7]+psf_vert[8])/9#for 60kvp control
        mean_psf_vert[x] = (psf_vert[0]+ psf_vert[1]+ psf_vert[2]+psf_vert[3]+psf_vert[4]+psf_vert[5]+psf_vert[6]+psf_vert[7])/8.#for 60kvp control
#============================================================================================================
#                   Plots and data output
#============================================================================================================
#watts = np.asarray([5,10,20,50,100,150])
watts = np.asarray([5,10,15,20,25,30,50,100,150])## for 60 kvp control
# bring data in right shape for output
if plots:
    x_max_h = (1.1*np.max(watts))
    y_max_h = (1.1*np.max(psf_vert))
    z = [ 0., x_max_h,  290., y_max_h]
    x_max_v = (1.1*np.max(watts))
    #y_max_v = (1.05*np.max(psf_hor))
    y_max_v = (1.05*470.29635215683805)##of similar shape as fir 60kvp
    v = [ 0., x_max_v,  320., y_max_v]
    plt.figure('Results Detector PSF @ 60 kVp')
    plt.subplots_adjust(wspace = 0.3)   
    plt.subplot(121,aspect = 1)
    plt.axis(z)
    plt.plot(watts, psf_vert, lw=1.5) ###for 60kvp control 470.29635215683805
    #plt.plot(watts, psf_vert,'g', lw=1.5)
    plt.plot(watts, mean_psf_vert, lw=1, ls = ':', color = 'k')
    plt.plot(watts, mean_psf_vert, lw=0 )
    #plt.title('Vertical PSF for 60kvp',fontsize = 10)
    plt.xticks(fontsize  = 12)
    plt.yticks(fontsize  = 12)
    plt.ylabel('PSF [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend(('Detector PSF vertical','Mean Detector PSF',str(round(mean_psf_vert[0]))+'$\mu$m'),fontsize = 10,loc = 0)
    plt.subplot(122, aspect = 1)
    #plt.autoscale(False)
    plt.axis(v)
    plt.plot(watts, psf_hor, lw=1.5 ) ## for 60 kvp control
    #plt.plot(watts, psf_hor,'g', lw=1.5 )
    plt.plot(watts, mean_psf_hor, lw=1, ls = ':', color = 'k')
    plt.plot(watts, mean_psf_hor, lw=0 )
    #plt.title('Horizontal PSF 60kvp',fontsize = 10)
    plt.xticks(fontsize  = 12)
    plt.yticks(fontsize  = 12)
    plt.ylabel('PSF [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend(('Detector PSF horizontal','Mean Detector PSF', str(round(mean_psf_hor[0]))+'$\mu$m'),fontsize = 10,loc = 0)
    #plt.savefig(filepath_pic+'Results Detector PSF for 60kvp_nodrop.pdf', format = 'pdf',dpi = 300)
    plt.savefig(filepath_pic+'Results Detector PSF for 60kvp_control_nodrop.pdf', format = 'pdf',dpi = 300)
if saveon:
    output_data = np.ones((3,powersteps))
    for y in range(powersteps):
        output_data[0,y] = watts[y]
        output_data[1,y] = psf_hor[y]
        output_data[2,y] = psf_vert[y]
    output_data = np.swapaxes(output_data,0,1)
    try:
        os.makedirs(filepath_data)    
    except:
        print('Folder already exist')
    try:
        np.savetxt(filepath_data+filename, output_data,fmt = '%10.5f',delimiter = ',',header = info1)
        np.savetxt(filepath_data+filenamemean, (mean_psf_vert[0],mean_psf_hor[0]),fmt = '%10.5f',delimiter = ',',header = info2)
    except:
        print('Error, could not save the file (does the script/terminal have write access?)')

#plt.figure('Edge for 150W 60 kV')
#plt.imshow(img_vert)
#plt.savefig('/users/Baier/analysis_pictures/'+'Edge for 150W 60kV Detector', format = 'tiff')    








