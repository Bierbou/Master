# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 12:34:59 2014

"""

# space for notations and other things

###tilt in projection yields from drift in flatfield and the correction is not exact

"""
@author: Markus Baier
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
import PIL
import pyfits as fit
#import fit_ellipse


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
data_folder = "linespread_acquisition_flat_measurement3_60kvphighpower"
"""   Vertical edges """
flatfield_v = 1950527 # picturenumber of first flat of vertical direction
sam_v_folder = data+data_folder+"/vertical_edge/pilatus/" # data folder 
sub_dirs_v=os.listdir(sam_v_folder)
flats_v =list()
data_v =list() ### the real pictures
dir_num=0
for sub_dir_v in sub_dirs_v:    ##### check the subfolders and sort all pictures in it in the data folder 
    dirs=np.sort(os.listdir(sam_v_folder+ sub_dir_v)) ################################daten löschen die die s zerhaut hat also die 
    for dir_name in dirs:                                       ##################### da wo die platte voll war
        if dir_num%4!=0:         
            data_v.append(sam_v_folder+ sub_dir_v+'/'+dir_name)
        else:    
            flats_v.append(sam_v_folder+ sub_dir_v+'/'+dir_name)
        dir_num=dir_num+1
data_vcopy = list(data_v)

"""   Horizontal edges """
flatfield_h = 1950531 # picturenumber of first flat
sam_h_folder = data+data_folder+"/horizontal_edge/pilatus/" # data folder of horizontal direction
sub_dirs_h=os.listdir(sam_h_folder)
flats_h=list()
data_h=list() ### the real pictures
dir_num=0
for sub_dir_h in sub_dirs_h:   ##### check the subfolders and sort all pictures in it in the data folder  
    dirs=np.sort(os.listdir(sam_h_folder+ sub_dir_h))
    for dir_name in dirs:
        if dir_num%4!=0:         
            data_h.append(sam_h_folder+sub_dir_h+'/'+dir_name)
        else:    
            flats_h.append(sam_h_folder+sub_dir_h+'/'+dir_name)
        dir_num=dir_num+1
data_hcopy = list(data_h)
"""   Diagonal edges """
flatfield_d = 1950535# picturenumber of first flat
sam_d_folder = data+data_folder+"/diagonal_edge/pilatus/" # data folder of diagonal direction
sub_dirs_d=os.listdir(sam_d_folder)
flats_d=list()
data_d=list() ### the real pictures
dir_num=0
for sub_dir_d in sub_dirs_d:    ##### check the subfolders and sort all pictures in it in the data folder 
    dirs=np.sort(os.listdir(sam_d_folder+ sub_dir_d))
    for dir_name in dirs:
        if dir_num%4!=0:         
            data_d.append(sam_d_folder+sub_dir_d+'/'+dir_name)
        else:    
            flats_d.append(sam_d_folder+sub_dir_d+'/'+dir_name)
        dir_num=dir_num+1
data_dcopy = list(data_d)


#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = data_folder+'.txt'
# write the header for txt file
info = 'contains at first column the different Watts, then vertical, horizontal, diagonal mean values of the edge measurement'

#============================================================================================================
#                       Analysis options
#============================================================================================================
# determine the exact projection angle (set True) of you know it (set False)
proj_angle = False
#choose your edge directions
fit_esf_v = True
fit_esf_h = True
fit_esf_d = False
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
power_range1 = 20
power_range2 = 40
stepsize1 = 5
stepsize2 = 10   
pic_number = 3 # Number of pictures analyzing per power
powersteps = int((power_range1/stepsize1)+(power_range2/stepsize2)+2)
angle_tilt = 2.64999 ## be carefull projector needs reflected angle v###tilt in projection yields from drift in flatfield and the correction is not exactalues
pixel_subdivision = 0.01
setup_len = 182.145
####### Parameter for the errorfit ######
a = 1. # Amplitude                                                                              #######without max_x_val 
b = 0.5 # offset in x-direction                                                                 ####### b = 20000
c = 0.5#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge      ####### c = 337
d = 0. # offset in y-direction
## Middle of the edge depth
edge_dist_z1 = 9.7
edge_dist_z2 = edge_dist_z1+19.
edge_dist_z3 = edge_dist_z1+38.
##### MAGNIFICATIONS FOR DIFFERENT DISTANCES #####
magnification_z0 = setup_len/edge_dist_z1
magnification_z1 = setup_len/edge_dist_z2
magnification_z2 = setup_len/edge_dist_z3
#### magnification correction ####
m0 = (magnification_z0/(magnification_z0-1))**2
m1 = (magnification_z1/(magnification_z1-1))**2
m2 = (magnification_z2/(magnification_z2-1))**2
#Projector settings
detector_pixel_size = 172.

#============================================================================================================
#                   find best projection angle 
#============================================================================================================
if proj_angle:
    start_angle = angle_tilt - .5
    end_angle = angle_tilt + .5
    proj_number = 1+(end_angle -start_angle)/0.1

    #for cbf pilatus format    
    flatfield,meta =  e17.io.cbfread(data+sam_v_folder+"/pilatus/ct/specadm_1_ct_%i.cbf"%(1950551))
    flatfield[:,242] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])/3 
    flatfield[:,243] = flatfield[:,241] +2*(flatfield[:,245]-flatfield[:,241])/3
    flatfield[:,244] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])
    img_vert,meta = e17.io.cbfread(data+sam_v_folder+"/pilatus/ct/specadm_1_ct_%i.cbf"%(1950553))
    #1/0    
    img_vert = (img_vert/(flatfield*1.))[215:385,120:290]
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
#        print fwhm[i]
#    plt.figure('fwhm_min')
#    plt.plot(angle,fwhm)
    best_angle = angle[np.argmin(fwhm)]
else:
    best_angle = angle_tilt
#============================================================================================================
#                    Analysis for vertical edge
#============================================================================================================
if fit_esf_v:
    #vertical_angle = angle_tilt + 90.
    vertical_angle = best_angle - 90.
    P = ct.P.Parallel_AWP()
    P.angles = vertical_angle
    P.is_rescale = True #for subdivision
    
    spot_size_horizontal = np.ones((powersteps,pic_number))
    #fwhm_vert = np.ones((powersteps,pic_number))
    mean_spot_size_horizontal = np.ones((powersteps))
    fit_params_v = np.ones((powersteps,3,4))
                                        
    for j in range(powersteps):         
        #for cbf pilatus format    
        flatfield,meta =  e17.io.cbfread(flats_v[j])
        flatfield[:,242] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])/3 
        flatfield[:,243] = flatfield[:,241] +2*(flatfield[:,245]-flatfield[:,241])/3
        flatfield[:,244] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])      
       
        fwhm_vert_pix = np.ones((pic_number))
        for i in np.arange(pic_number):
            img_vert,meta = e17.io.cbfread(data_vcopy.pop(0))
            #1/0
            img_vert = (img_vert/(flatfield*1.))[215:385,120:290]
            #1/0
            pic_shape = img_vert.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_vert = P.project(img_vert)[0]/np.sqrt(2)
            proj_norm_vert = (proj_vert/norm)
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = (detector_pixel_size * np.sqrt(2)) * pixel_subdivision
            x_val = np.arange(pic_shape*int(1/pixel_subdivision))*1.* projected_pixel_size
            #1/0
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;
            
            cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision)
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff], proj_norm_vert[cutoff :-cutoff], p0=start_params)
            esf_fit = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and adapt fitparameter             
            x_val = x_val*max_x_val
            #1/0
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
        #1/0
#        fwhm_vert[j,0] = fwhm_vert_pix[0]
#        fwhm_vert[j,1] = fwhm_vert_pix[1]
#        fwhm_vert[j,2] = fwhm_vert_pix[2]        
        spot_size_horizontal[j,0] = fwhm_vert_pix[0]/(magnification_z0-1.)
        spot_size_horizontal[j,1] = fwhm_vert_pix[1]/(magnification_z1-1.) 
        spot_size_horizontal[j,2] = fwhm_vert_pix[2]/(magnification_z2-1.)
        


        mean_spot_size_horizontal[j] = (spot_size_horizontal[j,0] + spot_size_horizontal[j,1] + spot_size_horizontal[j,2])/3.
#        print "FF Series", j
#        print 'FWHM horizontal ==>  ', spot_width_horizontalM[j,0], '[micron]'
#        print 'FWHM horizontal ==>  ', spot_width_horizontalM[j,1], '[micron]'
#        print 'FWHM horizontal ==>  ', spot_width_horizontalM[j,2], '[micron]'
#============================================================================================================
#                   Analysis for horizontal edge
#============================================================================================================
if fit_esf_h:
    horizontal_angle = best_angle
    P = ct.P.Parallel_AWP()
    P.angles = horizontal_angle
    P.is_rescale = True #for subdivision    
    
    spot_size_vertical = np.ones((powersteps,pic_number))
    mean_spot_size_vertical = np.ones((powersteps))
    fit_params_h = np.ones((powersteps,3,4))
    for j in range(powersteps):
        #for cbf pilatus format    
        flatfield,meta =  e17.io.cbfread(flats_h[j])
        flatfield[:,242] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])/3 
        flatfield[:,243] = flatfield[:,241] +2*(flatfield[:,245]-flatfield[:,241])/3
        flatfield[:,244] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])      

        fwhm_hor_pix = np.ones((pic_number))
        for i in np.arange(pic_number):
            img_hor,meta = e17.io.cbfread(data_hcopy.pop(0))
            #1/0
            img_hor = (img_hor/(flatfield*1.))[10:180,60:230]
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
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets 
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff ], proj_norm_hor[cutoff :-cutoff ], p0=start_params)
            esf_fit = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and adapt fitparameter 
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            if edge_fit:
                plt.figure('edge + fit, Gauss hor_edge')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff ],proj_norm_hor[cutoff :-cutoff ])
                plt.plot(x_val[cutoff :-cutoff ],esf_fit[cutoff :-cutoff ])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2] ), best_fitparam[3]))

            fwhm_hor_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])            
            fit_params_h[j,i] = best_fitparam
        #1/0
        spot_size_vertical[j,0] = (fwhm_hor_pix[0]-detector_pixel_size)/(magnification_z0-1.)
        spot_size_vertical[j,1] = (fwhm_hor_pix[1]-detector_pixel_size)/(magnification_z1-1.)
        spot_size_vertical[j,2] = (fwhm_hor_pix[2]-detector_pixel_size)/(magnification_z2-1.)
                
        mean_spot_size_vertical[j] = (spot_size_vertical[j,0] + spot_size_vertical[j,1] + spot_size_vertical[j,2])/3.
#============================================================================================================
#                   Analysis for diagonal edge
#============================================================================================================
if fit_esf_d:
    diagonal_angle = best_angle + 45.
    P = ct.P.Parallel_AWP()
    P.angles = diagonal_angle
    P.is_rescale = True #for subdivision
    
    spot_size_diagonal = np.ones((powersteps,pic_number))
    mean_spot_size_diagonal = np.ones((powersteps))
    fit_params_d = np.ones((powersteps,3,4))
    for j in range(powersteps):
        #for cbf pilatus format    
        flatfield,meta =  e17.io.cbfread(flats_d[j])
        flatfield[:,242] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])/3 
        flatfield[:,243] = flatfield[:,241] +2*(flatfield[:,245]-flatfield[:,241])/3
        flatfield[:,244] = flatfield[:,241] +(flatfield[:,245]-flatfield[:,241])                
            
        fwhm_diag_pix = np.ones((pic_number))
        for i in np.arange(pic_number):
            img_diag,meta = e17.io.cbfread(data_dcopy.pop(0))
            #1/0            
            img_diag = (img_diag/(flatfield*1.))[40:190,80:230]
            #1/0
            pic_shape = img_diag.shape[0]
            P.proj_area_size = (pic_shape*int(1/pixel_subdivision), 1)
            
            # projection of data and normation
            norm = P.project(np.ones((pic_shape,pic_shape)))[0]/np.sqrt(2)
            proj_diag = P.project(img_diag)[0]/np.sqrt(2)
            proj_norm_diag = np.nan_to_num(proj_diag/norm)
            #1/0
            # generation of x-values for plot and normation for fitting
            projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision
            x_val = projected_pixel_size * np.arange(pic_shape*int(1/pixel_subdivision))
            max_x_val = max(x_val)
            x_val = x_val/max_x_val;
            
            cutoff = int((pic_shape- pic_shape/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision)
            start_params = [a,b,c,d] # Startparameter for errorfit
            # fitting datasets 
            best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_val[cutoff :-cutoff ], proj_norm_diag[cutoff :-cutoff ], p0=start_params)
            esf_fit = error_fit_func(x_val, *best_fitparam)
            # convert x-values back to real and adapt fitparameter 
            x_val = x_val*max_x_val
            best_fitparam *= max_x_val
            if edge_fit:
                plt.figure('edge + fit, Gauss diag_edge')
                plt.subplot(211)
                plt.plot(x_val[cutoff :-cutoff ],proj_norm_diag[cutoff :-cutoff ])
                plt.plot(x_val[cutoff :-cutoff ],esf_fit[cutoff :-cutoff ])
                plt.subplot(212)
                plt.plot(x_val[cutoff :-cutoff ], gauss(x_val[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2]), best_fitparam[3]))

            fwhm_diag_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
            fit_params_d[j,i] = best_fitparam
        #1/0
        spot_size_diagonal[j,0] = (fwhm_diag_pix[0]-detector_pixel_size)/(magnification_z0-1.)# - (fwhm_det_psf_hor/magnification_z0)
        spot_size_diagonal[j,1] = (fwhm_diag_pix[1]-detector_pixel_size)/(magnification_z1-1.)# - (fwhm_det_psf_hor/magnification_z1)
        spot_size_diagonal[j,2] = (fwhm_diag_pix[2]-detector_pixel_size)/(magnification_z2-1.)# - (fwhm_det_psf_hor/magnification_z2)
        
        
        mean_spot_size_diagonal[j] = (spot_size_diagonal[j,0] + spot_size_diagonal[j,1] + spot_size_diagonal[j,2])/3.      
#============================================================================================================
#                   Plots and data output
#============================================================================================================
range1 = np.array(np.arange(80,(5*power_range1+  stepsize1),stepsize1))
range2 = np.array(np.arange(2*power_range2+ 3* stepsize2,150+ stepsize2,stepsize2))
watts = np.hstack((range1,range2))
# bring data in right shape for output
if plots:
    if fit_esf_v:
        x_max_v = (1.1*np.max(watts))
        y_max_v = (1.4*np.max(spot_size_horizontal))
        v = [ 75., x_max_v,  75., y_max_v]
        plt.figure('Results HORIZONTAL spot-size 60kvp')
        plt.axis(v)
        plt.plot(watts, spot_size_horizontal, lw=2 )
        plt.plot(watts, mean_spot_size_horizontal, lw=2, ls = '--',color = 'k')
        plt.title('Horizontal spot size 60kvp')
        plt.ylabel('Spot size [micron]')
        plt.xlabel('Power [W]')
        plt.legend((( magnification_z0, magnification_z1, magnification_z2, 'Mean spot size')) ,loc = 0,title = 'MAGNIFICATION')
    if fit_esf_h:
        x_max_h = (1.1*np.max(watts))
        y_max_h = (1.4*np.max(spot_size_vertical))
        h = [ 75., x_max_h,  75., y_max_h]
        plt.figure('Results VERTICAL spot-size 60kvp')
        plt.axis(h)
        plt.plot(watts, spot_size_vertical, lw=2)
        plt.plot(watts, mean_spot_size_vertical, lw=2, ls = '--', color = 'k')
        plt.title('Vertical spot size 60kvp')
        plt.ylabel('Spot size [micron]')
        plt.xlabel('Power [W]')
        plt.legend((( magnification_z0, magnification_z1, magnification_z2, 'Mean spot size')) ,loc = 0,title = 'MAGNIFICATION')
    if fit_esf_d:
        plt.figure('Results DIAGONAL spot-size 60kvp')
        plt.plot(watts, spot_size_diagonal, lw=2)
        plt.plot(watts, mean_spot_size_diagonal, lw=2, ls = '--', color = 'k')
        plt.title('Diagonal spot size 60kvp')
        plt.ylabel('Spot size [micron]')
        plt.xlabel('Power [W]')
        plt.legend((( magnification_z0, magnification_z1, magnification_z2, 'Mean spot size')) ,loc = 0,title = 'MAGNIFICATION')
        
#plt.figure('Edge for 80W 60 kV M=18.778')
#plt.imshow(img_vert)
#plt.savefig('/users/Baier/analysis_pictures/'+ 'Pilatus Edge for 80W 60kV_M=18.778', format = 'tiff')        
if saveon:
    output_data = np.ones((4,powersteps))
    for y in range(powersteps):
        output_data[0,y] = watts[y]
        output_data[1,y] = mean_spot_size_vertical[y]
        output_data[2,y] = mean_spot_size_horizontal[y]
        output_data[3,y] = mean_spot_size_diagonal[y]
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
        plt.imshow(img_vert)
        plt.savefig('/users/Baier/analysis_pictures/'+'Edge for 150W 60kV', format = 'tiff',dpi = 300)
