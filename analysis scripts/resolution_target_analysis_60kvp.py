# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 13:47:58 2015

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
sys.path.append("/vault/iSAXS/cSAXS_e14841_tomoSAXS_2014_05/processing_tooth_florian/")
import SAXS_Tools as ST
sys.path.append('/users/schaff/Python Scripts/')
"""
Script determing the witdth of focal size and shape
with line-spread-function analysis
"""
plt.close("all")
#============================================================================================================
#                       Analysis options
#============================================================================================================
data_input = True
proj_data = True
min_max = True
data = True
saveon = False
plots = True
#============================================================================================================
#                        Data input
#============================================================================================================
if data:
    """ Standard data folder"""
    data = "/data/DPC/local_setups/microfocus/samples/"
    """ Folder containing the data"""
    data_folder = "resolution_target_60kvp"
    """   Vertical resw target """
    flatfield_v = 195091 # picturenumber of first flat of vertical direction
    sam_v_folder = data+data_folder+"/vertical_target/paxscan/" # data folder 
    sub_dirs_v=os.listdir(sam_v_folder)
    flats_v =list()
    data_v =list() ### the real pictures
    dir_num=0
    for sub_dir_v in sub_dirs_v:    ##### check the subfolders and sort all pictures in it in the data folder 
        dirs=np.sort(os.listdir(sam_v_folder+ sub_dir_v)) 
        for dir_name in dirs:                                       
            if dir_num%9!=0:         
                data_v.append(sam_v_folder+ sub_dir_v+'/'+dir_name)
            else:    
                flats_v.append(sam_v_folder+ sub_dir_v+'/'+dir_name)
            dir_num=dir_num+1
    data_vcopy = list(data_v)
    
    """   Horizontal res target """
    flatfield_h = 195100 # picturenumber of first flat
    sam_h_folder = data+data_folder+"/horizontal_target/paxscan/" # data folder of horizontal direction
    sub_dirs_h=os.listdir(sam_h_folder)
    flats_h=list()
    data_h=list() ### the real pictures
    dir_num=0
    for sub_dir_h in sub_dirs_h:   ##### check the subfolders and sort all pictures in it in the data folder  
        dirs=np.sort(os.listdir(sam_h_folder+ sub_dir_h))
        for dir_name in dirs:
            if dir_num%9!=0:         
                data_h.append(sam_h_folder+sub_dir_h+'/'+dir_name)
            else:    
                flats_h.append(sam_h_folder+sub_dir_h+'/'+dir_name)
            dir_num=dir_num+1
    data_hcopy = list(data_h)

#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_pictures/restarget/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = data_folder+'.txt'
# write the header for txt file
info = 'contains at first column the different Watts, then vertical, horizontal, diagonal mean values of the edge measurement'

#============================================================================================================
#                       Function definitons
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
flatfield = e17.io.h5read(flats_v[0])["raw_data"]
pic_shape = flatfield.shape[0]
pixel_subdivision = 1.
power_range = 35
stepsize = 5
powersteps = int(power_range/stepsize + 2)
start_pictures = 8
usefullpics_v = 7
usefullpics_h = 6
### startparams
a = 1. # Amplitude                                                                              #######without max_x_val 
b = 0. # offset in x-direction                                                                 ####### b = 20000
c = 0.9#(detector_pixel_size)/((2*(2*np.log(2))**.5)*magnification[i]) # width of the edge      ####### c = 337
d = 0. # offset in y-direction

angle_tilt = 0.01
overlap_h = 27
cuttoff_h = 529
overlap_v = 23
cuttoff_v = 600 
detector_pixel_size = 127.
#============================================================================================================
#                       stiching different target pictures to one big pic
#============================================================================================================
if data_input:
    img_v_comp = np.ones((powersteps,pic_shape*8-(cuttoff_v+(overlap_v*7)),pic_shape))
    img_h_comp = np.ones((powersteps,pic_shape,pic_shape*8-(cuttoff_h+(overlap_h*7))))
    pic_change_v = np.ones((powersteps,5639,5639))
    pic_change_h = np.ones((powersteps,5682,5682))
    verticalpic_fill = np.zeros((5639,5418))
    horizontalpic_fill = np.zeros((5461,5682))
    for j in range(powersteps):
        img_v = np.ones((start_pictures,pic_shape,pic_shape))
        img_h = np.ones((start_pictures,pic_shape,pic_shape))
        merge_v = np.ones((start_pictures-1,overlap_v,pic_shape))
        merge_h = np.ones((start_pictures-1,pic_shape,overlap_h))
        flatfield_v = e17.io.h5read(flats_v[j])["raw_data"]
        flatfield_h = e17.io.h5read(flats_h[j])["raw_data"]
        for i in range(start_pictures):
            ###read in the vertical targets
            img_raw_v = e17.io.h5read(data_vcopy.pop(0))["raw_data"]
            img_v[i] = (img_raw_v/(flatfield_v*1.))
            ###read in the horizontal targets
            img_raw_h = e17.io.h5read(data_hcopy.pop(0))["raw_data"]
            img_h[i] = (img_raw_h/(flatfield_h*1.))
        
        #1/0
        #####stack all pics togehter  
        for l in range(start_pictures-1):
            merge_v[l] = (img_v[l][pic_shape-overlap_v:,:]+img_v[l+1][:overlap_v,:])/2.
            merge_h[l] = (img_h[l][:,:overlap_h]+img_h[l+1][:,pic_shape-overlap_h:])/2.
        img_v_comp[j] = np.vstack((img_v[0][:-overlap_v,:],merge_v[0],img_v[1][overlap_v:-overlap_v,:],merge_v[1],
                                   img_v[2][overlap_v:-overlap_v,:],merge_v[2],img_v[3][overlap_v:-overlap_v,:],merge_v[3],
                                   img_v[4][overlap_v:-overlap_v,:],merge_v[4],img_v[5][overlap_v:-overlap_v,:],merge_v[5],
                                   img_v[6][overlap_v:-overlap_v,:],merge_v[6],img_v[7][overlap_v:-cuttoff_v,:]))
        pic_change_v[j] = np.hstack((verticalpic_fill,-np.log(img_v_comp[j][:,250:471])))
        
        
        img_h_comp[j] = np.hstack((img_h[7][:,cuttoff_h:-overlap_h],merge_h[6],img_h[6][:,overlap_h:-overlap_h],merge_h[5],
                               img_h[5][:,overlap_h:-overlap_h],merge_h[4],img_h[4][:,overlap_h:-overlap_h],merge_h[3],
                               img_h[3][:,overlap_h:-overlap_h],merge_h[2],img_h[2][:,overlap_h:-overlap_h],merge_h[1],
                               img_h[1][:,overlap_h:-overlap_h],merge_h[0],img_h[0][:,overlap_h:]))
        pic_change_h[j] = np.vstack((horizontalpic_fill,-np.log(img_h_comp[j][290:511,:])))
#============================================================================================================
#                  Data projection  
#============================================================================================================
if proj_data:
   
    pic_shape_v_change = pic_change_v[0].shape[0]
    pic_shape_h_change = pic_change_h[0].shape[0]
    proj_vert = np.ones((powersteps,pic_shape_v_change*int(1/pixel_subdivision)))
    proj_hor = np.ones((powersteps,pic_shape_h_change*int(1/pixel_subdivision)))
    for k in np.arange(powersteps): # loop over all angles in 0.1° steps 
        # projection of data and normation for the vertical target
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = angle_tilt -180.# + 90. # projetion angle
        P.is_rescale = False#for subdivision      
        P.proj_area_size = (pic_shape_v_change*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        proj_vert[k] = P.project(pic_change_v[k])[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_vert_cut = proj_vert[:,25:5605]
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        #x_val_v = np.arange(pic_shape_v_change*int(1/pixel_subdivision)) # projected_pixel_size *  x-values for fitting and plotting onto the data
        #x_val_v_cut = x_val_v[25:5605]
        
        # projection of data and normation for the horizontal target
        Q = ct.P.Parallel_AWP()  # projector 
        Q.angles = angle_tilt + 90. # projetion angle
        Q.is_rescale = False #for subdivision      
        Q.proj_area_size = (pic_shape_h_change*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        proj_hor[k] = Q.project(pic_change_h[k])[0]/np.sqrt(2) # project the data set under the angle the edge was tilted
        proj_hor_cut = proj_hor[:,70:5620]#[:,70:5200]
        # generation of x-values for plot and normation for fitting
        #x_val_h = np.arange(pic_shape_h_change*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data
        #x_val_h_cut = x_val_h[70:5620]#[70:5200]
#============================================================================================================
#                  Min and Max determination for vertical and horizontal res target
#============================================================================================================
if min_max:
###### vertical res target                         
    maxlist_v_comp =[]
    end_list_max_v_comp =[]
    maxima_v_comp = []
    minima_v_comp = []
    for j in range(usefullpics_v):  ## code for the vertical res targets 
        max_val_v = []
        for i in range(len(proj_vert_cut[j])-1):
            if proj_vert_cut[j][i]> 46. and np.abs(proj_vert_cut[j][i] -proj_vert_cut[j][i+1]) <7. :
                max_val_v.append(proj_vert_cut[j][i])
            else:
                max_val_v.append(0)
        maxlist_v_comp.append(max_val_v)
        end_max_v = [0]
        for i in range(len(max_val_v)-1):
            if max_val_v[i+1]== 0 and max_val_v[i] != 0:
                end_max_v.append(i)
        end_list_max_v_comp.append(end_max_v)
        
        maxima_v = []
        minima_v = []
        for i in range(0,len(end_list_max_v_comp[j])-1,): #### find min and maxima  of the original curve vert
            arr_v = proj_vert_cut[j][end_list_max_v_comp[j][i]:end_list_max_v_comp[j][i+1]]
            maxima_v.append(np.max(arr_v))
            minima_v.append(np.min(arr_v))            
        maxima_v_comp.append(maxima_v)
        minima_v_comp.append(minima_v)   
    
    contrast_v_comp = []
    norm_contrast_v = []
    for j in range(usefullpics_v):
        contrast_v = []
        vnorm = []
        for i in range(len(maxima_v_comp[j])):
            contrast_v.append((maxima_v_comp[j][i]- minima_v_comp[j][i])/(maxima_v_comp[j][i]+ minima_v_comp[j][i]))
        while len(contrast_v) <= 149 :
            contrast_v.append (0.)
        contrast_v_comp.append(contrast_v)
        for i in range(len(contrast_v_comp[j])):
            vnorm.append((contrast_v_comp[j][i]/max(contrast_v_comp[j])))
        norm_contrast_v.append(vnorm)
###### horizontal res target 
    maxlist_h_comp =[]
    end_list_max_h_comp =[]
    maxima_h_comp = []
    minima_h_comp = []
    for j in range(usefullpics_v):  ## code for the horizontal res targets 
        max_val_h = []
        for i in range(len(proj_hor_cut[j])-1):
            if proj_hor_cut[j][i]>= 44.25 and np.abs(proj_hor_cut[j][i] -proj_hor_cut[j][i+1]) <7. :
                max_val_h.append(proj_hor_cut[j][i])
            else:
                max_val_h.append(0)
        maxlist_h_comp.append(max_val_h)
        end_max_h = [0]
        for i in range(len(max_val_h)-1):
            if max_val_h[i+1]== 0 and max_val_h[i] != 0:
                end_max_h.append(i)
        end_list_max_h_comp.append(end_max_h)
        
        maxima_h = []
        minima_h = []
        for i in range(0,len(end_list_max_h_comp[j])-1,): #### find min and maxima  of the original curve vert
            arr_h = proj_hor_cut[j][end_list_max_h_comp[j][i]:end_list_max_h_comp[j][i+1]]
            maxima_h.append(np.max(arr_h))
            minima_h.append(np.min(arr_h))            
        maxima_h_comp.append(maxima_h)
        minima_h_comp.append(minima_h)
    
    contrast_h_comp = []
    norm_contrast_h = []
    for j in range(usefullpics_v):
        contrast_h = []
        hnorm = []
        for i in range(len(maxima_h_comp[j])):
            contrast_h.append((maxima_h_comp[j][i]- minima_h_comp[j][i])/(maxima_h_comp[j][i]+ minima_h_comp[j][i]))
        while len(contrast_h) <= 149 :
            contrast_h.append (0.)
        contrast_h_comp.append(contrast_h)
        for i in range(len(contrast_h_comp[j])):
            hnorm.append((contrast_h_comp[j][i]/max(contrast_h_comp[j])))
        norm_contrast_h.append(hnorm)
#============================================================================================================
#                 Fwhm determination start_params = [a,b,c,d] # Startparameter for errorfit
#============================================================================================================
fwhm_v = np.ones((usefullpics_v))
fwhm_h = np.ones((usefullpics_v))
fitparam_v = np.ones((usefullpics_v,4))
fitparam_h = np.ones((usefullpics_v,4))
up_bound = 1000/4.
low_bound = 1000/32.
x_lefthalf = np.linspace(-32,-4,num = 75)
x_righthalf = np.linspace(4,32,num = 75)
x_val = np.hstack((x_lefthalf,x_righthalf))
line_pairs_1 = np.linspace(low_bound , up_bound, num = 150)
line_pairs_norm = line_pairs_1/max(line_pairs_1)
for i in range (usefullpics_v):
    start_params = [a,b,c,d] # Startparameter for errorfit
    trafo_v = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_v[i])))
    trafo_norm_v = trafo_v/max(trafo_v)
    trafo_h = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_h[i])))
    trafo_norm_h = trafo_h/max(trafo_h)
    #1/0
    # fitting datasets
#    x_val_v = np.arange(trafo_norm_v.shape[0])*1.
#    max_x_val_v = max(x_val_v)
#    x_val_norm_v = x_val_v/max_x_val_v
    x_val_norm = x_val/max(x_val) 
    #1/0
    best_fitparam_v, cov_esf = opti.curve_fit(gauss, x_val_norm, 
                                            trafo_norm_v, p0=start_params)
    #best_fitparam_v, cov_esf = opti.curve_fit(gauss, line_pairs_norm[:-3], 
     #                                       norm_contrast_v[0], p0=start_params)
    fitparam_v[i] = best_fitparam_v
    best_fitparam_v*=max(x_val)
    gauss_fit_v = gauss(x_val,*best_fitparam_v)
    plt.figure(1)
    plt.plot(x_val,gauss_fit_v)
    fwhm_v[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam_v[2])
    start_params = [a,b,c,d] # Startparameter for errorfit
    # fitting datasets 
#    x_val_h = np.arange(trafo_norm_h.shape[0])*1.
#    max_x_val_h = max(x_val_h)
#    x_val_norm_h = x_val_h/max_x_val_h
    
    best_fitparam_h, cov_esf = opti.curve_fit(gauss, x_val_norm, 
                                            trafo_norm_h, p0=start_params)
    fitparam_h[i] = best_fitparam_h
    best_fitparam_h*=max(x_val)
    gauss_fit_h = gauss(x_val,*best_fitparam_h)
    plt.figure(2)
    plt.plot(x_val,gauss_fit_h)
    fwhm_h[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam_h[2])
    #1/0
smooth_v = np.ones((usefullpics_v,(len(norm_contrast_v[0]))))
smooth_h = np.ones((usefullpics_v,(len(norm_contrast_h[1]))))
for i in range(usefullpics_v):
    smooth_v[i] = ST.savitzky_golay(norm_contrast_v[i],85,7)
    smooth_h[i] = ST.savitzky_golay(norm_contrast_h[i],65,7)
#    plt.figure('smoothv')
#    plt.plot(line_pairs_1, smooth_v[i])
#    plt.figure('smoothh')
#    plt.plot(line_pairs_1, smooth_h[i])
### check if smoothing is necessary for vertical 
#a = np.fft.fftshift(np.absolute(np.fft.fft(smooth_v[6])))
#b = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_v[6])))
#plt.plot(a),plt.plot(b)
#### check if smoothing is necessary for horizontal      
#c = np.fft.fftshift(np.absolute(np.fft.fft(smooth_h[6])))
#d = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_h[6])))
#plt.plot(c),plt.plot(d)
#### check if necessary avoiding the zeros      
#e = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_v[0])))
#f = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_v[0])))
#plt.plot(e),plt.plot(f)
#g = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_h[0])))
#h = np.fft.fftshift(np.absolute(np.fft.fft(norm_contrast_h[0])))
#plt.plot(g),plt.plot(h)
#
#
#plt.figure()
#for i in range(usefullpics_v):
#    plt.plot(smooth_h[i])
#    plt.plot(norm_contrast_h[i])
#    plt.axhline(y = .10, ls = '--', color = 'k')
#    

#============================================================================================================
#                 plots 
#============================================================================================================
up_bound = 1000/4.
low_bound = 1000/32.
line_pairs_1 = np.linspace(low_bound , up_bound, num = 150)
if plots:
    plt.figure('norm contrast vertical resolution target')
    plt.subplot(111)
    plt.plot(line_pairs_1, norm_contrast_v[0],'b+',line_pairs_1, smooth_v[0],
             line_pairs_1, norm_contrast_v[1],'g+',line_pairs_1, smooth_v[1],
             line_pairs_1, norm_contrast_v[2],'r+',line_pairs_1, smooth_v[2],
             line_pairs_1, norm_contrast_v[3],'c+',line_pairs_1, smooth_v[3],
             line_pairs_1, norm_contrast_v[4],'m+',line_pairs_1, smooth_v[4],
             line_pairs_1, norm_contrast_v[5],'y+',line_pairs_1, smooth_v[5],
             line_pairs_1, norm_contrast_v[6],'k+',line_pairs_1, smooth_v[6],lw = 2 )         
    plt.axis([low_bound,up_bound,0,1.10])
    plt.xticks( np.arange((low_bound-0.2), (up_bound+2),20), fontsize = 12)
    plt.yticks(np.arange(0,np.max(norm_contrast_v)*1.1,.1),fontsize = 12)
    #plt.title('MTFs of the vertical direction @ 60 kVp variing power ')
    plt.axhline(y = .1, ls = '--', color = 'k')
    plt.ylabel('MTF [%]')
    plt.xlabel('[lp/mm]')
    plt.legend(('3 ','' ,'5','','10','','15','','20','','25','','30','','10 % MTF'),loc = 0, title = 'Power [W]',fontsize = 9)
    #plt.savefig(filepath+'mtf vertical target 60kv.pdf', format = 'pdf', dpi = 300)
#    plt.subplot(122)
#    plt.plot(line_pairs_1, smooth_v[0],
#             line_pairs_1, smooth_v[1],
#             line_pairs_1, smooth_v[2],
#             line_pairs_1, smooth_v[3],
#             line_pairs_1, smooth_v[4],
#             line_pairs_1, smooth_v[5],
#             line_pairs_1, smooth_v[6])           
#    plt.axis([low_bound,up_bound,0,1.10])
#    plt.xticks( np.arange((low_bound-0.2), (up_bound+2),20), fontsize = 8)
#    plt.yticks(np.arange(0,np.max(smooth_v)*1.1,.1),fontsize = 8)
#    plt.title('Smoothed MTFs of the vertical direction @ 60 kVp variing power ')
#    plt.axhline(y = .1, ls = '--', color = 'k')
#    plt.ylabel('MTF [%]')
#    plt.xlabel('[lp/mm]')
#    plt.legend(('3 ','5','10','15','20','25','30','10 % MTF'),loc = 0, title = 'Power [W]')
    
    plt.figure('norm contrast horizontal resolution target')        
    plt.plot(line_pairs_1, norm_contrast_h[0],'b+',line_pairs_1, smooth_h[0],
             line_pairs_1, norm_contrast_h[1],'g+',line_pairs_1, smooth_h[1],
             line_pairs_1, norm_contrast_h[2],'r+',line_pairs_1, smooth_h[2],
             line_pairs_1, norm_contrast_h[3],'c+',line_pairs_1, smooth_h[3],
             line_pairs_1, norm_contrast_h[4],'m+',line_pairs_1, smooth_h[4],
             line_pairs_1, norm_contrast_h[5],'y+',line_pairs_1, smooth_h[5],
             line_pairs_1, norm_contrast_h[6],'k+',line_pairs_1, smooth_h[6],lw = 2)
    plt.axis([low_bound,up_bound,0,1.10])
    plt.xticks( np.arange((low_bound-0.2), (up_bound+2),20) )
    plt.yticks(np.arange(0,np.max(norm_contrast_h)*1.1,.1))
    #plt.title('MTFs of the horizontal direction @ 60 kVp variing power')
    plt.ylabel('MTF  [%]')
    plt.xlabel('[lp/mm]')
    plt.axhline(y = .1, ls = '--', color = 'k')
    plt.legend(('3 ','','5','','10','','15','','20','','25','','30','',' 10 % MTF'),loc = 0, title = 'Power [W]',fontsize = 9)
    #plt.savefig(filepath+'mtf horizontal target 60kv.pdf', format = 'pdf', dpi = 300)
##### eventually the right version for the verticatl target
#    plt.figure('vertical target for 3 first picture')
#    plt.imshow(img_v[0][:,150:750])
#    plt.xlabel('x-direction [pixel]')
#    #plt.xticks(np.arange(0,800,400))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+' first pic res target vert 3w 60kvp.pdf', format  = 'pdf', dpi = 300)
#    
#    plt.figure('horizontal target for 3 first picture')
#    plt.imshow(img_h[0][150:750,:])
#    plt.xlabel('x-direction [pixel]')
#    #plt.xticks(np.arange(0,800,400))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+' first pic res target hor 3w 60kvp.pdf', format  = 'pdf', dpi = 300) 
    
#    plt.figure('vertical target for 3 watt')
#    plt.imshow(img_v_comp[0])
#    plt.xlabel('x-direction [pixel]')
#    plt.xticks(np.arange(0,801,400))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target vert 3w 60kvp.pdf', format  = 'pdf', dpi = 300)  
#    
#    plt.figure('horizontal target for 3 watt')
#    plt.imshow(img_h_comp[0],vmin = 0,vmax = 1.5)
#    plt.xlabel('x-direction [pixel]')
#    plt.yticks(np.arange(0,800,200))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target hor 3w 60kvp.pdf', format  = 'pdf', dpi = 300)
#    
#    plt.figure('vertical target for 25 watt')
#    plt.imshow(img_v_comp[5])
#    plt.xlabel('x-direction [pixel]')
#    plt.xticks(np.arange(0,801,400))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target vert 25w 60kvp.pdf', format  = 'pdf', dpi = 300)  
#    
#    plt.figure('horizontal target for 25 watt')
#    plt.imshow(img_h_comp[5])#,vmin = 0,vmax = 1.5)
#    plt.xlabel('x-direction [pixel]')
#    plt.yticks(np.arange(0,800,200))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target hor 25w 60kvp.pdf', format  = 'pdf', dpi = 300)
#    
#    plt.figure('vertical target for 30 watt')
#    plt.imshow(img_v_comp[6])
#    plt.xlabel('x-direction [pixel]')
#    plt.xticks(np.arange(0,801,400))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target vert 30w 60kvp.pdf', format  = 'pdf', dpi = 300)  
#    
#    plt.figure('horizontal target for 30 watt')
#    plt.imshow(img_h_comp[6])#,vmin = 0,vmax = 1.5)
#    plt.xlabel('x-direction [pixel]')
#    plt.yticks(np.arange(0,800,200))
#    plt.ylabel('y-direction [pixel]')
#    plt.savefig(filepath+'whole res target hor 30w 60kvp.pdf', format  = 'pdf', dpi = 300)  
    
#    plt.figure('proj res traget vert')
#    plt.plot(proj_vert_cut[0])
#    plt.title('projection of res target vert for 3 watt')
#    plt.ylabel('Mean Intensity [arb.unit]')
#    plt.xlabel('Projected pixels')
#    plt.legend(('projected target lines',))
#    plt.savefig(filepath+'projection of res target vert for 3 watt.pdf', format = 'pdf',dpi = 300)
#    
#    plt.figure('proj res traget hor')
#    plt.plot(proj_hor_cut[0])
#    plt.title('projection of res target hor for 3 watt')
#    plt.ylabel('Mean Intensity [arb.unit]')
#    plt.xlabel('Projected pixels')
#    plt.legend(('projected target lines',))
#    plt.savefig(filepath+'projection of res target hor for 3 watt.pdf', format = 'pdf',dpi = 300)
    
    plt.figure('proj res traget vert')
    plt.plot(proj_vert_cut[6])
    plt.title('projection of res target vert for 30 watt')
    plt.ylabel('Mean Intensity [arb.unit]')
    plt.xlabel('Projected pixels')
    plt.legend(('projected target lines',))
    plt.savefig(filepath+'projection of res target vert for 30 watt.pdf', format = 'pdf',dpi = 300)
    
    plt.figure('proj res traget hor')
    plt.plot(proj_hor_cut[6])
    plt.title('projection of res target hor for 30 watt')
    plt.ylabel('Mean Intensity [arb.unit]')
    plt.xlabel('Projected pixels')
    plt.legend(('projected target lines',))
    plt.savefig(filepath+'projection of res target hor for 30 watt.pdf', format = 'pdf',dpi = 300)
#line_pairs_2 = np.linspace(low_bound , up_bound, num = 147)        
#line_pairs_3 = np.linspace(low_bound , up_bound, num = 148)                
#    plt.plot(line_pairs_2, norm_contrast_v[0],'+',line_pairs_2, norm_contrast_v[1],'+',
#             line_pairs_2, norm_contrast_v[2],'+',line_pairs_3, norm_contrast_v[3],'+',
#             line_pairs_3, norm_contrast_v[4],'+',line_pairs_1, norm_contrast_v[5],'+',
#             line_pairs_3, norm_contrast_v[6],'+')
#============================================================================================================
#                  Save pics and data 
#============================================================================================================
#if saveon:        







