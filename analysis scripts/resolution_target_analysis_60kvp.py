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
import pyCT as ct
import os
sys.path.append("/vault/iSAXS/cSAXS_e14841_tomoSAXS_2014_05/processing_tooth_florian/")
import SAXS_Tools as ST
#import fit_ellipse
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
data = True
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
filepath = '/users/Baier/analysis_files/'
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
angle_tilt = -0.01
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
        #####stack all pics togehter
  ### you have to stack the piccs together but first skip 27 of the right side of the left pic and 27 from the left side of the right
            #### pic
            ### merge them together and take the mean
            # after that stack the cutted pictures with inbetzween the te mean
        #1/0   
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
        pic_change_h[j] = np.vstack((horizontalpic_fill,(img_h_comp[j][290:511,:])))         

    
     

#============================================================================================================
#                  Data projection  
#============================================================================================================
if proj_data:
   
    pic_shape_v_change = pic_change_v[0].shape[0]
    pic_shape_h_change = pic_change_h[0].shape[0]
    proj_vert = np.ones((powersteps,pic_shape_v_change*int(1/pixel_subdivision)))
    proj_norm_v = np.ones((powersteps,pic_shape_v_change*int(1/pixel_subdivision)))
    
    proj_hor = np.ones((powersteps,pic_shape_h_change*int(1/pixel_subdivision)))
    proj_norm_h = np.ones((powersteps,pic_shape_h_change*int(1/pixel_subdivision)))
    for k in np.arange(powersteps): # loop over all angles in 0.1° steps 
        # projection of data and normation for the vertical target
        P = ct.P.Parallel_AWP()  # projector 
        P.angles = angle_tilt# + 90. # projetion angle
        P.is_rescale = False#for subdivision      
        P.proj_area_size = (pic_shape_v_change*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = P.project(np.ones((pic_shape_v_change,pic_shape_v_change)))[0]/np.sqrt(2) # project ones as normation
        proj_vert[k] = P.project(pic_change_v[k])[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm_v[k] = (proj_vert[k]/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size*np.arange(pic_shape_v_change*int(1/pixel_subdivision)) # projected_pixel_size *  x-values for fitting and plotting onto the data
        
#        plt.figure('angle -0.01 subdiv 0.1')
#        plt.plot(x_val,proj_norm)

        # projection of data and normation for the horizontal target
        Q = ct.P.Parallel_AWP()  # projector 
        Q.angles = angle_tilt + 90. # projetion angle
        Q.is_rescale = True #for subdivision      
        Q.proj_area_size = (pic_shape_h_change*int(1/pixel_subdivision), 1) # size after the projected data are streched to 
        norm = Q.project(np.ones((pic_shape_h_change,pic_shape_h_change)))[0]/np.sqrt(2) # project ones as normation
        proj_hor[k] = Q.project(pic_change_h[k])[0]/np.sqrt(2) # project the data set under the angle the edge was tilted 
        proj_norm_h[k] = (proj_hor[k]/norm) # normation of the data 
        # generation of x-values for plot and normation for fitting
        projected_pixel_size = detector_pixel_size * np.sqrt(2) * pixel_subdivision # pixelsize after the oversampling 
        x_val = projected_pixel_size * np.arange(pic_shape_h_change*int(1/pixel_subdivision)) # x-values for fitting and plotting onto the data
        
#        plt.figure('angle -0.01 subdiv 0.1')
#        plt.plot(x_val,proj_norm)

#        
#        cutoff = int((pic_shape_h_change- pic_shape_h_change/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision) # cause of projecting a square data set 
                                                        # left and right zeros are added get rid of them      
#transform_v = fft(proj_vert[0])
#ST.savitzky_golay(zeroabs,55,3)     



                                              
#### Plots and other output stuff                                                   
#plt.close('all')
#########pics for vertical res target for testing#############################
#plt.figure('test_v')
#plt.imshow(img_v_comp[0])  
#plt.figure('change_v')
#plt.imshow(pic_change_v[0],vmin =0)
#########pics for horizontal res target for testing################################
#plt.figure('test_h')
#plt.imshow(img_h_comp[0],vmin = 0)
#plt.figure('change_h')
#plt.imshow(pic_change_h[0],vmin =0)


#project them like all the time before
#do some fourier trafo





