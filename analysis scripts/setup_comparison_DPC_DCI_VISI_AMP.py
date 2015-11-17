# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 11:45:15 2015

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
from PIL import Image
sys.path.append('/users/schaff/Python Scripts/')



"""
Script providing the analysis of AMP,DPC,DCI and VISI signal for comparison between different setup-configurations
"""

def select_point(image,maxval,minval):
    """
    Function to click in an image and return indices
    """
    msg = "Select point with left mousebutton, confirm with left mousebutton. Delete selection with right mouse button. "
    plt.close(1)    
    plt.figure(1,figsize=(12,12))
    plt.imshow(image,cmap=cm.Greys_r,vmax = maxval,vmin = minval)
    plt.title(msg)
    plt.show()
    print msg
    p1, p2 = plt.ginput(2, timeout = 0)
    p1, p2 = (int(p1[0]+.5),int(p1[1]+.5)),(int(p2[0]+.5),int(p2[1]+.5))
    plt.close(1)
    return p1
#plt.close("all")
#============================================================================================================
#                        Analysis options
#============================================================================================================
saveon = True
plots = False
show = False
g0 = False
find_roi = False
#============================================================================================================
#                        Data input
#============================================================================================================
""" Standard data folder"""
all_data = "/data/DPC/local_setups/microfocus/samples/"
    #data_folder = "watertube_wood_z%i/" %(i+1)
""" Folder containing the data"""
#data = "/data/DPC/local_setups/microfocus/samples/Symm_G0_long/" ## contains images of the standard setup
data = "/data/DPC/local_setups/microfocus/samples/Symm_noG0_long/" ## contains images of the standard setup without G0
AMP = []
DPC = []
DCI = []
VISI = []


for i in range(6):
    if g0:
        data_folder = "watertube_wood_SiO2_05_z%i/"%(i+1) ## contains images of the standard setup
        data_files = "projection_watertube_wood_SiO2_05_z%i_info.tif"%(i+1)  ## contains images of the standard setup
    else:
        data_folder = "noG0_watertube_wood_SiO2_05_z%i/"%(i+1) ## contains images of the standard setup without G0    
        data_files = "projection_noG0_watertube_wood_SiO2_05_z%i_info.tif"%(i+1) ## contains images of the standard setup without G0
    
    data_file_AMP = "AMP_"+data_files
    data_file_DPC = "DPC_"+data_files
    data_file_DCI = "DCI_"+data_files
    data_file_VISI = "Flats_visi_"+data_files
    
    AMP.append(np.array(Image.open(data+data_folder+data_file_AMP))*1.)
    DPC.append(np.array(Image.open(data+data_folder+data_file_DPC))*1.)
    DCI.append(np.array(Image.open(data+data_folder+data_file_DCI))*1.)
    VISI.append(np.array(Image.open(data+data_folder+data_file_VISI))*1.)
    
    if show:
        plt.figure('AMP')
        plt.imshow(AMP)
        plt.figure('DPC')
        plt.imshow(DPC)
        plt.figure('DCI')
        plt.imshow(DCI)
        plt.figure('VISI')
        plt.imshow(VISI)
#1/0    
data_AMP = np.array(AMP)/1.
data_DPC = np.array(DPC)/1.
data_DCI = np.array(DCI)/1.
data_VISI = np.array(VISI)/1.
# rescale data to measured values
#(img.copy()-np.float(scale[0]))/(np.float(scale[1])-np.float(scale[0]))*(2**16-1)
original_data_DPC = (data_DPC/(2**16-1))*(2*np.pi) -np.pi
original_data_VISI= (data_VISI/(2**16-1))*(1.1) -0
original_data_DCI = (data_DCI/(2**16-1))*(1.1) -0
if find_roi:
    roi_list_DCI=np.zeros((6,4)) 
    for i in range(6):
        x_tl,y_tl=select_point(original_data_DCI[i],2,0)
        x_br,y_br=select_point(original_data_DCI[i],2,0) 
        roi_list_DCI[i,0]=x_tl
        roi_list_DCI[i,1]=y_tl
        roi_list_DCI[i,2]=x_br
        roi_list_DCI[i,3]=y_br
    roi_list_DPC=np.zeros((6,4)) 
    for i in range(6):
        x_tl,y_tl=select_point(original_data_DPC[i],np.pi,-np.pi)
        x_br,y_br=select_point(original_data_DPC[i],np.pi,-np.pi) 
        roi_list_DPC[i,0]=x_tl
        roi_list_DPC[i,1]=y_tl
        roi_list_DPC[i,2]=x_br
        roi_list_DPC[i,3]=y_br
#1/0
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/setupcomp/'
try:
    os.makedirs(filepath_pic)    
except:
    print('Folder already exist')
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = data_folder[:len(data_folder)-4]+'.txt'
# write the header for txt file
info = 'contains at first column the different relative positions between source/G0<-->G1 phase scale factor, then the phase shifts,then DCI scale factro, then the darkfield, then visibility '
"""
*(2**16-1)) divide the DPC data by this factor maybe thats it for the phase retrieval
"""

## define cut points for different sample positions
x_start  =np.array([300,300,300,300,300,300])
x_end= np.array([550,550,550,550,550,550])
y_start = np.array([300,335,360,370,380,390])
y_end = np.array([400,410,410,430,440,430]) 
cut_corrig_DPC = original_data_DPC#[:,245:550,200:560] #z1[245:550,50:560]

lineplot_pic = np.mean(original_data_DPC[5,400:500,30:220], axis = 1)
plt.figure('Phase shift line-plot')
plt.plot(lineplot_pic, label = 'averaged phase-shift')
#plt.xticks(np.arange(0,1.1,.10),fontsize = 12)
plt.yticks(np.arange(-1.,.1,.10),fontsize = 12)
#plt.xlim(0.,1.)
plt.ylim(-1,.1)
plt.xlabel('pixel vertial-direction',fontsize = 12)
plt.ylabel('Pahse-shift [pi]',fontsize = 12)
plt.legend(loc = 0,fontsize = 12 )
plt.savefig(filepath_pic+'phase_shift_example_nog0.pdf', format = 'pdf',dpi = 300 )


plt.figure('Raw phase image')
plt.imshow(data_AMP[5])
plt.xlabel('pixel horizontal-direction',fontsize = 12)
plt.ylabel('pixel vertial-direction',fontsize = 12)
plt.savefig(filepath_pic+'raw_amp_image_z5_400,500_30,220_nog0.pdf', format = 'pdf',dpi = 300 )
    #test_DPC = data_DPC[:,245:550,50:560] ### slicing for watertube_wood_info.tif
#lineplot_test = np.mean(test_DPC,axis = 1)#/(2**16-1)#,np.mean(test_DPC[34:,0:80],axis = 1),np.mean(test_DPC[34:,0:80],axis = 1)#/np.shape(test_DPC)[0]
#mean_VISI = np.zeros(np.shape(data_AMP)[0])
#mean_DCI = np.zeros(np.shape(data_AMP)[0])
#shift_max = np.zeros(np.shape(data_AMP)[0])
#shift_min = np.zeros(np.shape(data_AMP)[0])
#phase_shift_symm_G0_long =np.zeros((np.shape(data_AMP)[0],250))
#phase_shift_symm_noG0_long =np.zeros((np.shape(data_AMP)[0],250))
#source_sampledist =np.array([0.253,0.353,0.453,0.553,0.653,0.753])
#p1 = 5*10**-6
#wavelength = (12.398/45)*10**-10
##### scale shifts with ratio between go<-->g1 and source<-->g1, respectively
####distances with go
#if g0:
#    sourceg0 = 0.086
#    g0g1= 0.925
#    g0g1g2 = g0g1
#    g0g2_sampledist = 2*g0g1+sourceg0-source_sampledist
#    DCI_scalefactor = (wavelength/p1)*(2*g0g1-g0g2_sampledist)
#    phase_ratio = (source_sampledist-sourceg0)/g0g1
#else:
#    ### distances no g0
#    sourceg1 = 0.975
#    g1g2 = 0.967
#    g2_sampledist = 2*sourceg1 -source_sampledist
#    DCI_scalefactor = (wavelength/p1)*(sourceg1+g1g2-g2_sampledist)*(g1g2/sourceg1)
#    phase_ratio = source_sampledist/sourceg1 
#
#plt.figure('phase-'+data_folder[:len(data_folder)-2])
#for i in range(6):
#    mean_VISI[i] = np.mean(original_data_VISI[i,200:550,200:550]) 
#    #lineplot = np.mean(cut_corrig_DPC[i,300:550,y_start[i]:y_end[i]], axis = 1)
#    lineplot = np.mean(original_data_DPC[i,350:510,roi_list_DPC[i,0]:roi_list_DPC[i,2]], axis = 1)
#    shift_min[i] = np.min(lineplot)
#    shift_max[i] = np.max(lineplot)
#    mean_DCI[i] = np.mean(original_data_DCI[i,roi_list_DCI[i,1]:roi_list_DCI[i,2],roi_list_DCI[i,0]:roi_list_DCI[i,2]], axis = (0,1)) 
##    if g0:                                
##        phase_shift_symm_G0_long[i]= np.mean(cut_corrig_DPC[i,300:550,100:180], axis = 1)
##    else:
##        phase_shift_symm_noG0_long[i]= np.mean(cut_corrig_DPC[i,300:550,y_start[i]:y_end[i]], axis = 1)
#    plt.plot(lineplot,label = 'position z'+str(i+1))
#plt.legend(loc = 0)
#
##lineplot = np.mean(cut_corrig_DPC[5,300:550,100:180], axis = 1)
##plt.figure('test')
##plt.plot(lineplot)        
#    
##plt.figure('uncorrig')
##plt.plot(lineplot_test)
##test_DPC = data_DPC[396:492,59:700]### slicing for watertube_wood_SiO2_05_z1_info.tif
##lineplot = np.mean(test_DPC[34:,0:500],axis = 1)#,np.mean(test_DPC[34:,0:80],axis = 1),np.mean(test_DPC[34:,0:80],axis = 1)#/np.shape(test_DPC)[0]
##lineplot1  = test_DPC[:,40:41]
##plt.figure('test')
##plt.plot(phase_shift_symm_noG0_long[5])
##plt.plot(phase_shift_symm_G0_long[5])
#    
#    
#if saveon:
#    output_data = np.ones((5,np.shape(data_AMP)[0]))
#    for y in range(np.shape(data_AMP)[0]):
#        output_data[0,y] = phase_ratio[y]
#        output_data[1,y] = shift_min[y]
#        output_data[2,y] = DCI_scalefactor[y]        
#        output_data[3,y] = mean_DCI[y]
#        output_data[4,y] = mean_VISI[y]
#
#    output_data = np.swapaxes(output_data,0,1)
#    try:
#        os.makedirs(filepath)    
#    except:
#        print('Folder already exist')
#    try:
#       # np.savetxt(filepath+filename,((watts,mean_spot_size_vertical,mean_spot_size_horizontal,mean_spot_size_diagonal)),fmt = '%5.5f',delimiter = ',',header = info)
#        np.savetxt(filepath+filename, output_data,fmt = '%10.9f',delimiter = ',',header = info)
#        np.save(filepath+filename[:len(filename)-4]+'.npy',output_data)
#        
#        
#    except:
#        print('Error, could not save the file (does the script/terminal have write access?)')





