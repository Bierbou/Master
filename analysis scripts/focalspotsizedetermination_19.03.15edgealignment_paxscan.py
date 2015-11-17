# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:03:30 2015

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import scipy.special as scs
import scipy.optimize as opti
sys.path.append('/users/schaff/Python Scripts/')
#import pyFlorian as pyF
"""
Script determing the witdth of focal size and shape
with line-spread-function analysis
"""
plt.close("all")
#Picture
pic_numb = 20
"""
Function definitions
"""
#errorfunction fiting the lineedge
def error_fit_func(x,a,b,c,d): 
    return (a*scs.erf((x-b)/(math.sqrt(2)*c))+d)
#gaussfunction determin the FWHM
def gauss (x,a,b,c,d,):
    return (a*np.exp(-0.5*((x-b)/c)**2) + d)
    
"""
data reading
find angle for which the slope is max 
"""
# readin flatfield wihtout a sample and correct vertical black detector edge by linear regression
#for paxscan
flatfield1 = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_control/paxscan/ct/paximage_ct_479459.h5")["raw_data"]
#flatfield2 = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_measurement_paxscan_60kvp/horizontal_edge/paxscan/ct/paximage_ct_077369.h5")["raw_data"]
#flatfield3 = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_measurement_paxscan_60kvp/diagonal_edge/paxscan/ct/paximage_ct_077373.h5")["raw_data"]
#1/0



flatfield = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_control2/paxscan/ct/paximage_ct_479479.h5")["raw_data"]
flatfield = flatfield*1.0
flatfield = np.nan_to_num(flatfield)


# data with edge
data_norm = np.zeros((pic_numb,300,700))
linesum = np.zeros((pic_numb,700), dtype = "float64")
pixel_number = np.arange(0,700)
fwhm = np.zeros((pic_numb), dtype = "float64")
"""
Sum over all 440 horizontal lines 
and write it in a (41,440) array
""" 
# for vertical edge set axis = 0
# for horizontal edge set axis = 1
direction = 0
for i in range(pic_numb):
    filename = "/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_control2/paxscan/ct/paximage_ct_%i.h5"%(479480+i)
    print "Loading file", filename
    # for paxscan
    img = e17.io.h5read(filename)["raw_data"]
    #1/0
    img = img[400:700,0:700]*1.0

    data_norm[i] = img /flatfield[400:700,0:700] #divide data by flatfield 
    linesum[i] = np.mean(data_norm[i],axis = direction)
    linesum[i] = linesum[i][::-1]
    #linesum[i] = np.sum(data_norm[i],axis = direction)
"""
Evaluate the best fitparameter for given startvalues
"""
best_fitparam = np.zeros((pic_numb,4), dtype = "float64")  
# strongly depending on startparametervalue for c 


pixel_number  = pixel_number*1.
max_pixel_number = max(pixel_number)
pixel_number = pixel_number/(max_pixel_number*1.) 
start_params = [-1.,0.5,0.5,1.]
#1/0
for i in range(pic_numb):
    #opti.curve_fit(error_fit_func,pixel_number,linesum[0],p0 = params)
    best_fitparam[i],pcov = opti.curve_fit(error_fit_func,pixel_number,linesum[i,:],p0 = start_params)
##### problems with last 10 angles start parameter do not converge
#start_params = [1.,10.,350.,10.]
#popt,pcov = opti.curve_fit(error_fit_func,pixel_number,linesum[32],p0 = start_params)
#print popt
#print pcov
#plt.figure(3)
#plt.plot(pixel_number, error_fit_func(pixel_number,83.32279932,  81111.52052572,  -8290.19821163 ,   148.32969235))  
#####
pixel_number = pixel_number*max_pixel_number
best_fitparam*=max_pixel_number 
"""
Take fitparams and give it to gaussfit and determine FWHM
""" 
plt.figure(1)
for i in range(pic_numb):
    fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[i,2])
    plt.plot(pixel_number, gauss(pixel_number,best_fitparam[i,0],best_fitparam[i,1], best_fitparam[i,2], best_fitparam[i,3]))
plt.set_cmap('autumn')
print "Position of minimum ==> Value of minimum" 
print "       ",np.argmin(fwhm)+1 ,"     ==>     ", np.min(fwhm)
print "       ","right angle is at postion:" ,"         ==>     ", "curphi-startphi+", np.argmin(fwhm)+1
"""
check the tilt of the edge to the detector lines, 
comparing different picture lines and different angles
"""
#plt.figure(4)
#plt.plot(pixel_number, data_norm[0,0,:])
#plt.figure(4)
#plt.plot(pixel_number, data_norm[0,129,:])
#plt.figure(5)
#plt.plot(pixel_number, data_norm[40,0,:])
#plt.figure(5)
#plt.plot(pixel_number, data_norm[40,129,:])
"""
plot results for errorfitfunction
"""
plt.figure(2) #plot all 41 angle steps in one plot
for i in range(pic_numb):
    plt.plot(pixel_number, linesum[i]) 
plt.figure(3) #plot all 41 fitfunctions in one plot
for i in range(pic_numb):
    plt.plot(pixel_number, error_fit_func(pixel_number,best_fitparam[i,0], best_fitparam[i,1], best_fitparam[i,2], best_fitparam[i,3]), )

plt.figure('fwhm')
plt.plot(fwhm)


