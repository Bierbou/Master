# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 15:35:45 2014

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import scipy.special as scs
import scipy.optimize as opti
import pyfits as fit
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
#for pilatus data 
flatfield,meta =  e17.io.cbfread("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_fine/pilatus/ct/specadm_1_ct_1946493.cbf")
flatfield[:,242] = flatfield[:,241]+ (flatfield[:,245]-flatfield[:,241])/3
flatfield[:,243] = flatfield[:,241]+ 2*(flatfield[:,245]-flatfield[:,241])/3
flatfield[:,244] = flatfield[:,241]+ (flatfield[:,245]-flatfield[:,241])
#for ccd data
#flatfielt = fit.getdata("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_fine/pilatus/ct/specadm_1_ct_1946493.cbf")
flatfield = flatfield*1.0
# data with edge
data_norm = np.zeros((pic_numb,100,400))
data = np.zeros((pic_numb,100,400))
linesum = np.zeros((pic_numb,400), dtype = "float64")
pixel_number = np.asarray(range(50,450))
fwhm = np.zeros((pic_numb), dtype = "float64")
"""
Sum over all 440 horizontal lines 
and write it in a (41,440) array
""" 
# for vertical edge set axis = 0
# for horizontal edge set axis = 1
direction = 0
for i in range(pic_numb):
    print "Loading file", (1945978+i)
    # for pilatus data 
    img,meta = e17.io.cbfread("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_fine/pilatus/ct/specadm_1_ct_%i.cbf"%(1946494+i))
    img[:,242]=img[:,241]+(img[:,245]-img[:,241])/3
    img[:,243]=img[:,241]+2*(img[:,245]-img[:,241])/3
    img[:,244]=img[:,241]+(img[:,245]-img[:,241]) 
    #for ccd data 
    #img = fit.getdata("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_fine/pilatus/ct/specadm_1_ct_%i.cbf"%(1946494+i))
    img = img[200:300,50:450]*1.0

    data_norm[i] = img /flatfield[200:300,50:450] #divide data by flatfield 
    linesum[i] = np.mean(data_norm[i],axis = direction)
#    linesum[i] = data_norm[i,20]
"""
Evaluate the best fitparameter for given startvalues
"""
best_fitparam = np.zeros((pic_numb,4), dtype = "float64")  
# strongly depending on startparametervalue for c  
start_params = [1.,220.,0.7,1.]
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
"""
Take fitparams and give it to gaussfit and determine FWHM
""" 
plt.figure(1)
for i in range(pic_numb):
    fwhm[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[i,2]) 
    plt.plot(pixel_number, gauss(pixel_number,best_fitparam[i,0],best_fitparam[i,1], best_fitparam[i,2], best_fitparam[i,3]))
plt.set_cmap('autumn')
print "Position of minimum ==> Value of minimum" 
print "       ",np.argmin(fwhm)+1 ,"         ==>     ", np.min(fwhm)
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
