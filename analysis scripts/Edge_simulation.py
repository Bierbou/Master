# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 12:34:59 2014

@author: ga56pan
"""
import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import scipy.special as scs
import scipy.optimize as opti
import scipy.ndimage as ndi
import pyCT as ct
sys.path.append('/users/schaff/Python Scripts/')
#import pyFlorian as pyF
"""
Script determing the witdth of focal size and shape
with line-spread-function analysis
"""
plt.close("all")

#errorfunction fiting the lineedge
def error_fit_func(x,a,b,c,d): 
    return (a*scs.erf((x-b)/(math.sqrt(2)*c))+d)
#gaussfunction determin the FWHM
def gauss (x,a,b,c,d,):
    return (a*np.exp(-0.5*((x-b)/c)**2) + d)
#============================================================================================================
# Analysis parameter
#============================================================================================================
vertical_edge = True
horizontal_edge = False
diagonal_edge = False
detector_pixel = 800
pix_size = 20
#detector_diag = detector_pixel*np.sqrt(2)
well_size = 100.
angle_tilt = -3. ## be carefull projector needs reflected angle values
pixel_subdivision = 0.01 # for oversampling 
edge_dist = 8.
setup_len = 195.6
magnification = setup_len/edge_dist
cutoff = int((detector_pixel- detector_pixel/np.sqrt(2) )/2+.5)*int(1/pixel_subdivision)
#fwhm_vert_pix[i] = np.abs(2*(2*np.log(2))**.5*best_fitparam[2])
#======================
# parameter for pahntom
#======================
a_edge = -.5 # Amplitude 
b_edge = detector_pixel/2 # offset in x-direction
c_edge = (well_size*(magnification-1))/(pix_size*2*(2*np.log(2))**.5) # width of the edge
d_edge = 0.5 # offset in y-direction
start_params_edge = [a_edge,b_edge,c_edge,d_edge] # Startparameter for edgegeneration
#======================
# parameter for fit
#======================
a = 1. 
b = 0.5
c = 0.5#(pix_size*2*(2*np.log(2))**.5)*(magnification-1)# /(magnification-1.)
d = 0.
start_params = [a,b,c,d] # Startparameter for errorfit
#============================================================================================================
# Phantom generation
#============================================================================================================
x_val = np.arange(detector_pixel)
test_esf = error_fit_func(x_val, *start_params_edge)
test_edge = np.tile(test_esf,(detector_pixel,1))
#============================================================================================================
# Analysis 
#============================================================================================================
rot_edge = ndi.rotate(test_edge, angle_tilt, reshape = False , mode = 'nearest' )
if vertical_edge:
    angle = angle_tilt+ 90.
if horizontal_edge:
    angle = angle_tilt
if diagonal_edge:
    angle = angle_tilt + 45.
P = ct.P.Parallel_AWP()
P.angles = angle
P.is_rescale = True #for subdivision

P.proj_area_size = (detector_pixel*int(1/pixel_subdivision), 1)
norm = P.project(np.ones((detector_pixel,detector_pixel)))[0]/np.sqrt(2)
norm = norm[cutoff :-cutoff]
proj = P.project(rot_edge)[0]/np.sqrt(2)
proj = proj[cutoff :-cutoff]
proj_norm = np.nan_to_num(proj/norm)
#proj = np.nan_to_num(proj)

projected_pixel_size = (pix_size*np.sqrt(2)) * pixel_subdivision
x_values = 1. * np.arange(detector_pixel*int(1/pixel_subdivision))*projected_pixel_size

max_x_values = max(x_values)
x_values = x_values/max_x_values


#runtime test
best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_values[cutoff :-cutoff], proj_norm, p0=start_params)
#normal code
#best_fitparam, cov_esf = opti.curve_fit(error_fit_func, x_values[cutoff :-cutoff], proj_norm[cutoff :-cutoff], p0=start_params)
esf_fit = error_fit_func(x_values, *best_fitparam)

x_values = x_values*max_x_values
best_fitparam*= max_x_values
 
fwhm_vert = (np.abs(2*(2*np.log(2))**.5*best_fitparam[2]) )/(magnification-1)

plt.figure('Results Edge+ Fit')
plt.plot(x_values[cutoff :-cutoff],proj_norm)
plt.plot(x_values[cutoff :-cutoff],esf_fit[cutoff :-cutoff ])
plt.figure('Gauss')
plt.plot(x_values[cutoff :-cutoff ], gauss(x_values[cutoff :-cutoff ],best_fitparam[0],best_fitparam[1], (best_fitparam[2]/(magnification-1)), best_fitparam[3]))

print 'FWHM horizontal ==>  ', fwhm_vert, '[micron]'

plt.axis()





#================================
# Data directory vertical edge
#"/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_measurement2_60kvp/vertical_edge/pilatus/ct/specadm_1_ct_%i.cbf"%(1947170+j*12+i))





