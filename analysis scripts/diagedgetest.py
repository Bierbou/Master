# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 14:32:42 2014

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
sys.path.append('/users/schaff/Python Scripts/')
"""
Script determing the witdth of focal size and shape
with line-spread-function analysis
"""
plt.close("all")
# Picture number to analyzing
pic_numb = 0
"""
Function definitions
"""
#errorfunction fiting the lineedge
def error_fit_func(x,a,b,c,d): 
    return (a*scs.erf((x-b)/(math.sqrt(2)*c))+d)
#gaussfunction determin the FWHM
def gauss (x,a,b,c,d,):
    return (a*np.exp(-0.5*((x-b)/c)**2) + d)
# linear function
def line(x, m, t):
    return m* x + t   
# ******** ANALYSIS PARAMETERS ******** #
pixel_subdivision = .1
start_params = [1.,0.,3.5,1.]
fit_esf = True

"""
data reading
find angle for which the slope is max 
"""
# readin flatfield wihtout a sample and correct vertical black detector edge by linear regression
flatfield,meta =  e17.io.cbfread("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_measurement2_60kvp/vertical_edge/pilatus/ct/specadm_1_ct_1947169.cbf")
flatfield[:,242] = flatfield[:,241]+ (flatfield[:,245]-flatfield[:,241])/3
flatfield[:,243] = flatfield[:,241]+ 2*(flatfield[:,245]-flatfield[:,241])/3
flatfield[:,244] = flatfield[:,241]+ (flatfield[:,245]-flatfield[:,241])
flatfield = flatfield*1.0
# data with edge
img=np.zeros((407,487))
data_norm = np.zeros((pic_numb,100,400))
pixel_number = np.asarray(range(50,450))
fwhm = np.zeros((pic_numb), dtype = "float64")
"""
Sum over all 440 horizontal lines 172
and write it in a (41,440) array
""" 
print "Loading file", (1945808)
img,meta = e17.io.cbfread("/data/DPC/local_setups/microfocus/samples/linespread_acquisition_flat_measurement2_60kvp/vertical_edge/pilatus/ct/specadm_1_ct_%i.cbf"%(1947170))
img[:,242]=img[:,241]+(img[:,245]-img[:,241])/3
img[:,243]=img[:,241]+2*(img[:,245]-img[:,241])/3
img[:,244]=img[:,241]+(img[:,245]-img[:,241])  
img = img[215:405,50:450]*1.0
data_norm = img /flatfield[215:405,50:450] #divide data by flatfield   
data_norm = data_norm[:,65:255]    
# ******** DETERMINE EDGE ANGLE ******** #
img_grad_raw = ct.U.delxf(data_norm, axis=1)
img_grad = img_grad_raw.copy()
#img_grad[img_grad < img_grad.max()*0.5] = 0

edge = (img_grad * np.arange(img_grad.shape[1])).sum(axis = 1) / img_grad.sum(axis = 1) # obtain the weighted edge postion for each image row
p_line, cov_line = opti.curve_fit(line, np.arange(edge.size - 1)[np.isfinite(edge[:-1])], edge[:-1][np.isfinite(edge[:-1])])

angle = np.arctan(p_line[0])/np.pi * 180.
angle_err = np.arctan(1./(p_line[0]**2 + 1.) * np.sqrt(cov_line[0,0]))/np.pi * 180. # Gaussian error propagation with determined variance of the fit

#angle =0.
# ******** RE-PROJECTION OF EDGE / DETERMINATION OF LSF ******** #
P = ct.P.Parallel_AWP()
P.angles = [angle+90.]
P.is_rescale = True
P.proj_area_size = (data_norm.shape[1]*int(1/pixel_subdivision), 1)

esf_full = P.project(data_norm)[0]
#esf_full = esf_full[esf_full > 0]

lsf_full = np.gradient(esf_full)
center = np.argmin(lsf_full)
p = 10

span = p*int(1/pixel_subdivision)

esf = esf_full[center-span:center+span]
x_pix = np.linspace(-p,p,2*span)
#x_pix = np.asarray(range(100,200))

if fit_esf:
    best_fitparam_d, cov_esf = opti.curve_fit(error_fit_func, x_pix, esf, p0=start_params)
    esf_fit = error_fit_func(x_pix, *best_fitparam_d)

plt.figure('Results')
plt.plot(x_pix, esf, lw=2, label = 'Measured ESF')
if fit_esf:
    plt.plot(x_pix, esf_fit, lw=2, label = 'Fitted ESF')
plt.legend()

