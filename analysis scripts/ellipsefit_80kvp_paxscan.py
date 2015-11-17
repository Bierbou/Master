# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 14:24:09 2015

@author: ga56pan
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os 

import sys
sys.path.append('/users/Baier/Master/analysis_scripts')
import Ellipsefitclass as Ell
plt.close('all')
#============================================================================================================
#                        Data input
#============================================================================================================
data = np.loadtxt('/users/Baier/analysis_files/ellipse_xycoords_80kvp_paxscan.txt', delimiter = ', ')
data = np.swapaxes(data,0,1)
data = np.delete(data,0,1)    
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/80kvp_spot_sizes/'
filename = 'ellipse_data_points_80kvp_paxscan.txt'
filename1 = 'ellipse_data_points_80kvp_paxscan.npy'
# write the header for txt file
info = 'contains at first column the different Watts, then 2x2 vertical xy, then 2x2 horizontal xy, then 2x2 diagonal xy value colums,respectively. The left one of the two is always the x direction '

#============================================================================================================
#                        Data handling
#============================================================================================================
watts = data[0]
x_val = np.ones((np.shape(watts)[0],6))
y_val  =np.ones((np.shape(watts)[0],6))
for j in range(6):
    for i in range(np.shape(watts)[0]):
        x_val[i,j] = data[1+j*2,i]
        y_val[i,j] = data[2+j*2,i]
watts = watts[5:]
x_val = x_val[5:,:]
y_val = y_val[5:,:]
#============================================================================================================
#                        Main 
#============================================================================================================
alpha_shift = 3.1564
## --> different phis 
phi_v = alpha_shift + 90
phi_h = alpha_shift 
phi_d = alpha_shift + 90 +45 

arc = 2.
R = np.arange(0,arc*np.pi, 0.01)
length = np.shape(R)[0]
xx = np.ones((np.shape(watts)[0],length))
yy = np.ones((np.shape(watts)[0],length))
for i in range(np.shape(watts)[0]):#np.shape(watts)[0]
    plt.figure('ell_fit')
    try:
        fit_object = Ell.Ellipsefit()
        a = fit_object.fit_ellipse(x_val[i],y_val[i])
        center = fit_object.ellipse_center(a)
        phi = fit_object.ellipse_angle_of_rotation2(a)
        axes = fit_object.ellipse_semi_axis_length(a)
        print ("watts = ", watts[i])
        print("center = ",  center)
        print("angle of rotation = ",  phi)
        print("axes = ", axes)
        
        
        half_a, half_b = axes
        xx[i,:] = np.nan_to_num(center[0] + half_a*np.cos(R)*np.cos(phi) - half_b*np.sin(R)*np.sin(phi))
        yy[i,:] = np.nan_to_num(center[1] + half_a*np.cos(R)*np.sin(phi) + half_b*np.sin(R)*np.cos(phi))
        
        
        x_max_e = (1.2*np.max(x_val))
        y_max_e = (1.2*np.max(y_val))
        e = [ -x_max_e, x_max_e, -y_max_e, y_max_e]
        plt.axis(e)
        plt.plot(x_val[i],y_val[i], 'b+',lw = 1.)
        plt.plot(xx[i],yy[i], color = 'red')
        #1/0
    except:
        xx[i,:] = (center[0] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.cos(phi_h) - np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.sin(phi_h))
        yy[i,:] = (center[1] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.sin(phi_h) + np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.cos(phi_h))
        plt.plot(x_val[i],y_val[i], 'b+',lw = 1.)
        plt.plot(xx[i],yy[i], color = 'red')
        ++i

Axes3D(plt.figure())
for i in range(np.shape(watts)[0]):
    plt.plot(xx[i],yy[i],watts[i],lw = 10)
#Axes3D(plt.figure())
#for i in range(np.shape(watts)[0]):
#    plt.contour(xx[i],yy[i])#,watts[i])
#    
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#for i in range(np.shape(watts)[0]):
#    ax.plot_surface(xx[i], yy[i], watts[i], cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=1)

        
        
        
        



