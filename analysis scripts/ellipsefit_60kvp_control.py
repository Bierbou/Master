# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 11:41:15 2015

@author: ga56pan
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.collections import LineCollection
#from matplotlib.lines import Line2D
import matplotlib.cm as cm

import os 

import sys
sys.path.append('/users/Baier/Master/analysis_scripts')
import Ellipsefitclass as Ell
plt.close('all')
#============================================================================================================
#                        Data input
#============================================================================================================
data = np.loadtxt('/users/Baier/analysis_files/ellipse_xycoords_60kvp_control.txt', delimiter = ', ')
data = np.swapaxes(data,0,1)    
data = np.delete(data,0,1)
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/ellipse_fits/'
filename = 'ellipse_data_points_60kvp_control.txt'
filename1 = 'ellipse_data_points_60kvp_control.npy'
# write the header for txt file
info = 'contains at first column the different Watts, then 2x2 vertical xy, then 2x2 horizontal xy, then 2x2 diagonal xy value colums,respectively. The left one of the two is always the x direction '
#============================================================================================================
#                        func def
#============================================================================================================
#def ellipse(ra,rb,ang,x0,y0,Nb=50):
#    """ra - major axis length
#    rb - minor axis length
#    ang - angle
#    x0,y0 - position of centre of ellipse
#    Nb - No. of points that make an ellipse
#
#    based on matlab code ellipse.m written by D.G. Long,
#    Brigham Young University, based on the
#    CIRCLES.m original
#    written by Peter Blattner, Institute of Microtechnology,
#    University of
#    Neuchatel, Switzerland, blattner@imt.unine.ch
#    """
#    xpos,ypos=x0,y0
#    radm,radn=ra,rb
#    an=ang
#
#    the=linspace(0,2*np.pi,Nb)
#    X=radm*np.cos(the)*np.cos()-np.sin*radn*sin(the)+xpos
#    Y=radm*cos(the)*si+co*radn*sin(the)+ypos
#    return X,Y
#============================================================================================================
#                        Data handling
#============================================================================================================
watts = data[0]
x_val = np.ones((np.shape(watts)[0],6))
y_val =np.ones((np.shape(watts)[0],6))
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
alpha_shift = 3.05
## --> different phis 
phi_v = alpha_shift + 90
phi_h = alpha_shift 
phi_d = alpha_shift + 90 +45 

arc = 2.
R = np.arange(0,arc*np.pi, 0.01)
length = np.shape(R)[0]
xx = np.ones((np.shape(watts)[0],length))
yy = np.ones((np.shape(watts)[0],length))
for i in range(np.shape(watts)[0]-1):#np.shape(watts)[0]
    #plt.figure('ell_fit')
    try:
        fit_object = Ell.Ellipsefit()
        #1/0
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
    
        1/0
        #x_max_e = (2*np.max(x_val))
        #y_max_e = (2*np.max(y_val))
        #e = [ -x_max_e, x_max_e, -y_max_e, y_max_e]
        #plt.axis(e)
        #plt.plot(x_val[i],y_val[i], 'b+',lw = 1.)
        #plt.plot(xx[i],yy[i], color = 'red')
        #1/0
    except:
#        xx[i,:] = (center[0] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R) - np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R))
#        yy[i,:] = (center[1] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R) + np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R))
        #xx[i,:] = (center[0] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.cos(-phi_h) - np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.sin(-phi_h))
        #yy[i,:] = (center[1] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.sin(-phi_h) + np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.cos(-phi_h))
        xx[i,:] = (center[0] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.cos(0) - np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.sin(0))
        yy[i,:] = (center[1] + np.sqrt((x_val[i,2]**2.)+(y_val[i,2]**2.))*np.cos(R)*np.sin(0) + np.sqrt((x_val[i,0]**2.)+(y_val[i,0]**2.))*np.sin(R)*np.cos(0))
        
        ## without rotation
        
        #plt.plot(x_val[i],y_val[i], 'b+',lw = 1.)
        #plt.plot(xx[i],yy[i], color = 'red')
        ++i


lines_cut = [[xx[0], yy[0]],[xx[1], yy[1]],[xx[2], yy[2]],
         [xx[4], yy[4]],[xx[9], yy[9]],[xx[14], yy[14]],
         [xx[18], yy[18]],[xx[21], yy[21]]]
        
        
lines3D = [[xx[0], yy[0],watts[0]],[xx[1], yy[1],watts[1]],[xx[2], yy[2],watts[2]],
         [xx[3], yy[3],watts[3]],[xx[4], yy[4],watts[4]],[xx[5], yy[5],watts[5]],
         [xx[6], yy[6],watts[6]],[xx[7], yy[7],watts[7]],[xx[8], yy[8],watts[8]],
         [xx[9], yy[9],watts[9]],[xx[10], yy[10],watts[10]],[xx[11], yy[11],watts[11]],
         [xx[12], yy[12],watts[12]],[xx[13], yy[13],watts[13]],[xx[14], yy[14],watts[14]],
         [xx[15], yy[15],watts[15]],[xx[16], yy[16],watts[16]],[xx[17], yy[17],watts[17]],
         [xx[18], yy[18],watts[18]],[xx[19], yy[19],watts[19]],[xx[20], yy[20],watts[20]],
         [xx[21], yy[21],watts[21]],[xx[22], yy[22],watts[22]],[xx[23], yy[23],watts[23]],[xx[24], yy[24],watts[24]]]

# Reformat it to what `LineCollection` expects:
#lines = [zip(x, y) for x, y in lines]
w = np.array(np.arange(0,25,1.))
z = np.array(np.arange(0,8,1.))
def normalize(z):
    z = z.copy()
    z -= z.min()
    z /= z.max()
    return z


powers =['30 W','35 W','40W','50W','75W','100W','125W','135W']

cmap = plt.get_cmap('autumn_r')
plt.figure('ell_fit 60 kvp control')
for (x, y), color, label in zip(lines_cut, normalize(z), powers):
    plt.plot(x, y, label=label, color=cmap(color), lw=2)
m = cm.ScalarMappable(cmap=cmap)
m.set_array(z)
plt.plot(x_val[0],y_val[0], 'b+', label = 'fit points')
plt.plot(x_val[1],y_val[1], 'b+')
plt.plot(x_val[2],y_val[2], 'b+')
plt.plot(x_val[4],y_val[4], 'b+')
plt.plot(x_val[9],y_val[9], 'b+')
plt.plot(x_val[14],y_val[14], 'b+')
plt.plot(x_val[18],y_val[18], 'b+')
plt.plot(x_val[21],y_val[21], 'b+')
#plt.colorbar(m)
plt.xticks(np.arange(-180,181,30),fontsize = 11)
plt.yticks(np.arange(-90,91,10),fontsize = 11)
plt.xlim(-190,190)
plt.ylim(-95,95)
plt.xlabel('x-axis',fontsize = 12)
plt.ylabel('y-axis',fontsize = 12)
plt.grid(color = 'k',ls = ':')
plt.legend(fontsize = 10,title = 'Power',loc = 0)
plt.savefig(filepath_pic+'ellipse_fit_60kvp_control_2d_selection_phi0.pdf', format = 'pdf',dpi = 300 )


Axes3D(plt.figure('3d 60kvp control'))
for (x, y,z), color, label in zip(lines3D, normalize(w), watts):
    plt.plot(x, y,z, label=label, color=cmap(color), lw=2)
m = cm.ScalarMappable(cmap=cmap)
m.set_array(w)

plt.xlabel('x-axis',fontsize = 12)
plt.ylabel('y-axis',fontsize = 12)

#Axes3D(plt.figure('plot'))
#for i in range(np.shape(watts)[0]):
#    plt.plot(xx[i],yy[i],watts[i],lw = 10)
#Axes3D(plt.figure())
#for i in range(np.shape(watts)[0]):
#    plt.contour(xx[i],yy[i]),watts[i])
    
#fig = plt.figure('surface')
#ax = plt.axes(projection='3d')
#for i in range(np.shape(watts)[0]):
#    ax.plot_surface(xx[i], yy[i], watts[i], cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=10)




