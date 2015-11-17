# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:14:36 2015

@author: ga56pan
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.collections import LineCollection
#from matplotlib.lines import Line2D

plt.close('all')

def fit_linear(y_data, x_data=None,plot=False):
    """
    fits linear ax + b to data    
    """
    
    if x_data == None:
        x_data = np.arange(y_data.shape[0])
        
    A = np.array([np.ones(x_data.shape[0]), x_data])
    x = np.linalg.lstsq(A.T,y_data)[0] # obtaining the parameters
    
    if plot == True:
        # plotting the line
        line = x[0] + x[1]*A[1] 
        #plt.plot(A[1],line,'k-')#,A[1],y_data,'o')
        #plt.axis([0,1,0,1])
        print(("%s + %s * x" %(np.round(x[0],3),np.round(x[1],3))))
        #plt.draw()
    
    return x,line        
#============================================================================================================
#                        Data input
#============================================================================================================
data_g0 = np.load('/users/Baier/analysis_files/watertube_wood_SiO2_05.npy')
data_nog0 = np.load('/users/Baier/analysis_files/noG0_watertube_wood_SiO2_05.npy')    
#data_g0 = np.loadtxt('/users/Baier/analysis_files/watertube_wood_SiO2_05.txt', delimiter = ', ')
#data_g0 = np.swapaxes(data_g0,0,1)
#data_nog0 = np.loadtxt('/users/Baier/analysis_files/noG0_watertube_wood_SiO2_05.txt', delimiter = ', ')
#data_nog0 = np.swapaxes(data_nog0,0,1)
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/setupcomp/'
filename = 'ellipse_data_points_60kvp_control.txt'
filename1 = 'ellipse_data_points_60kvp_control.npy'
# write the header for txt file
info = 'agadfg '
#============================================================================================================
#                        Data handling
#============================================================================================================

## data for setup with g0
scale_factor_phase_g0 = data_g0[:,0]
phase_shift_g0 = data_g0[:,1]*-1.
DCI_scale_g0 = data_g0[:,2]*1000000
DCI_g0 = data_g0[:,3]
visibility_g0 = data_g0[:,4]

### data for setup without g0
scale_factor_phase_nog0 = data_nog0[:,0]
phase_shift_nog0 = data_nog0[:,1]*-1.
DCI_scale_nog0 = data_nog0[:,2]*1000000
DCI_nog0 = data_nog0[:,3]
visibility_n0g0 = data_nog0[:,4]


### linear fits 
xg0,ling0 = fit_linear(phase_shift_g0, x_data =scale_factor_phase_g0,plot=True)
xnog0,linnog0 =fit_linear(phase_shift_nog0, x_data =scale_factor_phase_nog0,plot=True)
x = np.asanyarray(np.arange(.10,.9,0.1)) 
linfit_g0 = xg0[0]+ xg0[1]*x
linfit_nog0 = xnog0[0]+ xnog0[1]*x

plt.figure('scaled phase')
plt.plot(scale_factor_phase_g0,phase_shift_g0,'bx',markersize = 10,label = 'values with G0')
plt.plot(x,linfit_g0,'b--',label = 'linear fit: %s + %s * x' %(np.round(xg0[0],3),np.round(xg0[1],3)))
plt.plot(scale_factor_phase_nog0,phase_shift_nog0,'rx',markersize = 10, label = 'values without G0')
plt.plot(x,linfit_nog0,'r--',label = 'linear fit: %s + %s * x' %(np.round(xnog0[0],3),np.round(xnog0[1],3)))
plt.xticks(np.arange(0,1.1,.10),fontsize = 12)
plt.yticks(np.arange(0,1.1,.10),fontsize = 12)
plt.xlim(0.,1.)
plt.ylim(0,1.)
plt.xlabel('Relative sample position',fontsize = 12)
plt.ylabel('Pahse-shift',fontsize = 12)
plt.legend(loc = 0, title = 'Induced Phase-shift',fontsize = 12 )
plt.savefig(filepath_pic+'phase_shifts_standard_setup_nog0_g0.pdf', format = 'pdf',dpi = 300 )


plt.figure('scaled DCI')
plt.plot(DCI_scale_g0,DCI_g0,'b+',markersize = 10,label = 'DCI with G0')
plt.plot(DCI_scale_g0,DCI_g0,'b--',)
plt.plot(DCI_scale_nog0,DCI_nog0,'r+',markersize = 10, label = 'DCI without G0')
plt.plot(DCI_scale_nog0,DCI_nog0,'r--')
#plt.xticks(np.arange(0,1.1,.10),fontsize = 12)
plt.yticks(np.arange(0,1.1,.10),fontsize = 12)
#plt.xlim(0.,1.)
plt.ylim(0,1.)
plt.xlabel('Correlation length $\zeta$ [um]',fontsize = 12)
plt.ylabel('DCI-singal',fontsize = 12)
plt.legend(loc = 0, title = 'DCI-signal',fontsize = 12 )
plt.savefig(filepath_pic+'DCI_standard_setup_nog0_g0.pdf', format = 'pdf',dpi = 300 )
