# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 15:24:16 2015

@author: ga56pan
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os 
def rect(r, theta):
    """theta in degrees

    returns tuple; (float, float); (x,y)
    """
    x = r * math.cos(math.radians(theta))
    y = r * math.sin(math.radians(theta))
    return x,y
#============================================================================================================
#                        Data input
#============================================================================================================
##### Mean spotsizes for vertical,horizontal and diagonal @ 60 kvp #########
f = np.loadtxt('/users/Baier/analysis_files/linespread_acquisition_flat_measurement_paxscan_60kvp.txt', delimiter = ', ')
f = np.swapaxes(f,0,1)    

#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/ellipse_fits/'
try:
    os.makedirs(filepath_pic)    
except:
    print('Folder already exist')
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
filename = 'ellipse_xycoords_60kvp.txt'
filename1 = 'ellipse_xycoords_60kvp.npy'
# write the header for txt file
info = 'contains at first column the different Watts, then 2x2 vertical xy, then 2x2 horizontal xy, then 2x2 diagonal xy value colums,respectively. The left one of the two is always the x direction '
#============================================================================================================
#                        Data handling
#============================================================================================================
'''
Be carefull with the theta because the tilt of the 'right' vertical and horizontal ellipse axis of the resulting spot
'''
alpha_shift = 3.05
## --> different phis 
phi_v = alpha_shift + 90
phi_h = alpha_shift 
phi_d = alpha_shift + 90 +45 
   
watts = f[0]
watts = np.insert(watts ,0,0)
vertical = f[1]
vertical = np.insert(vertical,0,0)
horizontal = f[2]
horizontal = np.insert(horizontal,0,0)
diagonal = f[3]
diagonal = np.insert(diagonal,0,0)
a = watts.shape[0]

""" vertical data points in cartesian coordinates """
### vertical points at phi_v
vert_phi_v1 = np.ones((2,a))
vert_phi_v1[0,:],vert_phi_v1[1,:]= rect(vertical/2., phi_v)
### vertical points at phi_v +180
vert_phi_v2 = -vert_phi_v1
""" horizontal data points in cartesian coordinates """
##horizontal points at phi_h 
hor_phi_h1 = np.ones((2,a))
hor_phi_h1[0,:],hor_phi_h1[1,:]= rect(horizontal/2., phi_h)
##horizontal points at phi_h +180
hor_phi_h2 = -hor_phi_h1
""" diagonal data points in cartesian coordinates """
## diagonal points at phi_d 
diag_phi_d1 = np.ones((2,a))
diag_phi_d1[0,:],diag_phi_d1[1,:]= rect(diagonal/2., phi_d)
## diagonal points at phi_d+ 180 
diag_phi_d2 = -diag_phi_d1


#============================================================================================================
#                        Plots
#============================================================================================================
plt.figure('60kvp cartes coords',figsize = (8,7))
plt.grid(color = 'k',ls = ':',axis = 'both')
plt.plot(vert_phi_v1[0], vert_phi_v1[1],'r+',label ='vertical')#, lw = 2.,ls = ':')
plt.plot(vert_phi_v2[0], vert_phi_v2[1],'r+')#, lw = 2.,ls = '.')
plt.plot(hor_phi_h1[0], hor_phi_h1[1],'b+',label = 'horizontal')#, ls = '.')
plt.plot(hor_phi_h2[0], hor_phi_h2[1],'b+')#, ls = '.')
plt.plot(diag_phi_d1[0], diag_phi_d1[1],'g+',label = 'diagonal')#, ls = '.')
plt.plot(diag_phi_d2[0], diag_phi_d2[1],'g+')#, ls = '.')
plt.title('spot size in x-y coords 60kvp')
plt.xlim((-130,130))
plt.xticks( np.arange(np.round(1.007*np.min(hor_phi_h2)), np.round(1.007*np.max(hor_phi_h1))+1,25) ,fontsize  = 12)
plt.yticks( np.arange(np.round(1.083*np.min(vert_phi_v2)),np.round(1.083*np.max(vert_phi_v1)+1),10) ,fontsize  = 12)
plt.ylabel(' x-axis [$\mu$m]',fontsize = 12)
plt.xlabel('x-axis [$\mu$m]',fontsize = 12)
plt.legend( loc = 0,fontsize = 10,title = 'Spot-sizes')
plt.savefig(filepath_pic+'Results xy coords 60kvp.pdf', format = 'pdf',dpi = 300)

#============================================================================================================
#                        data saving
#============================================================================================================
output_data = np.ones((13,a))
output_data[0,:] = watts
output_data[1,:] = vert_phi_v1[0]
output_data[2,:] = vert_phi_v1[1]
output_data[3,:] = vert_phi_v2[0]
output_data[4,:] = vert_phi_v2[1]
output_data[5,:] = hor_phi_h1[0]
output_data[6,:] = hor_phi_h1[1]
output_data[7,:] = hor_phi_h2[0]
output_data[8,:] = hor_phi_h2[1]
output_data[9,:] = diag_phi_d1[0]
output_data[10,:] = diag_phi_d1[1]
output_data[11,:] = diag_phi_d2[0]
output_data[12,:] = diag_phi_d2[1]
try:
    os.makedirs(filepath)    
except:
    print('Folder already exist')
try:
   # np.savetxt(filepath+filename,((watts,mean_spot_size_vertical,mean_spot_size_horizontal,mean_spot_size_diagonal)),fmt = '%5.5f',delimiter = ',',header = info)
    np.save(filepath+filename1,output_data)
    output_data = np.swapaxes(output_data,0,1)
    np.savetxt(filepath+filename, output_data,fmt = '%10.5f',delimiter = ', ',header = info)
    
except:
    print('Error, could not save the file (does the script/terminal have write access?)')





























