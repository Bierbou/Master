# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:16:41 2015

@author: ga56pan
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os 

#============================================================================================================
#                        Data input
#============================================================================================================
##### Mean spotsizes for vertical,horizontal and diagonal @ 60 kvp #########
f = np.loadtxt('/users/Baier/analysis_files/linespread_acquisition_flat_measurement_paxscan_40kvp.txt', delimiter = ', ')
f = np.swapaxes(f,0,1)
g = np.loadtxt('/users/Baier/analysis_files/linespread_acquisition_flat_measurement_control_60kvp.txt', delimiter = ', ')
g = np.swapaxes(g,0,1)    
h = np.loadtxt('/users/Baier/analysis_files/linespread_acquisition_flat_measurement_paxscan_80kvp.txt', delimiter = ', ')
h = np.swapaxes(h,0,1)        

z = np.loadtxt('/users/Baier/analysis_files/linespread_acquisition_flat_measurement_paxscan_60kvp.txt', delimiter = ', ')
z = np.swapaxes(z,0,1)    
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath = '/users/Baier/analysis_files/'
filepath_pic = '/users/Baier/analysis_pictures/combi_spotsizes/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
# write the header for txt file
info = ''
#============================================================================================================
#                       analyis params
#============================================================================================================
fit_esf_v = True
fit_esf_h = True
fit_esf_d = True
comparison = True
plt.close('all')
#============================================================================================================
#                        Data handling
#============================================================================================================

"""
40 kVp data
"""
watts40 = f[0]
watts40 = np.insert(watts40 ,0,0)
vert_spotsize40 = f[1]
vert_spotsize40 = np.insert(vert_spotsize40,0,0)
hor_spotsize40 = f[2]
hor_spotsize40 = np.insert(hor_spotsize40,0,0)
diag_spotsize40 = f[3]
diag_spotsize40 = np.insert(diag_spotsize40,0,0)

"""
60 kVp data control
"""
watts60 = g[0]
watts60 = np.insert(watts60 ,0,0)
vert_spotsize60 = g[1]
vert_spotsize60 = np.insert(vert_spotsize60,0,0)
hor_spotsize60 = g[2]
hor_spotsize60 = np.insert(hor_spotsize60,0,0)
diag_spotsize60 = g[3]
diag_spotsize60 = np.insert(diag_spotsize60,0,0)

"""
80 kVp data
"""
watts80 = h[0]
watts80 = np.insert(watts80 ,0,0)
vert_spotsize80 = h[1]
vert_spotsize80 = np.insert(vert_spotsize80,0,0)
hor_spotsize80 = h[2]
hor_spotsize80 = np.insert(hor_spotsize80,0,0)
diag_spotsize80 = h[3]
diag_spotsize80 = np.insert(diag_spotsize80,0,0)

"""
60 kVp data old
"""
watts = z[0]
watts = np.insert(watts ,0,0)
vert_spotsize = z[1]
vert_spotsize = np.insert(vert_spotsize,0,0)
hor_spotsize = z[2]
hor_spotsize = np.insert(hor_spotsize,0,0)
diag_spotsize = z[3]
diag_spotsize = np.insert(diag_spotsize,0,0)

if fit_esf_h:
    ## horizontal spots of all three energies
    x_max_h = (1.1*np.max(watts80))
    y_max_h = (1.05*np.max(hor_spotsize60))
    h = [ 0., x_max_h,  0., y_max_h]
    plt.figure('Results HORIZONTAL spot-size 40,60,80kVp')
    plt.axis(h)
    plt.plot(watts40, hor_spotsize40, lw=1 )
    plt.plot(watts60, hor_spotsize60, lw=1,)
    plt.plot(watts80, hor_spotsize80, lw=1,)
    plt.title('Horizontal spot size 40,60,80kVp')
    plt.xticks( np.arange(0, x_max_h,20),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_h,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '40 kVp','60 kVp','80 kVp')) ,loc = 0,title = 'Horizontal spot-sizes')
    plt.savefig(filepath_pic+'Results Horizontal spot-sizes  40,60,80kvp.pdf', format = 'pdf',dpi = 300)
if fit_esf_v:
    ## vertical spots of all three energies
    x_max_v = (1.1*np.max(watts80))
    y_max_v = (1.05*np.max(vert_spotsize60))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results VERTICAL spot-size 40,60,80kVp')
    plt.axis(v)
    plt.plot(watts40, vert_spotsize40, lw=1)
    plt.plot(watts60, vert_spotsize60, lw=1)
    plt.plot(watts80, vert_spotsize80, lw=1)
    plt.title('Vertical spot size 40,60,80kVp')
    plt.xticks( np.arange(0, x_max_v,20),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,10),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '40 kVp','60 kVp','80 kVp')) ,loc = 0,title = 'Vertical spot-sizes')
    plt.savefig(filepath_pic+'Results Vertical spot-sizes  40,60,80kvp.pdf', format = 'pdf',dpi = 300)
if fit_esf_d:
    ## diagonal spots of all three energies
    x_max_d = (1.1*np.max(watts80))
    y_max_d = (1.05*np.max(diag_spotsize60))
    d = [ 0., x_max_d,  0., y_max_d]
    plt.figure('Results DIAGONAL spot-size 40,60,80kVp')
    plt.axis(d)
    plt.plot(watts40, diag_spotsize40, lw=1)
    plt.plot(watts60, diag_spotsize60, lw=1)
    plt.plot(watts80, diag_spotsize80, lw=1)
    plt.title('Diagonal spot size 40,60,80kVp')
    plt.xticks( np.arange(0, x_max_d,20) ,fontsize  = 12)
    plt.yticks( np.arange(0, y_max_d,20) ,fontsize  = 12)
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '40 kVp','60 kVp','80 kVp')) ,loc = 0,title = 'Diagonal spot-sizes')
    plt.savefig(filepath_pic+'Results Diagonal spot-sizes  40,60,80kvp.pdf', format = 'pdf',dpi = 300)


if comparison:
    ## horizontal spots of 60 kvp
    x_max_h = (1.1*np.max(watts60))
    y_max_h = (1.05*np.max(hor_spotsize60))
    h = [ 0., x_max_h,  0., y_max_h]
    plt.figure('Results HORIZONTAL spot-size 60kVp')
    plt.axis(h)
    plt.plot(watts60, hor_spotsize60, lw=1 )
    plt.plot(watts, hor_spotsize, lw=1,)
    plt.title('Horizontal spot size 60kVp')
    plt.xticks( np.arange(0, x_max_h,10),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_h,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '60 kVp control','60 kVp')) ,loc = 0,title = 'Horizontal spot-sizes')
    plt.savefig(filepath_pic+'comparison Horizontal spot-sizes  60kvp.pdf', format = 'pdf',dpi = 300)
    
    ## vertical spots of 60 kvp
    x_max_v = (1.1*np.max(watts60))
    y_max_v = (1.05*np.max(vert_spotsize60))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results VERTICAL spot-size 60kvp')
    plt.axis(v)
    plt.plot(watts60, vert_spotsize60, lw=1)
    plt.plot(watts, vert_spotsize, lw=1)
    plt.title('Vertical spot size 60kvp')
    plt.xticks( np.arange(0, x_max_v,10),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '60 kVp control','60 kVp')) ,loc = 0,title = 'Vertical spot-sizes')
    plt.savefig(filepath_pic+'comparison Vertical spot-sizes  60kvp.pdf', format = 'pdf',dpi = 300)
    
    ## diagonal spots of 60 kvp
    x_max_d = (1.1*np.max(watts60))
    y_max_d = (1.05*np.max(diag_spotsize60))
    d = [ 0., x_max_d,  0., y_max_d]
    plt.figure('Results DIAGONAL spot-size 60kvp')
    plt.axis(d)
    plt.plot(watts60, diag_spotsize60, lw=1)
    plt.plot(watts, diag_spotsize, lw=1)
    plt.title('Diagonal spot size 60kvp')
    plt.xticks( np.arange(0, x_max_d,10) ,fontsize  = 12)
    plt.yticks( np.arange(0, y_max_d,20) ,fontsize  = 12)
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( '60 kVp control','60 kVp')) ,loc = 0,title = 'Diagonal spot-sizes')
    plt.savefig(filepath_pic+'comparison Diagonal spot-sizes  60kvp.pdf', format = 'pdf',dpi = 300)

    ## spots of 60 kvp_control
    x_max_v = (1.1*np.max(watts60))
    y_max_v = (1.05*np.max(hor_spotsize60))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results 60kvp_control')
    plt.axis(v)
    plt.plot(watts60, vert_spotsize60, lw=1)
    plt.plot(watts60, hor_spotsize60, lw=1)
    plt.plot(watts60, diag_spotsize60, lw=1)
    plt.title('Spot sizes 60kVp_control')
    plt.xticks( np.arange(0, x_max_v,10),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( 'Vertical','horizontal','diagonal')) ,loc = 0,title = 'Spot-sizes @ 60kVp')
    plt.savefig(filepath_pic+'comparison spot-sizes  60kvp_control.pdf', format = 'pdf',dpi = 300)
    
    ## spots of 60 kvp
    x_max_v = (1.1*np.max(watts))
    y_max_v = (1.05*np.max(hor_spotsize))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results 60kvp')
    plt.axis(v)
    plt.plot(watts, vert_spotsize, lw=1)
    plt.plot(watts, hor_spotsize, lw=1)
    plt.plot(watts, diag_spotsize, lw=1)
    plt.title('Spot sizes 60kVp')
    plt.xticks( np.arange(0, x_max_v,10),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( 'Vertical','horizontal','diagonal')) ,loc = 0,title = 'Spot-sizes @ 60kVp')
    plt.savefig(filepath_pic+'comparison spot-sizes  60kvp.pdf', format = 'pdf',dpi = 300)
    
    ## spots of 60 kvp control/normal
    x_max_v = (1.1*np.max(watts))
    y_max_v = (1.05*np.max(hor_spotsize60))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results 60kvp')
    plt.axis(v)
    plt.plot(watts, vert_spotsize, lw=1,label = 'vertical', color = 'b')
    plt.plot(watts, hor_spotsize, lw=1, label = 'horizontal',color = 'm')
    plt.plot(watts, diag_spotsize, lw=1, label = 'diagonal', color = 'g')
    plt.plot(watts60, vert_spotsize60, lw=1, label = 'vertical_control', color = 'c')
    plt.plot(watts60, hor_spotsize60, lw=1, label = 'horizontal_control',color = 'r')
    plt.plot(watts60, diag_spotsize60, lw=1, label = 'diagonal_control',color = 'y')
    plt.title('Spot sizes 60kVp')
    plt.xticks( np.arange(0, x_max_v,10),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend(fontsize = 10,loc = 0,title = 'Spot-sizes @ 60kVp')
    plt.savefig(filepath_pic+'comparison spot-sizes  60kvpcombi.pdf', format = 'pdf',dpi = 300)
    
    ## spots of 40 kvp
    x_max_v = (1.1*np.max(watts40))
    y_max_v = (1.05*np.max(hor_spotsize40))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results 40kvp')
    plt.axis(v)
    plt.plot(watts40, vert_spotsize40, lw=1)
    plt.plot(watts40, hor_spotsize40, lw=1)
    plt.plot(watts40, diag_spotsize40, lw=1)
    plt.title('Spot sizes 40kVp')
    plt.xticks( np.arange(0, x_max_v,5),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,10),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( 'Vertical','horizontal','diagonal')) ,loc = 0,title = 'Spot-sizes @ 40kVp')
    plt.savefig(filepath_pic+'comparison spot-sizes  40kvp.pdf', format = 'pdf',dpi = 300)
    
    ## spots of 80 kvp
    x_max_v = (1.1*np.max(watts80))
    y_max_v = (1.05*np.max(hor_spotsize80))
    v = [ 0., x_max_v,  0., y_max_v]
    plt.figure('Results 80kvp')
    plt.axis(v)
    plt.plot(watts80, vert_spotsize80, lw=1)
    plt.plot(watts80, hor_spotsize80, lw=1)
    plt.plot(watts80, diag_spotsize80, lw=1)
    plt.title('Spot sizes 80kVp')
    plt.xticks( np.arange(0, x_max_v,20),fontsize  = 12 )
    plt.yticks( np.arange(0, y_max_v,20),fontsize  = 12 )
    plt.ylabel('Spot size [$\mu$m]',fontsize = 12)
    plt.xlabel('Power [W]',fontsize = 12)
    plt.legend((( 'Vertical','horizontal','diagonal')) ,loc = 0,title = 'Spot-sizes @ 80kVp')
    plt.savefig(filepath_pic+'comparison spot-sizes  80kvp.pdf', format = 'pdf',dpi = 300)