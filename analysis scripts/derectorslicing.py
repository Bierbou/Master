# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 13:54:25 2014

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/users/schaff/Python Scripts/')
import pyFlorian as pyF

"""

Script for generatig sliced detectorpicures
checking constantness of pixelcounts

"""
#Meancounts vs picturenumber
#detector sliced in 4 parts
plt.close("all")

"""
four pieces detectorslicing

"""
data = np.zeros((938,500,500))

for i in range(938):
    print "Loading file", (1865795+i)
    data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_15102014/paxscan/ct/paximage_ct_%i.h5"%(1865795+i))["raw_data"]
   
picnumb = np.asarray(range(1,939,))
meancounts1 = np.mean(data[:,0:250,0:250],axis=(1,2))
meancounts2 = np.mean(data[:,0:250,250:500],axis=(1,2))
meancounts3 = np.mean(data[:,250:500,0:250],axis=(1,2))
meancounts4 = np.mean(data[:,250:500,250:500],axis=(1,2))
#    
#left top
plt.figure(1)
#plt.title("Meancounts vs fps") 
pyF.fit_linear(meancounts1, x_data = picnumb,plot=True)  
##pyF.fit_linear(meancounts[:-4], x_data = frames[:-4],plot=True)

#right top
plt.figure(2)
#plt.title("Meancounts vs fps") 
pyF.fit_linear(meancounts2, x_data = picnumb,plot=True)  
##pyF.fit_linear(meancounts[:-4], x_data = frames[:-4],plot=True)

#left bottom
plt.figure(3)
#plt.title("Meancounts vs fps") 
pyF.fit_linear(meancounts3, x_data = picnumb,plot=True)  
##pyF.fit_linear(meancounts[:-4], x_data = frames[:-4],plot=True)

#right bottom
plt.figure(4)
#plt.title("Meancounts vs fps") 
pyF.fit_linear(meancounts4, x_data = picnumb,plot=True)  
##pyF.fit_linear(meancounts[:-4], x_data = frames[:-4],plot=True)

e17.utils.plot_3d_array(data)

"""
16 pieces detectorslicing

"""

#data = np.zeros((938,500,500))
#
#for i in range(938):
#    print "Loading file", (1865795+i)
#    data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_15102014/paxscan/ct/paximage_ct_%i.h5"%(1865795+i))["raw_data"]
#   
#picnumb = np.asarray(range(1,939,))
#meancounts1  = np.mean(data[:,0:125,0:125],axis=(1,2))
#meancounts2  = np.mean(data[:,0:125,125:250],axis=(1,2))
#meancounts3  = np.mean(data[:,0:125,250:375],axis=(1,2))
#meancounts4  = np.mean(data[:,0:125,375:500],axis=(1,2))
#meancounts5  = np.mean(data[:,125:250,0:125],axis=(1,2))
#meancounts6  = np.mean(data[:,125:250,125:250],axis=(1,2))
#meancounts7  = np.mean(data[:,125:250,250:375],axis=(1,2))
#meancounts8  = np.mean(data[:,125:250,375:500],axis=(1,2))
#meancounts9  = np.mean(data[:,250:375,0:125],axis=(1,2))
#meancounts10 = np.mean(data[:,250:375,125:250],axis=(1,2))
#meancounts11 = np.mean(data[:,250:375,250:375],axis=(1,2))
#meancounts12 = np.mean(data[:,250:375,375:500],axis=(1,2))
#meancounts13 = np.mean(data[:,375:500,0:125],axis=(1,2))
#meancounts14 = np.mean(data[:,375:500,125:250],axis=(1,2))
#meancounts15 = np.mean(data[:,375:500,250:375],axis=(1,2))
#meancounts16 = np.mean(data[:,375:500,375:500],axis=(1,2))
##    
##left top
#plt.figure(1)
#pyF.fit_linear(meancounts1, x_data = picnumb,plot=True)  
#
##left middle top
#plt.figure(2)
#pyF.fit_linear(meancounts2, x_data = picnumb,plot=True)  
#
##right middle
#plt.figure(3)
#pyF.fit_linear(meancounts3, x_data = picnumb,plot=True)  
#
##right top
#plt.figure(4) 
#pyF.fit_linear(meancounts4, x_data = picnumb,plot=True)  
#
##left mid top
#plt.figure(5)
#pyF.fit_linear(meancounts5, x_data = picnumb,plot=True)  
#
##left mid middle top
#plt.figure(6)
#pyF.fit_linear(meancounts6, x_data = picnumb,plot=True)  
#
##right mid middle top
#plt.figure(7)
#pyF.fit_linear(meancounts7, x_data = picnumb,plot=True)  
#
##right mid top
#plt.figure(8)
#pyF.fit_linear(meancounts8, x_data = picnumb,plot=True)  
#
##left mid bottom
#plt.figure(9)
#pyF.fit_linear(meancounts9, x_data = picnumb,plot=True)  
#
## left mid middle bottom
#plt.figure(10)
#pyF.fit_linear(meancounts10, x_data = picnumb,plot=True)  
#
##right mid middle bottom
#plt.figure(11)
#pyF.fit_linear(meancounts11, x_data = picnumb,plot=True)  
#
##right mid bottom
#plt.figure(12)
#pyF.fit_linear(meancounts12, x_data = picnumb,plot=True)  
#
##left bottom
#plt.figure(13) 
#pyF.fit_linear(meancounts13, x_data = picnumb,plot=True)  
#
##left middle bottom
#plt.figure(14) 
#pyF.fit_linear(meancounts14, x_data = picnumb,plot=True)  
#
##right middle bottom
#plt.figure(15)
#pyF.fit_linear(meancounts15, x_data = picnumb,plot=True)  
#
##right bottom
#plt.figure(16)
#pyF.fit_linear(meancounts16, x_data = picnumb,plot=True)  
#
## plot of all picutres in a row
#e17.utils.plot_3d_array(data)






"""
Code for  focalsizedetermination
error_fit_func = lambda p, x: (p[0]*scs.erf((x-p[1])/(math.sqrt(2)*p[2]))+p[3])
errorfunction = lambda p, x, y: error_fit_func(p, x) - y

p0 = [1.,0.,1.,0.]
p,cov_p,infodict,mesg,ier = opti.leastsq(errorfunction, p0[:] , args = (pixel_number, linesum[0]), full_output = True, maxfev = 10000)
print p 
print cov_p
"""


