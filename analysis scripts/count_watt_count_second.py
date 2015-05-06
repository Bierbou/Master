# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/ga56pan/.spyder2/.temp.py
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.path.append('/users/schaff/Python Scripts/')
#import pyFlorian as pyF


"""
Script for generating plots with linear fits for 
meancounts vs power, meancounts vs time data 
and meancounts vs framerate
"""
plt.close("all")

short_m = False
long_m = False
oneweek = True
plots = True
## Meancounts vs Power
## output in figure 1
##data = np.zeros((51,800,800))
#data = np.zeros((14,800,800))
##for slicing
##data = np.zeros((14,500,500))
#for i in range(14):
#    print "Loading file", (1867997+i)
#    data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/power_test_22102014/paxscan/ct/paximage_ct_%i.h5"%(1867997+i))["raw_data"]
#
##power = np.array([1,5,10,15,20,25,30,35,40,45,50,75,100,150])
##power = np.asarray(range(0,102,2))
#power = np.asarray(range(0,28,2))
##meancounts = np.mean(data[:,250:600,200:600],axis=(1,2))
#meancounts = np.mean(data[:,:,:],axis=(1,2))
##plt.figure(1)    
##plt.plot(power,meancounts)
##plt.title("Meancounts vs Power")    
#pyF.fit_linear(meancounts, x_data = power,plot=True)
##plt.figure(4)  
##plt.imshow(np.mean(data, axis =0), vmax = 10000)
#e17.utils.plot_3d_array(data, vmax = 80000)


## Meancounts vs ExposureTime
## output in figure 2
#data = np.zeros((18,800,800))
#
#for i in range(18):
#    print "Loading file", (1865707+i)
#    data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/only_flat_fields/paxscan/ct/paximage_ct_%i.h5"%(1865707+i))["raw_data"]
#
#time = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,30,45,60])
#meancounts = np.mean(data[:,250:600,200:600],axis=(1,2))
#    
#
#plt.figure(2)
##plt.title("Meancounts vs Power") 
#pyF.fit_linear(meancounts, x_data = time,plot=True)  
##pyF.fit_linear(meancounts[:-4], x_data = time[:-4],plot=True)



#Meancounts vs picturenumber
#output in figure 3
"""
short measurement on 28102014
"""
if short_m:
    ##data = np.zeros((2999,800,800)) # including pictures without ilumination 1886749
    data = np.zeros((2958,800,800)) # excluding pictures without ilumination 1886790
    for i in range(2958):
        print "Loading file", (1886790+i)
        data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_28102014/paxscan/ct/paximage_ct_%i.h5"%(1886790+i))["raw_data"]
        
    picnumb = np.asarray(range(41,2999))
    meancount = np.mean(data[:,:,:],axis=(1,2))
    meancounts = meancount/meancount[0]
    plt.figure(3)
    plt.plot(picnumb,meancounts)
    e17.utils.plot_3d_array(data, vmax = 80000)
    #pyF.fit_linear(meancounts, x_data = picnumb,plot=True) 
    """
long measurement on 29102014
"""
if long_m:
    data = np.zeros((4199,800,800)) # including pictures without ilumination 1890127
    data = np.zeros((4181,800,800)) # excluding pictures without ilumination 1890145
    for i in range(4181):
        print "Loading file", (1890145+i)
        data[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_29102014/paxscan/ct/paximage_ct_%i.h5"%(1890145+i))["raw_data"]
        
    picnumb = np.asarray(range(19,4200))
    meancount = np.mean(data[:,:,:],axis=(1,2))
    meancounts = meancount/meancount[0]
    plt.figure(4)
    plt.plot(picnumb, meancounts)
    #pyF.fit_linear(meancounts, x_data = picnumb,plot=True)  
    #pyF.fit_linear(meancounts[:-4], x_data = frames[:-4],plot=True)
    

"""
measurment over a hole week sponsored by yashs xtt2 flatfields
"""
if oneweek:
    sam_ID ='carbonfibremeshXTT2'
    path_to_data = '/data/DPC/local_setups/microfocus/samples/'+sam_ID+'/paxscan/'
    sub_dirs=os.listdir(path_to_data)
    
    flat_dirs=list()
    
    dir_num=1
    for sub_dir in sub_dirs:    
        dirs=np.sort(os.listdir(path_to_data+sub_dir))
        #One FlatField Dir, One Data Dir, One FlatField Dir, One Data Dir........
        for dir_name in dirs:
            if dir_num%4!=0:         
                flat_dirs.append(path_to_data+sub_dir+'/'+dir_name)
            dir_num=dir_num+1
    
    
    relevant_folders = (len(flat_dirs)-8)
    meancounts_average_flatfields = np.zeros((relevant_folders))
    deltas=np.empty((relevant_folders))
    deltas[0]=0
    t = 0
    for i in range(relevant_folders):
        meancounts = np.zeros((8))
        print 'Foldernumber', i
        for j in range(8):
            ## for undertanding    ####################flatfieldfolder[i]||[j] flatfieldnumber
            imgs = e17.io.h5read (flat_dirs[i]+'/'+os.listdir(flat_dirs[i])[j]) ["raw_data"]
            meancounts[j]= np.mean(imgs[200:650,150:600], axis = (0,1))

        meancounts_average_flatfields[i] = (sum (meancounts)/8.)
        
        if i< relevant_folders-1:
            delta_t = os.stat(flat_dirs[i+1]).st_mtime - os.stat(flat_dirs[0]).st_mtime
            t=t+delta_t
            deltas[i+1]=delta_t/3600.
    meancounts_average_flatfields = (meancounts_average_flatfields/np.max(meancounts_average_flatfields))    
if plots:        
    plt.figure('Results')
    plt.plot(deltas, meancounts_average_flatfields, lw=2)
    plt.title('source time dependency')
    plt.ylabel('meancounts')
    plt.xlabel('time [hours]')

        
        
        