# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 17:29:38 2015

@author: ga56pan
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
import sys
import os
#sys.path.append('/users/schaff/Python Scripts/')
#import pyFlorian as pyF

'''
Script for picture generation for time power stability 
'''
#============================================================================================================
#                        Analysis options
#============================================================================================================
power = False
exp_time = False
measure281014 = False
measure291014 = False
oneweek = False
snakehead = False
limestone = True
plt.close("all")

#============================================================================================================
#                        Data input
#============================================================================================================

""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
""" Folder containing the data"""
data_files = "/users/Baier/analysis_pictures/time_power_stability/"
#============================================================================================================
#                        Data input
#============================================================================================================
if power:
    power_stuff = np.load(data_files+"stability_datapower.npy")
if exp_time:
    exp_stuff = np.load(data_files+"stability_datamean_vs_exptime.npy")
if measure281014:
    data281014 = np.load(data_files+"stability_datatime_stability_measurement_281014.npy")
if measure291014:
    data291014 = np.load(data_files+"stability_datatime_stability_measurement_291014.npy")
if oneweek:
    data_oneweek = np.load(data_files+"stability_datayashxttz_one_week.npy")
if snakehead:
    data_snake = np.load(data_files+"stability_datasnakehead_flo.npy")
if limestone:
    data_lime =np.load(data_files+"stability_datalimestone_fritz.npy")
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath_p = '/users/Baier/analysis_pictures/time_power_stability/'
filepath_d = '/users/Baier/analysis_files/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
# write the header for txt file
info = 'contains at first column the different times in hours, mean_visibility of the one week measurement' 
#plt.rc('text', usetex  = True)
#plt.rc('font', family='serif')   
if power:
    ### do linear regression for each dataset    
    regression1 = np.polyfit(power_stuff[0],power_stuff[1],1)
    fit_og_mp = np.poly1d(regression1)
    regression2 = np.polyfit(power_stuff[0],power_stuff[2],1)
    fit_mg_mp = np.poly1d(regression2)
    regression3 = np.polyfit(power_stuff[0],power_stuff[3]*2.,1)
    fit_og_mp_new = np.poly1d(regression3)
    regression4 = np.polyfit(power_stuff[0],power_stuff[4]*2.,1)
    fit_mg_mp_new = np.poly1d(regression4)
    
    #rc('text', usetex=True)
    plt.figure('Meancounts vs Power')               ###power_new  = power_stuff[0]
    plt.plot(power_stuff[0],power_stuff[1], 'go')   ###meancounts_og_mp = power_stuff[1]
    plt.plot(power_stuff[0],power_stuff[2], 'bo')   ###meancounts_mg_mp = power_stuff[2]
    plt.plot(power_stuff[0],power_stuff[3]*2.,'yo')    ###meancounts_og_mp_new = power_stuff[3]
    plt.plot(power_stuff[0],power_stuff[4]*2.,'ko')    ###meancounts_mg_mp_new = power_stuff[4]
    plt.plot(power_stuff[0],fit_og_mp(power_stuff[0]),'r',power_stuff[0],fit_mg_mp(power_stuff[0]),'r',
             power_stuff[0],fit_og_mp_new(power_stuff[0]),'r',power_stuff[0],fit_mg_mp_new(power_stuff[0]),'r', lw =.3 )
    #plt.title("Meancounts vs Power")
    plt.xticks( np.arange(0, np.max(power_stuff[0])+10,10) )
    plt.xticks(fontsize  = 8)
    plt.xlim(0,np.max(power_stuff[0])*1.01)
    plt.yticks( np.arange(0, np.max(power_stuff[1]),10000) )
    plt.yticks(fontsize  = 8)
    plt.ylim(0,np.max(power_stuff[1])*1.01)
    plt.ylabel(r'Mean intensity' r' $\bar{I}$ [arb.unit]',fontsize = 12)
    plt.xlabel(r'Power [W]',fontsize = 12)
    plt.legend(( 'results without gratings 22.10.2014','Results with gratings 22.10.2014',
    'Results without gratings 02.07.2015','results with gratings 02.07.215','linear regression of the data') ,loc = 0,fontsize =12)    
    plt.savefig(filepath_p+'Meancounts vs Power with regression.pdf', format = 'pdf',dpi = 300 )
        
if exp_time:
    ### regression for data 
    
    regression = np.polyfit( exp_stuff[0],exp_stuff[1],1)
    fit_expt = np.poly1d(regression)
    plt.figure('Meancounts vs exposure time')            ###time_s = exp_stuff[0]
    plt.plot(exp_stuff[0],exp_stuff[1],'bo')             ###meancounts_met = exp_stuff[1] 
    plt.plot(exp_stuff[0],fit_expt(exp_stuff[0]),'r',lw = .3)
#    plt.loglog(exp_stuff[0],exp_stuff[1],'bo')
#    plt.loglog(exp_stuff[0],fit_expt(exp_stuff[0]),'r',lw = .3)
    #plt.title("Meancounts vs exposure time")
    plt.xticks( np.arange(0, np.max(exp_stuff[0])+5,5) )
    plt.xticks(fontsize  = 12)
    plt.xlim(0,np.max(exp_stuff[0])*1.01)
    plt.yticks( np.arange(0, np.max(exp_stuff[1]),5000) )
    plt.yticks(fontsize  = 12)
    plt.ylim((0,np.max(exp_stuff[1])*1.01))
    plt.ylabel(r'Mean intensity $\bar{I}$ [arb.unit]',fontsize = 12)
    plt.xlabel(r'Time t [s]',fontsize = 12)
    plt.legend(('measured data', 'linear regression of the data') ,loc = 0,fontsize = 12)    
    plt.savefig(filepath_p+'Meancounts vs Exposure time.pdf', format = 'pdf',dpi = 300 )
    
if measure281014:
    plt.figure('Time stability measurement 28.10.2014 normed during 9 hours')   ###counttime = data281014[0]
    #plt.plot(data281014[0]/3600,data281014[1],lw = .5)                       ###meancounts = data281014[1]
    plt.plot(data281014[0]/3600,data281014[2],lw = .5)                      ###meancountsnorm first value = data281014[2]
    #plt.title("Meancounts vs Time during 9 hours")
    plt.xticks( np.arange(0, (np.max(data281014[0])/3600)+1,1) )
    plt.xticks(fontsize  = 8)
    plt.xlim(0,(np.max(data281014[0])/3600)*1.01)
    
    #plt.yticks( np.arange(900, np.max(data281014[1]),10) )
    plt.yticks( np.arange(1, np.max(data281014[2]),.01) )
    #plt.ylim(900,np.max(data281014[1])*1.01)
    plt.ylim(1,np.max(data281014[2])*1.01)            
    plt.yticks(fontsize  = 8)
    plt.ylabel(r'Mean intensity' r' $\bar{I}$ [arb.unit]',fontsize = 12)
    plt.xlabel('Measurement time [h]',fontsize = 12)
    plt.legend(('Meancounts per flatfield normed to initial value',) ,loc = 0,fontsize = 12)
    #plt.savefig(filepath_p+'Meancounts vs Time 9 hours period.pdf', format = 'pdf',dpi = 300 )
    plt.savefig(filepath_p+'Meancounts vs Time 9 hours period norm.pdf', format = 'pdf',dpi = 300 )

if measure291014:
    plt.figure('Time stability measurement 29.10.2014 during 14 hours')   ###counttime = data291014[0]
    plt.plot(data291014[0]/3600,data291014[1],lw = .5)                       ###meancounts = data291014[1]
    #plt.plot(data291014[0]/3600,data291014[2],lw = .5)                      ###meancountsnorm first value = data291014[2]
    #plt.title("Meancounts vs Time during 14 hours")
    plt.xticks( np.arange(0, (np.max(data291014[0])/3600)+1,1) )
    plt.xticks(fontsize  = 8)
    plt.xlim(0,(np.max(data291014[0])/3600)*1.01)
    
    plt.yticks( np.arange(890, np.max(data291014[1]),10) )
    #plt.yticks( np.arange(.93, np.max(data291014[2]),.01) )
    plt.ylim(890,np.max(data291014[1])*1.01)
    #plt.ylim(.93,np.max(data291014[2])*1.01)            
    plt.yticks(fontsize  = 8)
    plt.ylabel(r'Mean intensity' r' $\bar{I}$ [arb.unit]',fontsize = 12)
    plt.xlabel('Measurement time [h]',fontsize = 12)
    plt.legend(('Meancounts per flatfield',) ,loc = 0,fontsize = 12)
    plt.savefig(filepath_p+'Meancounts vs Time 14 hours period.pdf', format = 'pdf',dpi = 300 )
    #e17.utils.franzmap()
    #plt.savefig(filepath_p+'Meancounts vs Time 14 hours period norm.pdf', format = 'pdf',dpi = 300 )

if oneweek:
    plt.figure('Measurement over a hole week')                  ###deltas(time) = data_oneweek[0] 
    plt.plot(data_oneweek[0], data_oneweek[1], lw = .5)
    #plt.plot(data_oneweek[0], data_oneweek[1]/data_oneweek[1,0], lw = .5)          ###mean_visibility = data_oneweek[1]
    #plt.title('Mean visibility during one week ')                        ###mean_visibility_norm = data_oneweek[2] norm to max
    plt.xticks( np.arange(0, (np.max(data_oneweek[0])+1),10) )       ###mean a0 = data_oneweek[3]
    plt.xticks(fontsize  = 8)                                   ###mean a1 = data_oneweek[4]
    plt.xlim(0,np.max(data_oneweek[0])*1.01)                           ###meancounts_average_flatfields = data_oneweek[5] unprocessed
                     
    plt.yticks( np.arange(0.23, np.max(data_oneweek[1]),.01) )
    #plt.yticks( np.arange(0.85, np.max(data_oneweek[1])/data_oneweek[1,0],.01) )     ###meancounts_average_flatfields_norm = data_oneweek[6] unprocessed
    plt.ylim(0.23,np.max(data_oneweek[1])*1.01)            
    #plt.ylim(0.85,(np.max(data_oneweek[1])/data_oneweek[1,0])*1.01)
    plt.yticks(fontsize  = 8)                         
    plt.ylabel('Relative visibility [\%]',fontsize = 12)                      
    plt.xlabel('Measurement time [h]',fontsize = 12)                                    
    plt.legend(('Relative visibility',) ,loc = 0,fontsize = 12)
    plt.savefig(filepath_p+'Meanvisibility over a whole week.pdf', format = 'pdf',dpi = 300 )

if snakehead:
    plt.figure('Time stability measurement of snake head over 4 hours')   ###counttime = data_snake[0]
    #plt.plot(data_snake[0],data_snake[1],lw = .5)                       ###meancounts = data_snake[1]
    plt.plot(data_snake[0],data_snake[2],lw = .5)                      ###meancountsnorm first value = data_snake[2]
    #plt.title("Meancounts vs Time during 4 hours")
    plt.xticks( np.arange(0, np.max(data_snake[0])+1,.5) )
    plt.xticks(fontsize  = 8)
    plt.xlim(0,(np.max(data_snake[0]))*1.01)
    
    #plt.yticks( np.arange(11550, np.max(data_snake[1]),15) )
    plt.yticks( np.arange(.999, np.max(data_snake[2]),.001) )
    #plt.ylim(11550,np.max(data_snake[1])*1.001)
    plt.ylim(.999,np.max(data_snake[2])*1.001)            
    plt.yticks(fontsize  = 8)
    plt.ylabel('Meancounts [arb. unit]',fontsize = 12)
    plt.xlabel('Measurement time [h]',fontsize = 12)
    plt.legend(('Meancounts per flatfield normed to initial value',) ,loc = 0,fontsize = 12)
    #plt.savefig(filepath_p+'Meancounts vs Time 4 hours period.pdf', format = 'pdf',dpi = 300 )
    #e17.utils.franzmap()
    plt.savefig(filepath_p+'Meancounts vs Time 4 hours period norm .pdf', format = 'pdf',dpi = 300 )

if limestone:
    plt.figure('Measurement over two days')                          ###deltas(time) = data_lime[0] 
    plt.plot(data_lime[0], data_lime[1], lw = .5)
    #plt.plot(data_lime[0], data_lime[2], lw =.5)                   ###mean_visibility = data_lime[1]
    #plt.title('Mean visibility during two ')                        ###mean_visibility_norm = data_lime[2] norm to start
    plt.xticks( np.arange(0, (np.max(data_lime[0])+1),5) )          ### meancounts_average_flatfields = data_lime[3]
    plt.xticks(fontsize  = 12)                                       ###meancounts_average_flatfields_norm = data_lime[4]
    plt.xlim(0,np.max(data_lime[0])*1.01)
                     
    plt.yticks( np.arange(0.21, np.max(data_lime[1])*1.01,.001) )
    #plt.yticks( np.arange(0.99, np.max(data_lime[2])*1.005,.005) )     
    plt.ylim(0.21,np.max(data_lime[1])*1.01)            
    #plt.ylim(0.99,(np.max(data_lime[2]))*1.005)
    plt.yticks(fontsize  = 12)                         
    plt.ylabel('Normed visibility [$\%$]',fontsize = 12)                      
    plt.xlabel('Measurement time [h]',fontsize = 12)                                    
    plt.legend(('Normed visibility',) ,loc = 0,fontsize = 12)
    plt.savefig(filepath_p+'Relative visibility over two days.pdf', format = 'pdf',dpi = 300 )
    #plt.savefig(filepath_p+'Relative visibility over two days normed.pdf', format = 'pdf',dpi = 300 )
    
    plt.figure('Measurement over two days meancounts')                  ###deltas(time) = data_lime[0] 
    plt.plot(data_lime[0], data_lime[3], lw = .5)
    #plt.plot(data_lime[0], data_lime[4], lw =.5)                       ###mean_visibility = data_lime[1]
    #plt.title('Mean visibility during two ')                        ###mean_visibility_norm = data_lime[2] norm to start
    plt.xticks( np.arange(0, (np.max(data_lime[0])+1),5) )           ### meancounts_average_flatfields = data_lime[3]
    plt.xticks(fontsize  = 12)                                       ###meancounts_average_flatfields_norm = data_lime[4]
    plt.xlim(0,np.max(data_lime[0])*1.01) 
                    
    plt.yticks( np.arange(35600, np.max(data_lime[3])*1.001,50) )
    #plt.yticks( np.arange(0.995, np.max(data_lime[4]),.001) )     
    plt.ylim(35600,np.max(data_lime[3])*1.001)            
    #plt.ylim(0.995,(np.max(data_lime[4]))*1.001)
    plt.yticks(fontsize  = 12)                         
    plt.ylabel(r'Mean intensity normed ' r' $\bar{I}$ [arb.unit]',fontsize = 12)                      
    plt.xlabel('Measurement time [h]',fontsize = 12)                                    
    plt.legend(('Normed intensity',) ,loc = 0,fontsize = 12)
    plt.savefig(filepath_p+'Meancounts over two days.pdf', format = 'pdf',dpi = 300 )
    #plt.savefig(filepath_p+'Meancounts over two days normed.pdf', format = 'pdf',dpi = 300 )
    
    