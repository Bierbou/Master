# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/home/ga56pan/.spyder2/.temp.py
"""

import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
import ddfSetupProcessing.HarmonicAnalysis.lsqprocessing as dpcprocess
import sys
import os
sys.path.append('/users/schaff/Python Scripts/')
import pyFlorian as pyF


"""
Script for generating plots with linear fits for 
meancounts vs power, meancounts vs time data 
and meancounts vs framerate
"""
plt.close("all")
#============================================================================================================
#                        Analysis options
#============================================================================================================
data_input = False
power = False
exp_time = False
short_m = False
long_m = False
oneweek = False
limestone = True
snakehead = False
unprocessed = False
processed = True
saveon = False
plots = False
#============================================================================================================
#                        Data input
#============================================================================================================
""" Standard data folder"""
data = "/data/DPC/local_setups/microfocus/samples/"
""" Folder containing the data"""
data_folder_ng = "meancounts_vs_power_new_ng"
data_folder_g = "meancounts_vs_power_new_g"


#data_files = "/users/Baier/analysis_pictures/time_power_stability/"
#power_stuff = np.load(data_files+"stability_datapower.npy")


if data_input:
    """   pictures without gratings """
    folder_ng = data+data_folder_ng+"/paxscan/" # data folder 
    sub_dirs_ng=os.listdir(folder_ng)
    data_ng =list() ### the real pictures
    dir_num_og=0
    for sub_dir_ng in sub_dirs_ng:    ##### check the subfolders and sort all pictures in it in the data folder 
        dirs_ng=np.sort(os.listdir(folder_ng+ sub_dir_ng)) 
        for dir_name in dirs_ng:                                                
            data_ng.append(folder_ng+ sub_dir_ng+'/'+dir_name)
            dir_num_og=dir_num_og+1
    """   pictures without gratings """
    folder_g = data+data_folder_g+"/paxscan/" # data folder 
    sub_dirs_g=os.listdir(folder_g)
    data_g =list() ### the real pictures
    dir_num_g=0
    for sub_dir_g in sub_dirs_g:    ##### check the subfolders and sort all pictures in it in the data folder 
        dirs_g=np.sort(os.listdir(folder_g+ sub_dir_g)) 
        for dir_name in dirs_g:                                                
            data_g.append(folder_g+ sub_dir_g+'/'+dir_name)
            dir_num_g=dir_num_g+1
#============================================================================================================
#                        Data output
#============================================================================================================
# choose your folder
filepath_p = '/users/Baier/analysis_pictures/time_power_stability/'
filepath_d = '/users/Baier/analysis_files/'
#filepath = data+'linespread_acquisition_flat_measurement3_60kvphighpower/'
# write the header for txt file
info = 'contains at first column the different times in hours, mean_visibility of the one week measurement'
power_range = 100
stepsize = 2   
powersteps = int(power_range/stepsize)



# Meancounts vs Power
if power:
    # output in figure 1
    data_og_mp = np.zeros((51,800,800)) #### og means without gratings
    data_mg_mp = np.zeros((51,500,500)) #### mg means with gratings
    data_og_mp_new = np.zeros((50,800,770)) #### og means without gratings 02.07.15
    data_mg_mp_new = np.zeros((50,800,770)) #### mg means with gratings
    #data = np.zeros((14,800,800))
    #for slicing
    #data = np.zeros((14,500,500))
    for i in range(51):
        print "Loading file", (1867997+i)#### 22.10.2014 without gratings ###15.10.2014 with gratings
        data_og_mp[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/power_test_22102014/paxscan/ct/paximage_ct_%i.h5"%(1867997+i))["raw_data"]
        data_mg_mp[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/power_test_15102014/paxscan/ct/paximage_ct_%i.h5"%(1865744+i))["raw_data"]
        power = np.asarray(range(0,102,2))
    meancounts_og_mp = np.mean(data_og_mp[:,250:600,200:600],axis=(1,2))
    meancounts_og_mp = np.delete(meancounts_og_mp,0)
    meancounts_og_mp = np.insert(meancounts_og_mp,0,0)
    meancounts_mg_mp = np.mean(data_mg_mp[:,250:600,200:600],axis=(1,2))
    meancounts_mg_mp = np.delete(meancounts_mg_mp,0)
    meancounts_mg_mp = np.insert(meancounts_mg_mp,0,0)
    for i in range(powersteps): 
        data_og_mp_new[i] = e17.io.h5read(data_ng[i])["raw_data"]
        data_mg_mp_new[i] = e17.io.h5read(data_g[i])["raw_data"]
    power_new = np.asarray(range(2,102,2))
    meancounts_og_mp_new = np.mean(data_og_mp_new[:,250:600,200:600],axis=(1,2))
    meancounts_og_mp_new = np.insert(meancounts_og_mp_new,0,0)
    meancounts_mg_mp_new = np.mean(data_mg_mp_new[:,250:600,200:600],axis=(1,2))
    meancounts_mg_mp_new = np.insert(meancounts_mg_mp_new,0,0)
    out_put_power = np.zeros((5,np.shape(power)[0]))
    out_put_power[0,:] = power
    out_put_power[1,:] = meancounts_og_mp
    out_put_power[2,:] = meancounts_mg_mp
    out_put_power[3,:] = meancounts_og_mp_new
    out_put_power[4,:] = meancounts_mg_mp_new
    ## Meancounts vs ExposureTime
if exp_time:
    data_met = np.zeros((18,800,800))
    for i in range(18):
        print "Loading file", (1865707+i)
        data_met[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/only_flat_fields_old_2/paxscan/ct/paximage_ct_%i.h5"%(1865707+i))["raw_data"]
    
    time_s = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,30,45,60])
    meancounts_met = np.mean(data_met[:,250:600,200:600],axis=(1,2))
    out_put_exp_time = np.zeros((2,np.shape(data_met)[0]))
    out_put_exp_time[0,:] = time_s
    out_put_exp_time[1,:] = meancounts_met
    
#Meancounts vs picturenumber
"""
short measurement on 28102014
"""
if short_m:
    ##data = np.zeros((2999,800,800)) # including pictures without ilumination 1886749
    data_mps = np.zeros((2958,800,800)) # excluding pictures without ilumination 1886790
    for i in range(2958):
        print "Loading file", (1886790+i)
        data_mps[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_28102014/paxscan/ct/paximage_ct_%i.h5"%(1886790+i))["raw_data"]
        
    time_span = 2958*11.
    seconds_mps = np.asarray(range(0,time_span,11))
    meancounts_s = np.mean(data_mps[:,:,:],axis=(1,2))
    meancounts_mps = meancounts_s/meancounts_s[0]
    
    out_put_281014 = np.zeros((3,np.shape(data_mps)[0]))
    out_put_281014[0,:] = seconds_mps
    out_put_281014[1,:] = meancounts_s
    out_put_281014[2,:] = meancounts_mps
 
    """
long measurement on 29102014
"""
if long_m: #### with gratings
    #data = np.zeros((4199,800,800)) # including pictures without ilumination 1890127
    data_mpl = np.zeros((4181,800,800)) # excluding pictures without ilumination 1890145
    for i in range(4181):
        print "Loading file", (1890145+i)
        data_mpl[i] = e17.io.h5read("/data/DPC/local_setups/microfocus/samples/stability_test_29102014/paxscan/ct/paximage_ct_%i.h5"%(1890145+i))["raw_data"]

    meancounts_l = np.mean(data_mpl[:,:,:],axis=(1,2))
    meancounts_mpl = meancounts_l/meancounts_l[0]
    time_span= 4181*12
    seconds_mpl = np.asarray(range(0,time_span,12))
    
    out_put_291014 = np.zeros((3,np.shape(data_mpl)[0]))
    out_put_291014[0,:] = seconds_mpl
    out_put_291014[1,:] = meancounts_l
    out_put_291014[2,:] = meancounts_mpl

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
    deltas=np.empty((relevant_folders))
    deltas[0]=0.
    t = 0
    if processed:
        a0 = np.ones((relevant_folders,450,450))
        a1 = np.ones((relevant_folders,450,450))
        phi = np.ones((relevant_folders,450,450))
        mean_visibility = np.ones(relevant_folders)
        #mean_visibility_norm = np.ones(relevant_folders)
        visibility = np.ones((relevant_folders,450,450))
        for i in range(relevant_folders):
            print 'Foldernumber', i
            ff = np.ones((8,450,450))
            for j in range(8):
                ## for undertanding    ####################flatfieldfolder[i]||[j] flatfieldnumber
                ff[j] = e17.io.h5read (flat_dirs[i]+'/'+os.listdir(flat_dirs[i])[j]) ["raw_data"][200:650,150:600]
            coeff = dpcprocess.lsq_fit(ff, nb_periods = 1, order = 2)
            a0[i] = coeff[0,0]
            a1[i] = coeff[0,1]
            phi[i] = coeff[1,1]
            visibility[i]=a1[i]/a0[i]
            visibility[i][np.isnan(visibility[i])]=1.
            visibility[i][np.isinf(visibility[i])]=1.
            mean_visibility[i] = np.mean(visibility[i])
        mean_visibility_norm = mean_visibility/np.max(mean_visibility)
    if unprocessed:
        meancounts_average_flatfields = np.zeros((relevant_folders))
        for i in range(relevant_folders):
            print 'Foldernumber', i
            meancounts = np.zeros((8))
            for j in range(8):
                imgs = e17.io.h5read (flat_dirs[i]+'/'+os.listdir(flat_dirs[i])[j]) ["raw_data"] 
                meancounts[j]= np.mean(imgs[200:650,150:600], axis = (0,1))
            meancounts_average_flatfields[i] = (sum (meancounts)/8.)
        meancounts_average_flatfields_norm = (meancounts_average_flatfields/np.max(meancounts_average_flatfields))
    for i in range(relevant_folders):
        if i< relevant_folders-1:
            delta_t = os.stat(flat_dirs[i+1]).st_mtime - os.stat(flat_dirs[0]).st_mtime
            t=t+delta_t
            deltas[i+1]=delta_t/3600.
    out_put_yashxtt = np.zeros((7,np.shape(mean_visibility)[0]))
    out_put_yashxtt[0,:] = deltas
    out_put_yashxtt[1,:] = mean_visibility
    out_put_yashxtt[2,:] = mean_visibility_norm
    out_put_yashxtt[3,:] = np.mean(a0,(1,2))
    out_put_yashxtt[4,:] = np.mean(a1,(1,2))
    out_put_yashxtt[5,:] = meancounts_average_flatfields
    out_put_yashxtt[6,:] = meancounts_average_flatfields_norm

"""
measurment over two days sponsored by friedrichs limestone flatfields
"""
if limestone:
    sam_ID ='NDT_FreshCon_Limestone_TimeScan_#2'
    path_to_data = '/data/DPC/local_setups/microfocus/samples/'+sam_ID+'/paxscan/'
    sub_dirs=os.listdir(path_to_data)
        
    flat_dirs=list()
    
    dir_num=1
    for sub_dir in sub_dirs:    
        dirs=np.sort(os.listdir(path_to_data+sub_dir))
        #One FlatField Dir, One Data Dir, One FlatField Dir, One Data Dir........
        for dir_name in dirs:
            if dir_num%2 !=0:         
                flat_dirs.append(path_to_data+sub_dir+'/'+dir_name)
            dir_num=dir_num+1
        
    relevant_folders = (len(flat_dirs)-39)
    deltas=np.empty((relevant_folders))
    deltas[0]=0.
    t = 0
    if processed:
        a0 = np.ones((relevant_folders,800,800))
        a1 = np.ones((relevant_folders,800,800))
        mean_visibility = np.ones(relevant_folders)
        #mean_visibility_norm = np.ones(relevant_folders)
        visibility = np.ones((relevant_folders,400,450))
        meancounts_average_flatfields = np.zeros((relevant_folders))
        for i in range(relevant_folders):
            print 'Foldernumber', i
            ff = np.ones((8,800,800))
            for j in range(8):
                ## for undertanding    ####################flatfieldfolder[i]||[j] flatfieldnumber
                ff[j] = e17.io.h5read (flat_dirs[i]+'/'+os.listdir(flat_dirs[i])[j]) ["raw_data"]#[200:650,150:600]
            coeff = dpcprocess.lsq_fit(ff, nb_periods = 1, order = 2)
            a0[i] = coeff[0,0]
            a1[i] = coeff[0,1]
        for i in range(relevant_folders):
            visibility[i]=a1[i][200:600,150:600]/a0[i][200:600,150:600]
            visibility[i][np.isnan(visibility[i])]=1.
            visibility[i][np.isinf(visibility[i])]=1.
            mean_visibility[i] = np.mean(visibility[i])
            meancounts_average_flatfields[i] =(np.mean(a0[i][5:100,5:170])+np.mean(a0[1][5:80,585:700])+np.mean(a0[i][728:795,17:150]))/3.
        mean_visibility_norm = mean_visibility/mean_visibility[0]
        meancounts_average_flatfields_norm = meancounts_average_flatfields/meancounts_average_flatfields[0]
    for i in range(relevant_folders):
        if i< relevant_folders-1:
            delta_t = os.stat(flat_dirs[i+1]).st_mtime - os.stat(flat_dirs[0]).st_mtime
            t=t+delta_t
            deltas[i+1]=delta_t/3600.
    out_put_limestone = np.zeros((5,np.shape(visibility)[0]))
    out_put_limestone[0,:] = deltas
    out_put_limestone[1,:] = mean_visibility
    out_put_limestone[2,:] = mean_visibility_norm
    out_put_limestone[3,:] = meancounts_average_flatfields
    out_put_limestone[4,:] = meancounts_average_flatfields_norm


    
"""
measurement of a snakehead by Florian
"""            
if snakehead:
    sam_ID ='microCT_SnakeHead'
    path_to_data = '/data/DPC/local_setups/microfocus/samples/'+sam_ID+'/paxscan/'
    sub_dirs=os.listdir(path_to_data)
        
    flat_dirs=list()
    
    dir_num=1
    for sub_dir in sub_dirs:    
        dirs=np.sort(os.listdir(path_to_data+sub_dir))
        #One FlatField Dir, One Data Dir, One FlatField Dir, One Data Dir........
        for dir_name in dirs:
            if dir_num%2!=0:         
                flat_dirs.append(path_to_data+sub_dir+'/'+dir_name)
            dir_num=dir_num+1
    ntrblocks = 52
    numflats = 1
    sizey = 900
    sizex = 1536
    #deltas=np.empty((ntrblocks))
    #deltas[0]=0.
    #t = 0
    
    ff = np.zeros((ntrblocks,sizey,sizex))
    for i in range(ntrblocks):
        filename = "/S%05d/paximage_%05d_%05d_00.h5"%(826 + 2*i,826 + 2*i,0)
        print "Reading file", filename
        ff[i] += e17.io.h5read(flat_dirs[i]+'/'+os.listdir(flat_dirs[i])[0])["raw_data"]
    meancounts_snake= np.mean(ff[:,:,:],axis = (1,2))
    meancounts_snake_flats= meancounts_snake/meancounts_snake[0]
    time_span = 4*3600+13*60+58
    time = np.asarray(np.linspace(0,time_span,52))
    time = time/3600.
    
    out_put_snakehead = np.zeros((3,np.shape(meancounts_snake)[0]))
    out_put_snakehead[0,:] = time
    out_put_snakehead[1,:] = meancounts_snake
    out_put_snakehead[2,:] = meancounts_snake_flats
if plots:
    if power:
        regression1 = np.polyfit(power_new,meancounts_og_mp,1)
        fit_og_mp = np.poly1d(regression1)
        regression2 = np.polyfit(power_new,meancounts_mg_mp,1)
        fit_mg_mp = np.poly1d(regression2)
        regression3 = np.polyfit(power_new,meancounts_og_mp_new,1)
        fit_og_mp_new = np.poly1d(regression3)
        regression4 = np.polyfit(power_new,meancounts_mg_mp_new,1)
        fit_mg_mp_new = np.poly1d(regression4)
        plt.figure('Meancounts vs Power')    
        plt.plot(power_new,meancounts_og_mp, 'b.')
        plt.plot(power_new,meancounts_mg_mp, 'g.')
        plt.plot(power_new,meancounts_og_mp_new,'k.')
        plt.plot(power_new,meancounts_mg_mp_new,'c.')
        plt.plot(power_new,fit_og_mp(power_new),'r',power_new,fit_mg_mp(power_new),'r',
                 power_new,fit_og_mp_new(power_new),'r',power_new,fit_mg_mp_new(power_new),'r', lw =.3 )
#        pyF.fit_linear(meancounts_og_mp, x_data = power_new,plot=False)
#        pyF.fit_linear(meancounts_mg_mp, x_data = power_new,plot=True)
#        pyF.fit_linear(meancounts_og_mp_new, x_data = power_new,plot=True)
#        pyF.fit_linear(meancounts_mg_mp_new, x_data = power_new,plot=True)
        plt.title("Meancounts vs Power")
        plt.xticks( np.arange(0, np.max(power_new),20) )
        plt.xticks(fontsize  = 8)
        plt.yticks( np.arange(0, np.max(meancounts_og_mp),20) )
        plt.yticks(fontsize  = 8)
        plt.ylim((0,100100.))
        plt.ylabel('Meancounts $[arb.unit]$',fontsize = 12)
        plt.xlabel('Power $[W]$',fontsize = 12)
        plt.legend(( 'results without gratings 22.10.2014','Results with gratings 22.10.2014',
        'Results without gratings 02.07.2015','results with gratings 02.07.215','Red lines = Regression of the data') ,loc = 0)    
        plt.savefig(filepath_p+'Meancounts vs Power with regression.tiff', format = 'tiff',dpi = 300 )   
    if exp_time:
        plt.figure('Meancounts vs exposure time')
        #plt.plot(time_s,meancounts_met,'bx')
        pyF.fit_linear(meancounts_met, x_data = time_s,plot=True)  
        plt.title("Meancounts vs exposure time") 
        plt.ylabel('Meancounts $[arb. unit]$')
        plt.xlabel('Time $[s]$')
        plt.legend((( 'Green dots = measured data'), ('Red line = Regression of the data')) ,loc = 0)    
        plt.savefig(filepath_p+'Meancounts vs Exposure time.tiff', format = 'tiff',dpi = 300 )
    if short_m:    
        plt.figure('Short timestability measurement')
        plt.plot(seconds_mps,meancounts_mps)
        plt.title("Meancounts vs Time short time scale") 
        plt.ylabel('Meancounts $[arb. unit]$')
        plt.xlabel('Time $[s]$')
        plt.legend(('Meancounts of each picture',) ,loc = 0)
        plt.savefig(filepath_p+'Meancounts vs Time short time scale.tiff', format = 'tiff',dpi = 300 )
    if long_m:
        plt.figure('Long timestability measurement')
        plt.plot(seconds_mpl,meancounts_mpl)
        plt.title("Meancounts vs Time long time scale") 
        plt.ylabel('Meancounts $[arb. unit]$')
        plt.xlabel('Time $[s]$')
        plt.legend(('Meancounts of each picture',) ,loc = 0)
        plt.savefig(filepath_p+'Meancounts vs Time long time scale.tiff', format = 'tiff',dpi = 300 )
    if oneweek:        
        plt.figure('Measurement over a hole week')
        plt.plot(deltas, mean_visibility, lw=1, ls = '--')
        plt.title('source time dependency')
        plt.ylabel('meancounts $[arb. unit]$')
        plt.xlabel('time $[h]$')
        plt.legend(('mean visibility',) ,loc = 0)
        plt.savefig(filepath_p+'Timedependency whole week.tiff', format = 'tiff',dpi = 300 )
    if snakehead:    
        plt.figure('Measurement snakehead')
        plt.plot(time, meancounts_snake_flats, lw=1, ls = ':')
        plt.title('source time dependency')
        plt.ylabel('meancounts $[arb. unit]$')
        plt.xlabel('time $[h]$')
        plt.legend(('meancounts of each picture',) ,loc = 0)
        #plt.savefig(filepath_p+'Timedependency while snakemeasurement.tiff', format = 'tiff',dpi = 300 )
        
if saveon:
    filename = 'stability_data'
    #deltas = np.swapaxes(deltas,0,1)
    #mean_visibility = np.swapaxes(mean_visibility,0,1)
    
    
    try:
        os.makedirs(filepath_p)    
    except:
        print('Folder already exist')
    try:
        np.save(filepath_p+filename+'power', out_put_power)
        np.save(filepath_p+filename+'mean_vs_exptime', out_put_exp_time)
        np.save(filepath_p+filename+'time_stability_measurement_281014', out_put_281014)
        np.save(filepath_p+filename+'time_stability_measurement_291014', out_put_291014)
        np.save(filepath_p+filename+'yashxttz_one_week',out_put_yashxtt)
        np.save(filepath_p+filename+'limestone_fritz',out_put_limestone)
        np.save(filepath_p+filename+'snakehead_flo',out_put_snakehead)
        #np.savetxt(filepath_d+filename+'.txt',((deltas, mean_visibility)),fmt = '%10.5f',delimiter = ', ',header = info)
        #np.savez(filepath_d+filename+'from_yashxtt2',deltas = 'time',mean_visibility = 'visibility',meancounts_average_flatfields ='meancounts',
        #         a0 = 'a0',a1 = 'a1',phi = 'phi')
    except:
        print('Error, could not save the file (does the script/terminal have write access?)')

#    '/users/Baier/analysis_files/stability_datafrom_yashxtt2.npz'
        