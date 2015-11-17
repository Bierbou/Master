import numpy as np
import scipy.constants as cons
import re
from scipy import ndimage
import csv
import matplotlib.pyplot as plt
import ddfSetupProcessing.HarmonicAnalysis.lsqprocessing as dpcprocess
ff_file=open('/users/Prade/spectra/spectra/visimapzentral100kv100W60s/visimap comp.csv')
im_file=open('/users/Prade/spectra/spectra/sio2_2.79microspheresshaked/SIO2_2.79shakecomp.csv')

ff_csv=csv.reader(ff_file)
im_csv=csv.reader(im_file)

ff_data=[]
im_data=[]
p_two=10e-6
d=0.925-0.04

for row in ff_csv:
    ff_data.append(row[0])  
    
for row in im_csv:
    im_data.append(row[0])      
    
    
cal_ch_a=float(ff_data[7][0:4])
cal_en_a=float(ff_data[7][5:9])

cal_ch_b=float(ff_data[8][0:4])
cal_en_b=float(ff_data[8][5:9])

ch_en=(cal_en_b-cal_en_a)/(cal_ch_b-cal_ch_a)
ch_zero=cal_en_a-ch_en*cal_ch_a
ch_corr=np.empty((8206-14,1))

ff=np.empty((8,8206-14,1))
im=np.empty((8,8206-14,1))
energies=np.empty((8206-14,1))

for j in np.arange(0,8206-14,1):
    ff_line=map(int,re.findall(r'\d+',ff_data[j+14]))
    im_line=map(int,re.findall(r'\d+',im_data[j+14]))
    for i in np.arange(0,8,1):    
        ff[i,j,0]=float(ff_line[2*i+1])
        im[i,j,0]=float(im_line[2*i+1])
        energies[j]=j*ch_en+ch_zero
        
        if j==0:
            ch_corr[j]=100000000
        else:    
            ch_corr[j]=(d*cons.h*cons.c)/(1000.*energies[j]*cons.e*p_two)
            
smooth_ff=np.empty((8,8206-14,1))
smooth_im=np.empty((8,8206-14,1))
for i in np.arange(0,8,1):  
    smooth_ff[i]=ndimage.uniform_filter(ff[i],size=100)
    smooth_im[i]=ndimage.uniform_filter(im[i],size=100)
        
coeff_ff = dpcprocess.lsq_fit(smooth_ff, nb_periods = 1, order = 1)        
coeff_im = dpcprocess.lsq_fit(smooth_im, nb_periods = 1, order = 1) 
               
atc, dpc, dfc = dpcprocess.dpcprocess(smooth_ff, smooth_im, nb_periods = 1, order = 1)
plt.figure()
plt.plot(energies,coeff_ff[0,1]/coeff_ff[0,0])        
plt.figure()
plt.plot(ch_corr,dfc)            
            
            