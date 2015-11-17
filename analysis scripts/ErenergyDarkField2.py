import numpy as np
import scipy.constants as cons
import re
from scipy import ndimage
import csv
import matplotlib.pyplot as plt
import ddfSetupProcessing.HarmonicAnalysis.lsqprocessing as dpcprocess
ff_file=open('/data/DPC/local_setups/microfocus/CdTEspcectra/visimapzentral120kv10w30s512chan/visimapzentral_120kv10w30s.csv')
im_file=open('/data/DPC/local_setups/microfocus/CdTEspcectra/SiO2_7.38_settled_120kv10w30s512chan/SiO2_7.38_settled_120kv10w30s512chan.csv')

ff_csv=csv.reader(ff_file)
im_csv=csv.reader(im_file)

ff_data=[]
im_data=[]
p_two=14.54e-6
d=0.6-0.03

for row in ff_csv:
    ff_data.append(row[0])  
    
for row in im_csv:
    im_data.append(row[0])      
    
    
cal_ch_a=float(ff_data[7][0:3])
cal_en_a=float(ff_data[7][4:8])

cal_ch_b=float(ff_data[8][0:3])
cal_en_b=float(ff_data[8][4:8])

ch_en=(cal_en_b-cal_en_a)/(cal_ch_b-cal_ch_a)
ch_zero=cal_en_a-ch_en*cal_ch_a
ch_corr=np.empty((526-14,1))

ff=np.empty((8,526-14,1))
im=np.empty((8,526-14,1))
energies=np.empty((526-14,1))

for j in np.arange(0,526-14,1):
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

smooth_ff=np.empty((8,526-14,1))
smooth_im=np.empty((8,526-14,1))
            
for i in np.arange(0,8,1):  
    smooth_ff[i]=ndimage.uniform_filter(ff[i],size=5)
    smooth_im[i]=ndimage.uniform_filter(im[i],size=5)
        
coeff_ff = dpcprocess.lsq_fit(smooth_ff, nb_periods = 1, order = 1)        
coeff_im = dpcprocess.lsq_fit(smooth_im, nb_periods = 1, order = 1) 
               
atc, dpc, dfc = dpcprocess.dpcprocess(smooth_ff, smooth_im, nb_periods = 1, order = 1)
figure()
plt.plot(energies,coeff_ff[0,1]/coeff_ff[0,0])        
figure()
plt.plot(energies,dfc,'r+')    
figure()
plt.plot(ch_corr,dfc,'r+')               
figure()
plt.plot(energies,smooth_ff[0,:,0])             
            