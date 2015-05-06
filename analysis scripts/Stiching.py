# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 14:33:38 2015

@author: gu32hej
"""

path=''
import glob
import pyE17 as e17
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
files = glob.glob(path+'*.h5')
files.sort()
flat=e17.io.h5read(files[0])['raw_data']
a=e17.io.h5read(files[1])['raw_data']
b=e17.io.h5read(files[2])['raw_data']
a=a/(flat*1.)
b=b/(flat*1.)
#plt.imshow(flat,vmin=0)
#plt.figure('a')
#plt.imshow(a,vmin=0)
#plt.figure('b')
#plt.imshow(b,vmin=0)

#crop_a,coords_a=e17.utils.Interactive_rect_roi(a,False,False)
#crop_b,coords_b=e17.utils.Interactive_rect_roi(b,False,False)
c=a[800-40:,:]

d=b[:40,:]
#R=e17.register.Register()
#template=e17.register.RegisterData(crop_a)
#image=e17.register.RegisterData(crop_b)
#p,s=R.register(image,template,e17.register.model.Shift(),sampler = e17.register.sampler.spline)

distances,trans_vec,alpha = e17.utils.shift_best(d,c)

plt.figure('c')
plt.imshow(c)
plt.figure('d')
plt.imshow(d)
plt.figure('dist')
plt.imshow(distances)

f = a[800-25:,:]
g = b[:25,:]

h = (f+g)/2.
pic2_test25 = np.vstack((a[:-25,:],h,b[25:,:]))
 
plt.figure('pic2_test25')
plt.imshow(pic2_test25)



