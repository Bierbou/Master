# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 13:50:27 2012

@author:Guillaume Potdevin

        Edited by Chabior 22.08.2012 
        updated to microfocus setup 18.08.2014 Florian Schaff
"""

import pyE17 as e17
import numpy as np
#import pyqtgraph as pg
import matplotlib.pylab as pp
import os


process_sample = True
nb_periods = 1#0.95
saveon = True

filepath = "/data/DPC/local_setups/microfocus/samples/noG0_EM_stepping_bean/ct/"

ffs = np.zeros((19,800,770))
for i in range(19):
    ffs[i] = e17.io.h5read(filepath+"paximage_ct_4809%02d.h5"%(i+8))["raw_data"]


sort = np.array([18,17,16,15,14,4,3,2,1,0,5,6,7,8,9,10,11,12,13])
ffs = ffs[sort]

#pg.image(ffs.swapaxes(1,2))



import ddfSetupProcessing.HarmonicAnalysis.lsqprocessing as dpcprocess
import ddfSetupProcessing.PreprocessingTools as ppt

coeff = dpcprocess.lsq_fit(ffs, nb_periods = nb_periods, order = 2)

a0=coeff[0,0]
a1=coeff[0,1]
phi=coeff[1,1]
a0[a0<=1e-6]=1e-6
V=a1/a0

############################################################################
## Make Figures

mask = ppt.round_mask(coeff[0,0].shape, coeff[0,0].shape[0]*0.35)
roi_mask = 1-(mask-ppt.round_mask(coeff[0,0].shape, coeff[0,0].shape[0]*0.345))

a0=coeff[0,0]
a1=coeff[0,1]
phi=coeff[1,1]

a0[a0==0]=1
a0[np.isnan(a0)]=1
a1[np.isnan(a1)]=1
a0[np.isinf(a0)]=1
a1[np.isinf(a1)]=1

pp.figure(num=None, figsize=(10, 12.5), dpi=80, facecolor='w', edgecolor='k')
mean_int = a0.mean()
pp.subplot(3,2,1)
pp.imshow(a0, vmin=mean_int*0.65, vmax=mean_int*1.35, interpolation='none')
pp.colorbar()
pp.title('Flat Field Intensity')
pp.jet()
pp.grid()

pp.subplot(3,2,2)
pp.imshow(a1, interpolation='none',vmin=0,vmax=3000)      
pp.colorbar()
pp.title('Flat Field Amplitude')
pp.grid()

pp.subplot(3,2,3)
pp.imshow(phi, vmax=np.pi, vmin=-np.pi, interpolation='none')
cbar=pp.colorbar()
cbar.set_ticks([-np.pi, -0.75*np.pi, -0.5*np.pi, -0.25*np.pi, 0, 0.25*np.pi, 0.5*np.pi, 0.75*np.pi, np.pi])
cbar.set_ticklabels(['-pi', '-pi3/4', '-pi/2', '-pi/4', '0 pi', 'pi/4', 'pi/2', 'pi3/4', 'pi'])
pp.title('Flat Field Phase')
pp.grid()

pp.subplot(3,2,4)
mean_vis = (a1/a0)[mask].mean()
pp.imshow(a1/a0*roi_mask, vmin=mean_vis*0.5, vmax=mean_vis*1.5, interpolation='none')      
pp.colorbar()
pp.title('Flat Field Visibility')
pp.grid()
mean_vis_str = str(mean_vis*100)[:5] + '%'
pp.text(20,50, 'mean visibility: ' + mean_vis_str, bbox=dict(facecolor='white', alpha=1))
   
pp.subplot(3,1,3)
point_coord = [a0.shape[0]/2,a0.shape[1]/2]
steps = ffs[:,point_coord[0], point_coord[1] ]
dsteps = np.hstack((steps,steps))
N=np.linspace(1,ffs.shape[0]*2,ffs.shape[0]*2)
pp.plot(N,dsteps, 'bo--',label="measured")

Np=np.linspace(0.5,ffs.shape[0]+0.5,ffs.shape[0]*100)
Np2=np.linspace(0.5,ffs.shape[0]*2+0.5,ffs.shape[0]*2*100)
a0c = coeff[0,0,point_coord[0],point_coord[1]]
a1c = coeff[0,1,point_coord[0],point_coord[1]]
phc = coeff[1,1,point_coord[0],point_coord[1]]
curve1 = a0c + a1c * np.sin(2*np.pi*(Np-1)/(ffs.shape[0]) + phc + np.pi/2)
curve = np.hstack((curve1,curve1))
pp.plot(Np2,curve, 'r-',label="fitted")
pp.hlines(coeff[0,0,point_coord[0],point_coord[1]],0, ffs.shape[0]*2, 'k','--')
pp.vlines(ffs.shape[0]+0.5, coeff[0,0,point_coord[0],point_coord[1]]-coeff[0,1,point_coord[0],point_coord[1]],coeff[0,0,point_coord[0],point_coord[1]]+coeff[0,1,point_coord[0],point_coord[1]],'k','--')
pp.title('stepping curve of pixel (' + str(point_coord[0]) + ',' + str(point_coord[1]) + '). The curve is displayed twice.')
pp.ylabel('intensity [adu.]')
pp.xlabel('phase steps')
pp.legend()
pp.grid()

pp.show()


if process_sample:
    datm = np.zeros((19,800,770))
    for i in range(19):
        datm[i] = e17.io.h5read(filepath+"paximage_ct_4809%02d.h5"%(i+27))["raw_data"]
    
    
    # Compute basic parameters
    AMP, DPC, DCI = dpcprocess.dpcprocess(ffs, datm, nb_periods = nb_periods, order = 1)
    #DPC = e17.utils.rmphaseramp(DPC, mask)
    
    # Make Figures
    pp.figure(2)
    pp.imshow(AMP)
    pp.colorbar()
    pp.title('Absorption Signal')
    
    pp.figure(3)
    pp.imshow(DPC)
    pp.colorbar()
    pp.title('Phase Contrast Signal')
    #
    pp.figure(4)
    pp.imshow(DCI,vmax=1)
    pp.colorbar()
    pp.title('Dark Field Signal')

    if saveon:
        import PIL
        AMP_img=PIL.Image.fromarray(AMP)
        DPC_img=PIL.Image.fromarray(DPC)
        DCI_img=PIL.Image.fromarray(DCI)    
        AMP_img.save(os.getcwd()+'/AMP.tif')
        DPC_img.save(os.getcwd()+'/DPC.tif')
        DCI_img.save(os.getcwd()+'/DCI.tif')



















