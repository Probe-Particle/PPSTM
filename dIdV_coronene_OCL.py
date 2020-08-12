#!/usr/bin/python

import os
import numpy as np
import sys
sys.path.append('./')

import pyProbeParticle.GridUtils as GU
import pyPPSTM.ReadSTM           as RS
import pyPPSTM.OCL               as ocl

import time

# ======= setup

path  = './'

decay = -1.0    # WARRNING ... with positive decay you get bullshit

# ======= Prepare inputs

eigEn, coefs, Ratin = RS.read_FIREBALL_all(   
    name = path+'phik_example_', geom=path+'crazy_mol.xyz', 
    fermi=None, orbs='sp', pbc=(0,0),
    cut_min=-1.0, cut_max=+1.0, cut_at=-1, 
    lower_atoms=[], lower_coefs=[]
);

print("---------------")
print("Ratin ", Ratin)
print("---------------")
print("coefs ", coefs)
print("---------------")
print("eigEn.shape ", eigEn.shape)
print("coefs.shape ", coefs.shape)
print("Ratin.shape ", Ratin.shape)

atoms    = ocl.xyzq2float4( Ratin, np.ones(len(Ratin))*decay )   # 4th component per atom is decay, can be used to introduce effect of different local potential barrier
CAOs     = ocl.CAOsp2f4   (coefs, len(atoms) )
spectral = ocl.getSpectral( eigEn, Wf = 1.0, w=0.2 )

print("---------------")
print("CAOs", CAOs)
        
lvec=np.array([
    (0.0 ,0.0 ,7.0),
    (16.0,0.0 ,0.0),
    (0.0 ,16.0,0.0),
    (0.0 ,0.0 ,6.0)
])
rTips  = ocl.getPos_f4( lvec )

# ======= Run

t1 = time.clock()
kargs = ocl.initArgs(atoms, CAOs, spectral, rTips)
Gout  = ocl.run( kargs, rTips.shape[:3] )
t1 = time.clock() - t1

print(" Run time %g [s] " %t1) 

# =====================
# =====================  Plotting
# =====================

Ftmp = np.zeros(rTips.shape[:3]); 
Ftmp[:,:,:] = Gout[:,:,:]
#Ftmp[:,:,:] = rTips[:,:,:,0]

GU.saveXSF(path+"G_s_sp_ocl.xsf", Ftmp, lvec )

'''
import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

extent  = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[2,1])
'''
    
print("===== ALL DONE ==== ")

