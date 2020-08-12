#!/usr/bin/python

import os
import numpy as np

#import pyProbeParticle.GridUtils as GU # --- not needed at all for the rigid tip scans
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS

import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

# --- specification of paths to the STM input files, PP positions and stored df results & format in which the data are stored - xsf or npy

path=''
# --- not needed at all for the rigid tip scan
#path_pos='Q0.00K0.50/'
#path_df = path_pos+'Amp0.40/'
#data_format ="npy"

# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(1,1)		# 3x3 cell shifted so the original cell is in the middle
lvs = np.loadtxt(path+'input.lvs')	
#lvs = np.array([[23.661884370846174,-13.661203825532905,0.],[23.661884370846174,13.661203825532905,0.],[0.,0.,99.]]
WorkFunction = 5.0 #more or less standart.
fermi=None		# the Fermi from phik ... .dat file; !!! All energies are relative to Fermi !!!!
orbs= 'spd'		# 'sp' and 'spd' works now
cut_min=-2.5	# +/- 0.5 eV bigger than the maximal bias (+/- 2.0 V)
cut_max=+2.5	# 
cut_at=18		# Adatoms (12) and restatoms (6) only
eta = 0.1		# 0.1 - standart for semiconductors
# -- next two parameters will be specified in the main loop
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9			# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings
Vmin = -2.0	# Minimum for STM & dIdV scan
Vmax = +2.0	# Maximum for STM & dIdV scan
dV   =  0.1

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_0001_', geom=path+'input.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
                        cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lvs=lvs, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin = RS.read_AIMS_all(name = 'KS_eigenvectors_up.band_1.kpt_1.out', geom='geometry.in',fermi=fermi, orbs = 'sp', pbc=pbc,
#					imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
#					lower_atoms=lower_atoms, lower_coefs=lower_coefs)
#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; tip_r2 - uniform grid:

dz=0.1
dx=dy =0.25

xstart = -2.0
ystart = -15.0
zstart = 9.7  # the highest atom is at 4.7, this is only 5Ang above scan

xl = 50.0
yl = 30.0
zl = 0.01
extent = (xstart,xstart+xl,ystart,ystart+yl)

tip_r2 = RS.mkSpaceGrid(xstart,xstart+xl,dx,ystart,ystart+yl,dy,zstart,zstart+zl,dz)
sh = tip_r2.shape

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

STM0 = np.array([])
STM1 = np.array([])
dIdV0 = np.array([])
dIdV1 = np.array([])

for WorkFunction in [WorkFunction]:
    WF_decay = 0 # ( tunnelling barrienr doesn't change)
    STM0, dIdV0 = PS.MSTM( Vmin, Vmax, dV, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, WF_decay=WF_decay)
    WF_decay = 1 # ( tunnelling barrier changes 1:1 with applied sample bias)
    STM1, dIdV1 = PS.MSTM( Vmin, Vmax, dV, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, WF_decay=WF_decay)

# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='didV-7x7-5Ang'
# ploting part here:
plt.figure( figsize=(0.4* xl , 0.4*yl/2 ) )
ni=21; nn=2;
for i in range(ni):
        plt.subplot(3,7,i+1)
        plt.imshow( dIdV0[i*nn,0,:,:], origin='image', extent=extent , cmap='gray')
        if np.mod(i,7)==0 :
            plt.ylabel(r' Tip_y $\AA$')
        if i >= 14:
            plt.xlabel(r' Tip_x $\AA$')
        plt.title("Sample bias:"+str(round(Vmin+i*nn*dV,1))+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='didV-7x7-5Ang_test'
# ploting part here:
plt.figure( figsize=(0.4* xl , 0.4*yl/2 ) )
for i in range(ni):
        plt.subplot(3,7,i+1)
        plt.imshow( dIdV1[i*nn,0,:,:], origin='image', extent=extent , cmap='gray')
        if np.mod(i,7)==0 :
            plt.ylabel(r' Tip_y $\AA$')
        if i >= 14:
            plt.xlabel(r' Tip_x $\AA$')
        plt.title("Sample bias:"+str(round(Vmin+i*nn*dV,1))+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='STM-7x7-5Ang_no_changes_in_WF'
# ploting part here:
plt.figure( figsize=(0.4* xl , 0.4*yl/2 ) )
for i in range(ni):
        plt.subplot(3,7,i+1)
        plt.imshow( STM0[i*nn,0,:,:], origin='image', extent=extent , cmap='gray')
        if np.mod(i,7)==0 :
            plt.ylabel(r' Tip_y $\AA$')
        if i >= 14:
            plt.xlabel(r' Tip_x $\AA$')
        plt.title("Sample bias:"+str(round(Vmin+i*nn*dV,1))+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='STM-7x7-5Ang_WF_changes_w_V'
# ploting part here:
plt.figure( figsize=(0.4* xl , 0.4*yl/2 ) )
for i in range(ni):
        plt.subplot(3,7,i+1)
        plt.imshow( STM1[i*nn,0,:,:], origin='image', extent=extent , cmap='gray')
        if np.mod(i,7)==0 :
            plt.ylabel(r' Tip_y $\AA$')
        if i >= 14:
            plt.xlabel(r' Tip_x $\AA$')
        plt.title("Sample bias:"+str(round(Vmin+i*nn*dV,1))+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()


# --- the end

print() 
print()
print("Done")
