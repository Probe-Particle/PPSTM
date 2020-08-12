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
cut_min=-1.0	#cut_min=-2.5	# +/- 0.5 eV bigger than the maximal bias (+/- 2.0 V)
cut_max=+1.0	#cut_max=+2.5	# <- older from older script
cut_at=18		# Adatoms (12) and restatoms (6) only
eta = 0.1		# 0.1 - standart for semiconductors
# -- next two parameters will be specified in the main loop
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9			# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings

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

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

#Voltages=[-0.88,+0.88]
#namez=['HOMO','LUMO']

#Voltages=np.arange(-2.0,+2.0+0.01,0.5) # this part is important for scans over slabs at different voltages
# ABOVE - original long calculations, 
Voltages=np.arange(-0.5,+0.5+0.01,0.5) # this part is important for scans over slabs at different voltages
namez = []
for V in Voltages:
    namez.append(str(round(V,1)))

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

curr0 = np.array([])
curr1 = np.array([])
curr2 = np.array([])

for WorkFunction in [WorkFunction]:
    i=0;
    for V in Voltages:
        if (V != 0.0):
            nV = int(abs(V/eta)+1)
            print("Voltage: ", V, "; number of integration steps: ", nV)
            for eta in [eta]:
                current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0) 
                WF_decay = 1 # ( tunnelling barrier changes 1:1 with applied sample bias)
                current1 = PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, WF_decay=WF_decay)
                WF_decay = 0 # ( tunnelling barrienr doesn't change)
                current2 = PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, WF_decay=WF_decay)
                curr0 = np.append(curr0,current0)
                curr1 = np.append(curr1,current1)
                curr2 = np.append(curr2,current2)
    i = i+1

curre=np.reshape(curr0,(len(Voltages)-1,sh[0],sh[1],sh[2]))
# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='didV-7x7-5Ang.dat'
# ploting part here:
plt.figure( figsize=(0.2* xl , 0.2*yl/2 ) )
for i in range(2):
        j=i if i <=0 else i+1
        plt.subplot(1,2,i+1)
        plt.imshow( curre[i,0,:,:], origin='image', extent=extent , cmap='gray')
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title("Sample bias:"+namez[j]+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

curre=np.reshape(curr1,(len(Voltages)-1,sh[0],sh[1],sh[2]))
# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='STM_corresponding_barrier-7x7-5Ang.dat'
# ploting part here:
plt.figure( figsize=(0.2* xl , 0.2*yl/2 ) )
for i in range(2):
        j=i if i <=0 else i+1
        plt.subplot(1,2,i+1)
        plt.imshow( curre[i,0,:,:], origin='image', extent=extent , cmap='gray')
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title("Sample bias:"+namez[j]+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

curre=np.reshape(curr2,(len(Voltages)-1,sh[0],sh[1],sh[2]))
# --- plotting part here, plots calculated signal:
print(" plotting ")
name_file='STM_same_barrier-7x7-5Ang.dat'
# ploting part here:
plt.figure( figsize=(0.2* xl , 0.2*yl/2 ) )
for i in range(2):
        j=i if i <=0 else i+1
        plt.subplot(1,2,i+1)
        plt.imshow( curre[i,0,:,:], origin='image', extent=extent , cmap='gray')
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title("Sample bias:"+namez[j]+"V")

plt.savefig( name_file+'.png', bbox_inches='tight' )
#plt.show()
plt.close()

# --- the end

print() 
print()
print("Done")
