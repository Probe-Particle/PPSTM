#!/usr/bin/python

import os
import numpy as np

#import pyProbeParticle.GridUtils as GU
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS

import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

# --- specification of paths to the STM input files, PP positions and stored df results & format in which the data are stored - xsf or npy

path=''
# --- not needed at all for the rigid tip scan
#path_pos='Q0.00K0.24/'
#path_df = path_pos+'Amp1.00/'
#data_format ="npy"

# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(0,0)
lvs = None      # automatically taken from geometry.in (however not needed since 
WorkFunction = 5.0 #more or less standart.
fermi=None	# the Fermi from AIMS automatically 0 -- None means no shift!!! All energies are relative to Fermi !!!!
orbs= 'spd'	# 'sp' & 'spd' works now
cut_min=[-1.81,-0.52,1.76]	# To cut-out all the orbitals axcept for the wanted once
cut_max=[-0.95,-0.43,1.85]	# (HOMO,SOMO,LUMO1,LUMO2); their energies are at Voltages
cut_at=57	# All atoms of the molecule (MetalPc - 57 atoms, H2Pc - 58, NoPc - 56)
eta = 0.1	# very very high, to pronounce the together orbital symmetry
# --- these two not needed now (no STM in this script)
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

#eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
#					    cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);


#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; no tip_r1 - PP distored by the relaxation in the PPAFM code;  only tip_r2 - uniform grid:

#tip_r1, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)

dz=1.0
dx=dy =0.1

xstart = -12.0
ystart = -12.0
zstart =  1.0  # the purely theoretical scan showing appearance of the original orbitals0 at height - 1 - 5 Angstrom above the molecule

xl = 24.0
yl = 24.0
zl = 4.0
extent = (xstart,xstart+xl,ystart,ystart+yl)

tip_r2 = RS.mkSpaceGrid(xstart,xstart+xl,dx,ystart,ystart+yl,dy,zstart,zstart+zl,dz)
sh = tip_r2.shape

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

Voltages=[-1.38,-0.48,1.805]
namez=['HOMO-1','HOMO','LUMO']

#Voltages=np.arange(-1.0,+1.0+0.01,0.1) # this part is important for scans over slabs at different voltages
#namez = []
#for V in Voltages:
#    namez.append(str(round(V,1)))

# --- downloading the df data

#df, lvec2, nDim2 = GU.load_scal_field( path_df+'df' ,data_format=data_format)

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

curr0 = np.array([])
curr1 = np.array([])

for WorkFunction in [WorkFunction]:
	i=0;
	for V in Voltages:
		eigEn1, coefs1, Ratin = RS.read_AIMS_all(name = 'KS_exx_1_spin_up.out', geom='geometry.in',fermi=fermi, orbs = orbs, pbc=pbc,
												imaginary = False, cut_min=cut_min[i], cut_max=cut_max[i], cut_at=cut_at,
												lower_atoms=lower_atoms, lower_coefs=lower_coefs)
		eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = 'KS_exx_1_spin_dn.out', geom='geometry.in',fermi=fermi, orbs = orbs, pbc=pbc,
												imaginary = False, cut_min=cut_min[i], cut_max=cut_max[i], cut_at=cut_at,
												lower_atoms=lower_atoms, lower_coefs=lower_coefs)
		
		print Ratin
		#exit()

		eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
		print "eigEn: ", eigEn
		coefs = np.concatenate((coefs1, coefs2), axis=0)
		for eta in [eta]:
			current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
			current1 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=0.0, px=0.5, py=0.5, pz = 0.0)
			curr0 = np.append(curr0,current0)
			curr1 = np.append(curr1,current1)
		i += 1

curre=np.reshape(curr0,(len(Voltages),sh[0],sh[1],sh[2]))
curre1=np.reshape(curr1,(len(Voltages),sh[0],sh[1],sh[2]))
curre2=curre+curre1

# --- plotting part here, plots calculated signal:
print " plotting "
for zi in range(len(tip_r2)):
	name_file='didV-FePc-d_orb_decay_slater_height_%03dA' %(zi+1)
	# ploting part here:
	plt.figure( figsize=(0.45* xl , 0.45*yl ) )
	for i in range(3):
		plt.subplot(3,3,i+1)
		plt.imshow( curre[i,zi,:,:], origin='image', extent=extent , cmap='gray')
		plt.ylabel(r' Tip_y $\AA$'+"_s_tip")
		plt.title("dI/dV "+namez[i])
		plt.subplot(3,3,i+4)
		plt.imshow( curre1[i,zi,:,:], origin='image', extent=extent , cmap='gray')
		plt.ylabel(r' Tip_y $\AA$'+"_pxy_tip")
		plt.subplot(3,3,i+7)
		plt.imshow( curre2[i,zi,:,:], origin='image', extent=extent , cmap='gray')
		plt.ylabel(r' Tip_y $\AA$'+"_sp_tip")
		plt.xlabel(r' Tip_x $\AA$')
	plt.savefig( name_file+'.png', bbox_inches='tight' )
	plt.show()
	plt.close()


# --- the end

print 
print
print "Done"
