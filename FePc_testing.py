#!/usr/bin/python

import os
import numpy as np

import pyProbeParticle.GridUtils as GU
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS

import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
import timeit

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
#cut_min=[3.15897074,2.27956401,3.15461612,-1.80707498,-2.29990905]	# test only
#cut_min+=-0.0001
#cut_max=[3.15897074,2.27956401,3.15461612,-1.80707498,-2.29990905]	# test_only
#cut_max+=+0.0001
cut_min=-1.81
cut_max=+3.2
cut_at=57	# All atoms of the molecule (MetalPc - 57 atoms, H2Pc - 58, NoPc - 56)
eta = 1e-8		# test
# --- these two not needed now (no STM in this script)
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

#eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
#					    cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
eigEn1, coefs1, Ratin = RS.read_AIMS_all(name = 'KS_exx_1_spin_up.out', geom='geometry.in',fermi=fermi, orbs = orbs, pbc=pbc,
					imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
					lower_atoms=lower_atoms, lower_coefs=lower_coefs)
eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = 'KS_exx_1_spin_dn.out', geom='geometry.in',fermi=fermi, orbs = orbs, pbc=pbc,
					imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
					lower_atoms=lower_atoms, lower_coefs=lower_coefs)
eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
print "eigEn: ", eigEn
coefs = np.concatenate((coefs1, coefs2), axis=0)


#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; no tip_r1 - PP distored by the relaxation in the PPAFM code;  only tip_r2 - uniform grid:

#tip_r1, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)

dz=0.25
dx=dy =0.25

xstart = -12.0
ystart = -12.0
zstart =  -7.0  # the purely theoretical scan showing appearance of the original orbitals0 at height - 1 - 5 Angstrom above the molecule

xl = 24.0
yl = 24.0
zl = 14.0
extent = (xstart,xstart+xl,ystart,ystart+yl)

tip_r2 = RS.mkSpaceGrid(xstart,xstart+xl,dx,ystart,ystart+yl,dy,zstart,zstart+zl,dz)
sh = tip_r2.shape

lvec= np.array([[0.0,0.0,0.0],[xl,0.,0.],[0.,yl,0.],[0.,0.,zl]])
# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

Voltages=[3.15897074,2.27956401,3.15461612,-1.80707498]
namez=['148_dn_test2',"147_dn_test2","149_up_test2","145_up_test2"]

#Voltages=np.arange(-1.0,+1.0+0.01,0.1) # this part is important for scans over slabs at different voltages
#namez = []
#for V in Voltages:
#    namez.append(str(round(V,1)))

# --- downloading the df data

#df, lvec2, nDim2 = GU.load_scal_field( path_df+'df' ,data_format=data_format)

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

for WorkFunction in [WorkFunction]:
	i=0;
	for V in Voltages:
		for eta in [eta]:
			s = timeit.timeit()
			current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
			e = timeit.timeit()
			print "time in c:", e-s
			GU.saveXSF(namez[i]+".xsf", current0, lvec )
		i += 1
# --- the end

print 
print
print "Done"
