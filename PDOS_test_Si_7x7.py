#!/usr/bin/python

import os
import numpy as np

#import pyProbeParticle.GridUtils as GU
#import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import pyPPSTM.PreSTMutils       as SU

import matplotlib
#matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties


# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(0,0)		# just try
lvs = None      # automatically taken from geometry.in
WorkFunction = 5.0 #more or less standart.
fermi=None	# None - all energies relative to the Fermi
orbs= 'spd'	# 'sp' works now, 'spd' works as well
cut_min=-2.5	# 
cut_max=+2.5	# 
cut_at=18	# All atoms of the highest layer
eta = 0.1	# very low, to pronounce the single orbitals only
# -- next two not needed for IETS
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = 'phik_0001_', geom='input.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
                        cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- on which energies you want to plot pseudo-projected density of states

energies = np.arange(-2.0,2.0,0.01)

# --- getting P-PDOS for different atoms and spherical functions

#PDOS0 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[], orbs=orbs ,spherical='all') # all atoms = []
PDOS1 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[0], orbs=orbs ,spherical='all')
PDOS2 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[0], orbs=orbs ,spherical='dz2')
PDOS3 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[5], orbs=orbs ,spherical='all')
PDOS4 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[5], orbs=orbs ,spherical='pz')
PDOS5 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[6], orbs=orbs ,spherical='all')
PDOS6 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[6], orbs=orbs ,spherical='pz')
PDOS7  = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[10], orbs=orbs ,spherical='all')
PDOS8  = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[10], orbs=orbs ,spherical='pz')
PDOS9  = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[15], orbs=orbs ,spherical='all')
PDOS10 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[15], orbs=orbs ,spherical='pz')
PDOS11 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[14], orbs=orbs ,spherical='all')
PDOS12 = SU.pPDOS(eigEn, coefs, energies, eta=eta, atoms=[14], orbs=orbs ,spherical='pz')

# --- plotting P-PDOS for different atoms and spherical functions

#plt.plot(energies,PDOS1,ls='-',c='k',label='PDOS whole surface')
plt.plot(energies,PDOS1,ls='-' ,c='r',label='PDOS FC adatom')
plt.plot(energies,PDOS2,ls='--',c='r',label='PDOS FC adatom - pz orbital only')
plt.plot(energies,PDOS3,ls='-' ,c='m',label='PDOS FE adatom')
plt.plot(energies,PDOS4,ls='--',c='m',label='PDOS FE adatom - pz orbital only')
plt.plot(energies,PDOS5,ls='-' ,c='b',label='PDOS UC adatom')
plt.plot(energies,PDOS6,ls='--',c='b',label='PDOS UC adatom - pz orbital only')
plt.plot(energies,PDOS7,ls='-' ,c='c',label='PDOS UE adatom')
plt.plot(energies,PDOS8,ls='--',c='c',label='PDOS UE adatom - pz orbital only')
plt.plot(energies,PDOS9,ls='-' , c='g',label='PDOS F restatom')
plt.plot(energies,PDOS10,ls='--',c='g',label='PDOS F restatom - pz orbital only')
plt.plot(energies,PDOS11,ls='-' ,c='y',label='PDOS U restatom')
plt.plot(energies,PDOS12,ls='--',c='y',label='PDOS U restatom - pz orbital only')
plt.title("PDOS eta: "+str(eta)+' eV')
plt.xlabel("E-Ef [eV]")
plt.ylabel("DOS [arb.un.]")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.axis('image')
plt.savefig( 'PDOS_eta_'+str(eta)+'eV.png' , bbox_inches='tight' ) 
#plt.show()
plt.close()

# --- the end

print() 
print()
print("Done")
