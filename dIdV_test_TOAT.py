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
path_pos='Q0.00K0.50/'
path_df = path_pos+'Amp0.40/'
data_format ="npy"

# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(0,0)
lvs = None      # automatically taken from geometry.in (however not needed since 
WorkFunction = 5.0 #more or less standart.
fermi=None	# the Fermi from AIMS automatically 0 -- None means no shift!!! All energies are relative to Fermi !!!!
orbs= 'sp'	# 'sp' & 'spd' works now
cut_min=-0.6	# bellow HOMO-3
cut_max=+2.6	# above LUMO+1
cut_at=-1	# All atoms of the molecule (MetalPc - 57 atoms, H2Pc - 58, NoPc - 56)
eta = 0.01	# very low to pronounce single orbitals only
# --- these two not needed now (no STM in this script)
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[22,23,24]		# atoms 23-25 - oxygens will have lowered tunneling  !!! python numbering of atoms !!!
lower_coefs=[0.75,0.75,0.75]	# Lowering of the hoppings # NOTE: it seems, that for the CP2K input there is not suc a big need to lower down the contributions from
                # oxygens (on a TOAT), then from the other codes - 0.5

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

#eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
#					    cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin = RS.read_AIMS_all(name = 'KS_exx_1_spin_up.out', geom='geometry.in',fermi=fermi, orbs = orbs, pbc=pbc,
#					imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
#					lower_atoms=lower_atoms, lower_coefs=lower_coefs)
eigEn, coefs, Ratin  = RS.read_CP2K_all(name = 'TOAT', fermi=fermi, orbs = orbs, pbc=pbc,
                                        cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; no tip_r1 - PP distored by the relaxation in the PPAFM code;  only tip_r2 - uniform grid:

tip_r1, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)

dz=0.1
dx=dy =0.1

xl = lvec[1,0]
yl = lvec[2,1]
zl = lvec[3,2]
extent = (lvec[0,0],lvec[0,0]+xl,lvec[0,1],lvec[0,1]+yl)

tip_r2 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+xl,dx,lvec[0,1],lvec[0,1]+yl,dy,lvec[0,2],lvec[0,2]+zl,dz)

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

Voltages=[-0.338, -0.09406976,  0., 2.422684]
namez=['HOMO-3-2','HOMO-1','HOMO','LUMO++1']
#Voltages=np.arange(-1.0,+1.0+0.01,0.1) # this part is important for scans over slabs at different voltages
#namez = []
#for V in Voltages:
#    namez.append(str(round(V,1)))

# --- downloading the df data

df, lvec2, nDim2 = GU.load_scal_field( path_df+'df' ,data_format=data_format)

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

i=0;
for V in Voltages:
    current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
    current1 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
    current2 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px=0.5, py=0.5, pz = 0.0)
    # --- plotting part here, plots all calculated signals:
    print(" plotting ")
    for k in [3,11]:
        dff = np.array(df[k,:,:]).copy()
        curr0 = np.array(current0[k,:,:]).copy()
        curr1 = np.array(current1[k,:,:]).copy()
        curr2 = np.array(current2[k,:,:]).copy()
    
        name_file='didV-'+namez[i]+'_%03d.dat' %k
        name_plot_df='height:%03dA; df [Hz]' %k
        name_plot0=namez[i]+';height:%03dA; dIdV [G0] s-fixed-tip' %k
        name_plot1=namez[i]+';height:%03dA; dIdV [G0] s-tip' %k
        name_plot2=namez[i]+';height:%03dA; dIdV [G0] pxy-tip' %k
    
        # ploting part here:
        plt.figure( figsize=(1.5* xl , 1.5*yl/4 ) )
        plt.subplot(1,4,1)
        plt.imshow( dff, origin='image', extent=extent , cmap='gray')
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title(name_plot_df)
    
        # ploting part here:
        plt.subplot(1,4,2)
        plt.imshow( curr0, origin='image', extent=extent, cmap='gray' )
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title(name_plot0)
    
        # ploting part here:
        plt.subplot(1,4,3)
        plt.imshow( curr1, origin='image', extent=extent, cmap='gray' )
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title(name_plot1)
    
        # ploting part here:
        plt.subplot(1,4,4)
        plt.imshow( curr2, origin='image', extent=extent, cmap='gray' )
        plt.xlabel(r' Tip_x $\AA$')
        plt.ylabel(r' Tip_y $\AA$')
        plt.title(name_plot2)
        plt.savefig( 'didv_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
        #plt.show()
        plt.close()
    i += 1


# --- the end

print() 
print()
print("Done")
