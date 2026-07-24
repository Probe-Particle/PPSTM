#!/usr/bin/python
import math
import os
import numpy as np

import ppafm.io as io
import sys
sys.path.append('../../')
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS

import matplotlib
# matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

# --- specification of paths to the STM input files, PP positions and stored df results & format in which the data are stored - xsf or npy

path=''
path_pos='Q0.00K0.50/'
path_df = path_pos+'Amp0.40/'
data_format ="npy"

# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(0,0)
lvs = None
#lvs = np.array([[15., 0.,0.],[0.,15.,0],[0.,0.,15.]]
WorkFunction = 5.0 #more or less standart.
fermi=None    # the Fermi from phik ... .dat file; !!! All energies are relative to Fermi !!!! None -  means: -5.04612664712 eV
orbs= 'sp'    # 'sp' works now, 'spd' works for fireball as well
cut_min=-1.0    # HOMO -0.88 bellow the Fermi Level, other orbitals cut
cut_max=+1.0    # LUMO -0.88 above the Fermi Level
cut_at=-1    # All atoms of the molecule
eta = 0.01    # very low, to pronounce the single orbitals only
WF_decay=1.0    # for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
nV = 9        # for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]    # No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]    # Lowering of the hoppings

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
                        cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin = RS.read_AIMS_all(name = 'KS_eigenvectors_up.band_1.kpt_1.out', geom='geometry.in',fermi=fermi, orbs = 'sp', pbc=pbc,
#                    imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
#                    lower_atoms=lower_atoms, lower_coefs=lower_coefs)
#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#                    cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin  = RS.read_CP2K_all(name = 'crazy_mol', fermi=fermi, orbs = orbs, pbc=pbc,
#                    cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; tip_r1 - PP distored by the relaxation in the PPAFM code; tip_r2 - uniform grid:

tip_r1, lvec, nDim, atomic_info_or_head = io.load_vec_field( path_pos+'PPpos' ,data_format=data_format)

dz=0.1
dx=dy =0.1

xl = lvec[1,0]
yl = lvec[2,1]
zl = lvec[3,2]
extent = (lvec[0,0],lvec[0,0]+xl,lvec[0,1],lvec[0,1]+yl)

tip_r2 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+xl,dx,lvec[0,1],lvec[0,1]+yl,dy,lvec[0,2],lvec[0,2]+zl,dz)

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

Voltages=[-0.88,+0.88]
namez=['HOMO','LUMO']

#Voltages=np.arange(-1.0,+1.0+0.01,0.1) # this part is important for scans over slabs at different voltages
#namez = []
#for V in Voltages:
#    namez.append(str(round(V,1)))

# --- downloading the df data

df, lvec2, nDim2, atomic_info_or_head = io.load_scal_field( path_df+'df' ,data_format=data_format)

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

for wf in [WorkFunction]:
    for V, energy_orbital in zip(Voltages, namez):
        print ("Voltage:", V,"name:", energy_orbital)
        for eta in [eta]:
            scans = [
                df,
                PS.dIdV(V, wf, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0),
                PS.dIdV(V, wf, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0),
                PS.dIdV(V, wf, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px=1.0, py=1.0, pz = 0.0),
                PS.dIdV_tilt(V, wf, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, pz = 1.0, al =1.0),
                PS.dIdV_tilt(V, wf, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, dxyz = 1.0, al =1.0),
                PS.STM(V, nV, wf, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, px=0.5, py=0.5, WF_decay=WF_decay),
                PS.IETS_simple(V, wf, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px =0.5, py=0.5, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.02),
                PS.dIdV_tilt(V, wf, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, pxy = 1.0, al =1.0),
                PS.dIdV_tilt(V, wf, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, dz2 = 1.0, al =1.0),
            ]
            
            # --- plotting part here, plots all calculated signals:
            print(" plotting ")
            for tip_height in range(df.shape[0]):
                name_file=f'didv_{energy_orbital}_WF_{wf}_{eta}_{tip_height:03d}.png'

                subplot_titles = [
                    f'height:{tip_height:03d}A; df [Hz]',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] s-fixed-tip',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] s-tip',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] pxy-tip',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] pz-tip tilting',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] dxyz-tip tilting',
                    f'{energy_orbital};height:{tip_height:03d}A; STM [I] pxy-tip',
                    f'{energy_orbital};height:{tip_height:03d}A; Under develoment',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] pxy-tip tilting',
                    f'{energy_orbital};height:{tip_height:03d}A; dIdV [G0] dz2-tip tilting',
                ]

                plt.figure( figsize=(1.5* xl , 1.5*yl/2 ) )

                n_plot_rows = 2
                n_plot_cols = math.ceil(len(scans) / 2)

                for j, (scan_j, subtitle_j) in enumerate(zip(scans, subplot_titles)):
                    plt.subplot(n_plot_rows, n_plot_cols, j + 1)
                    plt.imshow(scan_j[tip_height], origin='lower', extent=extent, cmap='gray')
                    plt.xlabel(r' Tip_x $\AA$')
                    plt.ylabel(r' Tip_y $\AA$')
                    plt.title(subtitle_j, fontdict={'fontsize': 'medium'})
                
                plt.savefig(name_file, bbox_inches='tight' )
                plt.close()
                #plt.show()
                #
                # --- for saving WSxM format, if wanted
                #
                #tmp_curr=curr.flatten()
                #out_curr=np.zeros((len(tmp_curr),3))
                #out_curr[:,0]=tip_r[k,:,:,0].flatten()
                #out_curr[:,1]=tip_r[k,:,:,1].flatten()
                #out_curr[:,2]=tmp_curr.copy()
                #f=open(name_file,'w')
                #print >> f, "WSxM file copyright Nanotec Electronica"
                #print >> f, "WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al."
                #print >> f, "X[A]  Y[A]  dIdV[G0]"
                #print >> f, ""
                #np.savetxt(f, out_curr)
                #f.close()
                #
        
        #plt.show()


# --- the end

print() 
print()
print("Done")
