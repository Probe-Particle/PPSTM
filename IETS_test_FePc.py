#!/usr/bin/python

import os
import numpy as np

import pyPPSTM.GridUtils as GU
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt

# --- specification of paths to the STM input files, PP positions and stored df results & format in which the data are stored - xsf or npy

path=''
path_pos='Q0.00K0.24/'
path_df = path_pos+'Amp1.00/'
data_format ="npy"
wsxm      = False
save_npy  = False
Plot_iets = True
Plot_dI2  = True

# --- specification of PBC, cell, and All the important stuff concerning electrons tunneling:

pbc=(0,0)		# just try
lvs = None      # automatically taken from geometry.in
WorkFunction = 5.0 #more or less standart.
fermi=None	# the Fermi from AIMS automatically 0 !!! All energies are relative to Fermi !!!! None -  means: -5.04612664712 eV
orbs= 'sp'	# 'sp' works now, 'spd' works for fireball as well
cut_min=-1.0	# 
cut_max=+1.0	# 
cut_at=57	# All atoms of the molecule (MetalPc - 57 atoms, H2Pc - 58, NoPc - 56)
#eta = 0.1	# very low, to pronounce the single orbitals only
# -- next two not needed for IETS
#WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
#nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings
M = 16		# Effective mass of CO frustrated translation - only O - 16

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

#eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
#					    cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
eigEn1, coefs1, Ratin = RS.read_AIMS_all(name = "eigen_up.out", geom='geom-cube.in',fermi=fermi, orbs = orbs, pbc=pbc,
                    imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
                    lower_atoms=lower_atoms, lower_coefs=lower_coefs)
eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = "eigen_dn.out", geom='geom-cube.in',fermi=fermi, orbs = orbs, pbc=pbc,
                    imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
                    lower_atoms=lower_atoms, lower_coefs=lower_coefs)

eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
#print eigEn
coefs = np.concatenate((coefs1, coefs2), axis=0)

#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; tip_r1 - PP distored by the relaxation in the PPAFM code; tip_r2 - uniform grid:

tip_r1, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)
eigenEner, lvec, nDim = GU.load_vec_field( path_pos+'eigvalKs' ,data_format=data_format)
eigenVec1, lvec, nDim = GU.load_vec_field( path_pos+'eigvecK1' ,data_format=data_format)
eigenVec2, lvec, nDim = GU.load_vec_field( path_pos+'eigvecK2' ,data_format=data_format)
eigenVec3, lvec, nDim = GU.load_vec_field( path_pos+'eigvecK3' ,data_format=data_format)

dz=0.1
dx=dy =0.1

xl = lvec[1,0]
yl = lvec[2,1]
zl = lvec[3,2]
extent = (lvec[0,0],lvec[0,0]+xl,lvec[0,1],lvec[0,1]+yl)

tip_r2 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+xl,dx,lvec[0,1],lvec[0,1]+yl,dy,lvec[0,2],lvec[0,2]+zl,dz)
ddown=5
upp=15

#ddown=8
#upp=8

ki = ddown

tip_r1 = tip_r1[ddown:upp+1]
tip_r2 = tip_r2[ddown:upp+1]

eigenEner = eigenEner[ddown:upp+1]
eigenVec1 = eigenVec1[ddown:upp+1]
eigenVec2 = eigenVec2[ddown:upp+1]
eigenVec3 = eigenVec3[ddown:upp+1]

# --- downloading the df data

df, lvec2, nDim2 = GU.load_scal_field( path_df+'df' ,data_format=data_format)
df = df [ddown:upp+1]
for ii in range(ddown,upp+1):
    tmp=np.loadtxt(path_pos+'IETS_%03d.xyz' % (ii),skiprows=4) 
    #tmp2 = np.reshape()
    iets_afm = np.array([tmp]) if ii == ddown else np.append(iets_afm, np.array([tmp]),axis=0)

iets_afm = iets_afm.reshape((df.shape[0],df.shape[1],df.shape[2],3))

print("Debug: df.shape", df.shape)
print("Debug: iets_afm.shape", iets_afm.shape)

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

#Voltages=[0.0]
#namez=['0.0']

Voltages=np.arange(-0.1,0.1+0.01,0.1) # this part is important for scans over slabs at different voltages
namez = []
namez_der = []
for V in Voltages:
    namez.append(str(round(V,1)))
    namez_der.append(str(round(V-0.05,2)))



# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

lvec1 = lvec
lvec1 [3,2] = (upp+1-ddown)*dz
lvec1 [0,2] += ddown*dz

sh = tip_r1.shape

for WorkFunction in [WorkFunction]:
    i=0;
    for V in Voltages:
#	    if (-0.0001<V<0.0001):
        print("V:",V)
        curs = np.zeros((sh[0],sh[1],sh[2],1));
        curp = curs.copy();
        ietss = curs.copy(); ietsp = curs.copy();
        IETSs = curs.copy(); IETSp = curs.copy();
        j = 0;
        for eta in [1.0]:
            print("eta: ",eta)
            current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
            current1 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px=0.5, py=0.5, pz = 0.0)
            # next procedure is under development
            denomin1, current3, current4 = PS.IETS_complex( V, WorkFunction, eta, eigEn, tip_r1, eigenEner, eigenVec1, eigenVec2, eigenVec3, Ratin, coefs, orbs=orbs, s=1.0, px =0.0, py=0.0, pz=0.0, Amp=0.05, M=M)
            denomin1, current5, current6 = PS.IETS_complex( V, WorkFunction, eta, eigEn, tip_r1, eigenEner, eigenVec1, eigenVec2, eigenVec3, Ratin, coefs, orbs=orbs, s=0.0, px =0.5, py=0.5, pz=0.0, Amp=0.05, M=M)
            curs[:,:,:,j]= current0;
            curp[:,:,:,j]= current1;
            ietss[:,:,:,j]= current3;
            ietsp[:,:,:,j]= current5;
            IETSs[:,:,:,j]= current4;
            IETSp[:,:,:,j]= current6;
            j+=1;

        # --- saving part here
        if save_npy:
            print("saving, V:",V)
            j=0;
            for eta in [1.0]:
                GU.save_scal_field( path_pos+'curs'+namez[i], curs[:,:,:,j], lvec1, data_format=data_format )
                GU.save_scal_field( path_pos+'curp'+namez[i], curp[:,:,:,j], lvec1, data_format=data_format )
                GU.save_scal_field( path_pos+'ietss'+namez[i], ietss[:,:,:,j], lvec1, data_format=data_format )
                GU.save_scal_field( path_pos+'ietsp'+namez[i], ietsp[:,:,:,j], lvec1, data_format=data_format )
                GU.save_scal_field( path_pos+'IETSs'+namez[i], IETSs[:,:,:,j], lvec1, data_format=data_format )
                GU.save_scal_field( path_pos+'IETSp'+namez[i], IETSp[:,:,:,j], lvec1, data_format=data_format )
            j+=1

        # --- plotting part here, plots all calculated signals:
        if wsxm:
            print(" plotting wsxm")
            for eta in [1.0]:
            #import pyProbeParticle.GridUtils as GU
                print(" printing current into WSxM files :")
                GU.saveWSxM_3D(path_pos+"current_eta_"+str(eta)+"_"+str(ki)+"+" , 0.15*curs[:,:,:,0]+curp[:,:,:,0] , extent , slices=None)

                print(" printing IETS-stm into WSxM files :")
                GU.saveWSxM_3D(path_pos+"IETS-stm_part_eta_"+str(eta)+"_"+str(ki)+"+" , 0.15*ietss[:,:,:,0]+ietsp[:,:,:,0] , extent , slices=None)

                print(" printing IETS_amplitude into WSxM files :")
                GU.saveWSxM_3D(path_pos+"IETS_amplitude_eta_"+str(eta)+"_"+str(ki)+"+" , 0.15*IETSs[:,:,:,0]+IETSp[:,:,:,0] , extent , slices=None)

            
        if Plot_iets:
            print(" plotting IETS images")
            for k in range(len(current0)):
                name_plot0 =namez[i]+';height:%03d$\AA$; dIdV [G0] s-tip'              %(k+ki)
                name_plot1 =namez[i]+';height:%03d$\AA$; dIdV [G0] pxy-tip'            %(k+ki)
                name_plot2 =namez[i]+';height:%03d$\AA$; dIdV [G0] s-pxy-tip'          %(k+ki)

                name_plot3 =namez[i]+';height:%03d$\AA$; iets(stm) [G0/$\AA$^2] s-tip'              %(k+ki)
                name_plot4 =namez[i]+';height:%03d$\AA$; iets(stm) [G0/$\AA$^2] pxy-tip'            %(k+ki)
                name_plot5 =namez[i]+';height:%03d$\AA$; iets(stm) [G0/$\AA$^2] s-pxy-tip'          %(k+ki)

                name_plot6 =namez[i]+';height:%03d$\AA$; IETS [?] s-tip'              %(k+ki)
                name_plot7 =namez[i]+';height:%03d$\AA$; IETS [?] pxy-tip'            %(k+ki)
                name_plot8 =namez[i]+';height:%03d$\AA$; IETS [?] s-pxy-tip'          %(k+ki)

                name_plot9 =namez[i]+';height:%03d$\AA$; df [hz]'                     %(k+ki)
                name_plot10=namez[i]+';height:%03d$\AA$; iets(AFM only) [?]'          %(k+ki)
                name_plot11=namez[i]+';height:%03d$\AA$; denominators   [s]'          %(k+ki)

                # ploting part here:
                plt.figure( figsize=(3./3* xl , 3./3*yl ) )

                plt.subplot(4,3,1)
                plt.imshow(  curs[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.ylabel(r' Tip_y $\AA$; eta =1.00 eV')
                plt.title(name_plot0)

                plt.subplot(4,3,2)
                plt.imshow(  curp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot1)

                plt.subplot(4,3,3)
                plt.imshow(  0.15*curs[k,:,:,0]+curp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot2)

                plt.subplot(4,3,4)
                plt.imshow(  ietss[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.ylabel(r' Tip_y $\AA$; eta =1.00 eV')
                plt.title(name_plot3)

                plt.subplot(4,3,5)
                plt.imshow(  ietsp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot4)

                plt.subplot(4,3,6)
                plt.imshow(  0.15*ietss[k,:,:,0]+ietsp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot5)

                plt.subplot(4,3,7)
                plt.imshow(  IETSs[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.ylabel(r' Tip_y $\AA$; eta =1.00 eV')
                plt.title(name_plot6)

                plt.subplot(4,3,8)
                plt.imshow(  IETSp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot7)

                plt.subplot(4,3,9)
                plt.imshow(  0.15*IETSs[k,:,:,0]+IETSp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot8)

                plt.subplot(4,3,10)
                plt.imshow(  df[k,:,:], origin='image', extent=extent, cmap='gray' )
                plt.ylabel(r' Tip_y $\AA$;PP-AFM code')
                plt.title(name_plot9)
                plt.xlabel(r' Tip_x $\AA$')

                plt.subplot(4,3,11)
                plt.imshow(  iets_afm[k,:,:,2], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot10)
                plt.xlabel(r' Tip_x $\AA$')

                plt.subplot(4,3,12)
                plt.imshow(  denomin1[k,:,:], origin='image', extent=extent, cmap='gray' )
                plt.title(name_plot11)
                plt.xlabel(r' Tip_x $\AA$')



                plt.savefig('All_'+namez[i]+"_fermi_"+str(fermi)+'_%03d.png' %(k+ki) , bbox_inches='tight' )
                plt.close()
            
        if ((i>0)and Plot_dI2):
            print(" plotting d^2I/dV^2 ")
            for k in range(current6.shape[0]):

                name_plot0 =namez[i]+';height:%03d$\AA$; d^2I/dV^2 [G0] s-tip'              %(k+ki)
                name_plot1 =namez[i]+';height:%03d$\AA$; d^2I/dV^2 [G0] pxy-tip'            %(k+ki)
                name_plot2 =namez[i]+';height:%03d$\AA$; d^2I/dV^2 [G0] s-pxy-tip'          %(k+ki)
            
                # ploting part here:
                plt.figure( figsize=(xl , 4./3*yl ) )
            
                plt.subplot(1,3,1)
                plt.imshow(  curs[k,:,:,0] - o_curs[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.ylabel(r' Tip_y $\AA$; eta =1.0 eV')
                plt.xlabel(r' Tip_x $\AA$')

                plt.subplot(1,3,2)
                plt.imshow(  curp[k,:,:,0] - o_curp[k,:,:,0], origin='image', extent=extent, cmap='gray' )
                plt.xlabel(r' Tip_x $\AA$')

                plt.subplot(1,3,3)
                plt.imshow(  0.15*curs[k,:,:,0]+curp[k,:,:,0] - (0.15*o_curs[k,:,:,0]+o_curp[k,:,:,0]), origin='image', extent=extent, cmap='gray' )
                plt.xlabel(r' Tip_x $\AA$')

                plt.savefig( 'd2IdV2_'+namez_der[i]+"_fermi_"+str(fermi)+'_%03d.png' %(k+ki) , bbox_inches='tight' )
                plt.close()
        
        o_curs = curs.copy();
        o_curp = curp.copy();
        
        i = i+1

# --- the end

print() 
print()
print("Done")
