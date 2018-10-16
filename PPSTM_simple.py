#!/usr/bin/python
# !!!!!!!!!!!!!!!!!!!!!!!!!!!  TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################################################################################
#                                                                                                                        #
#                                  What follows are options of the PP-STM (dIdV) code                                    #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** Main informations ******
#
scan_type     = 'didv'               # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'fixed'              # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) #
V             = -2.0                 # !!!! V = Vmin for SCAN !!!! #
V_max         = +2.0                 # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1                 # voltage step , dV <= 0.1 V #
eta           =  0.1                 # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.0                 # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb       = 's'                  # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
sample_orbs   = 'spd'                # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
dft_code      = 'fireball'           # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
geometry_file = 'input.xyz'          # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW $
pbc           = (0,0)                # (0,0) = None = False -- only original geometry ; (0.5,0.5) -- 2x2 cell ; (1,1) -- 3x3 cell (around original) ; (2,2) -- 5x5 cell (around original) ... #
lvs           = None                 # None ; [[ax,ay,0],[bx,by,0]],[0,0,cz]] or [[ax,ay],[bx,by]] ; 'input.lvs' -- files with specified cell ; in FHI-AIMS & GPAW allready specified with geometry #
spin          = None                 # None=False ; for FHI-AIMS & CP2K: None, "both" , "up" or "down" (last 3 are for spin-polarizes or spin-unrestricted calculations #
















##########################################################################################################################
#                                                                                                                        #
#                                  DO NOT TOUCH LATER CODE (until you know what you are doing)                           #
#                                                                                                                        #
##########################################################################################################################

import os
import numpy as np
import pyProbeParticle.GridUtils as GU
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
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
fermi=None	# the Fermi from phik ... .dat file; !!! All energies are relative to Fermi !!!! None -  means: -5.04612664712 eV
orbs= 'sp'	# 'sp' works now, 'spd' works for fireball as well
cut_min=-1.0	# HOMO -0.88 bellow the Fermi Level, other orbitals cut
cut_max=+1.0	# LUMO -0.88 above the Fermi Level
cut_at=-1	# All atoms of the molecule
eta = 0.01	# very low, to pronounce the single orbitals only
WF_decay=1.0	# for STM only - how fast the exponential decay fall, with the applied bias ( if 1 - 1:1 correspondence with bias; if 0, it doesn't change)
nV = 9		# for STM only - number of STM integrational steps nV ~ V/eta
lower_atoms=[]	# No atoms has lowered hopping - be aware python numbering occurs here [0] - means lowering of the 1st atom
lower_coefs=[]	# Lowering of the hoppings

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = path+'phik_example_', geom=path+'crazy_mol.xyz', fermi=fermi, orbs = orbs, pbc=pbc,
					    cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin = RS.read_AIMS_all(name = 'KS_eigenvectors_up.band_1.kpt_1.out', geom='geometry.in',fermi=fermi, orbs = 'sp', pbc=pbc,
#					imaginary = False, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at,
#					lower_atoms=lower_atoms, lower_coefs=lower_coefs)
#eigEn, coefs, Ratin  = RS.read_GPAW_all(name = 'out_LCAO_LDA.gpw', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#eigEn, coefs, Ratin  = RS.read_CP2K_all(name = 'crazy_mol', fermi=fermi, orbs = orbs, pbc=pbc,
#					cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

# --- the grid on which the STM signal is calculated; tip_r1 - PP distored by the relaxation in the PPAFM code; tip_r2 - uniform grid:

tip_r1, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)

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

df, lvec2, nDim2 = GU.load_scal_field( path_df+'df' ,data_format=data_format)

# --- the Main Loop - for different WorkFunction (exponential z-decay of current), sample bias Voltages & eta - lorentzian FWHM

for WorkFunction in [WorkFunction]:
    i=0;
    for V in Voltages:
	for eta in [eta]:
	    current0 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r2, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
	    current1 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=1.0, px=0.0, py=0.0, pz = 0.0)
	    current2 = PS.dIdV( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px=1.0, py=1.0, pz = 0.0)
	    current3 = PS.dIdV_tilt( V, WorkFunction, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, pz = 1.0, al =1.0)
	    current4 = PS.dIdV_tilt( V, WorkFunction, eta, eigEn, tip_r1, tip_r2, Ratin, coefs, orbs=orbs, dxyz = 1.0, al =1.0)
	    current5 = PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, px=0.5, py=0.5, WF_decay=WF_decay)
	    # next procedure is under development
	    current6 = PS.IETS_simple( V, WorkFunction, eta, eigEn, tip_r1, Ratin, coefs, orbs=orbs, s=0.0, px =0.5, py=0.5, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.02)
	
	    # --- plotting part here, plots all calculated signals:
	    print " plotting "
	    for k in range(df.shape[0]):
		dff = np.array(df[k,:,:]).copy()
		curr0 = np.array(current0[k,:,:]).copy()
		curr1 = np.array(current1[k,:,:]).copy()
		curr2 = np.array(current2[k,:,:]).copy()
		curr3 = np.array(current3[k,:,:]).copy()
		curr4 = np.array(current4[k,:,:]).copy()
		curr5 = np.array(current5[k,:,:]).copy()
		curr6 = np.array(current6[k,:,:]).copy()
		
		name_file='didV-'+namez[i]+'_%03d.dat' %k
		name_plot_df='height:%03dA; df [Hz]' %k
		name_plot0=namez[i]+';height:%03dA; dIdV [G0] s-fixed-tip' %k
		name_plot1=namez[i]+';height:%03dA; dIdV [G0] s-tip' %k
		name_plot2=namez[i]+';height:%03dA; dIdV [G0] pxy-tip' %k
		name_plot3=namez[i]+';height:%03dA; dIdV [G0] pz-tip tilting' %k
		name_plot4=namez[i]+';height:%03dA; dIdV [G0] dxyz-tip tilting' %k
		name_plot5=namez[i]+';height:%03dA; STM [I] pxy-tip' %k
		name_plot6=namez[i]+';height:%03dA; Under develoment' %k
		
		# ploting part here:
		plt.figure( figsize=(1.5* xl , 1.5*yl/2 ) )
		plt.subplot(2,4,1)
		plt.imshow( dff, origin='image', extent=extent , cmap='gray')
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot_df)
		
		# ploting part here:
		plt.subplot(2,4,2)
		plt.imshow( curr0, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot0)

		# ploting part here:
		plt.subplot(2,4,3)
		plt.imshow( curr1, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot1)
		
		# ploting part here:
		plt.subplot(2,4,4)
		plt.imshow( curr2, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot2)
		
		plt.subplot(2,4,5)
		plt.imshow( curr3, origin='image', extent=extent , cmap='gray')
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot3)
		
		# ploting part here:
		plt.subplot(2,4,6)
		plt.imshow( curr4, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot4)
		
		# ploting part here:
		plt.subplot(2,4,7)
		plt.imshow( curr5, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot5)
		
		# ploting part here:
		plt.subplot(2,4,8)
		plt.imshow( curr6, origin='image', extent=extent, cmap='gray' )
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title(name_plot6)
		
		plt.savefig( 'didv_'+namez[i]+"_WF_"+str(WorkFunction)+"_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
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
	i = i+1
	

# --- the end

print 
print
print "Done"
