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
# ***** System information: path (absolute or relative) to your PPSTM code *****
#
ppstm_path = '$HOME/Program_Files/PPSTM/'
#
# ***** Main informations ******
#
scan_type     = 'didv'       # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'fixed'      # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
V             = -2.0         # !!!! V = Vmin for SCAN !!!! #
V_max         = +2.0         # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1         # voltage step , dV <= 0.1 V #
eta           =  0.1         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.0         # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb       = 's'          # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
sample_orbs   = 'spd'        # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
dft_code      = 'fireball'   # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
geometry_file = 'input.xyz'  # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW #
pbc           = (0,0)        # (0,0) = None = False -- only original geometry ; (0.5,0.5) -- 2x2 cell ; (1,1) -- 3x3 cell (around original) ; (2,2) -- 5x5 cell (around original) ... #
lvs           = None         # None ; [[ax,ay,0],[bx,by,0]],[0,0,cz]] or [[ax,ay],[bx,by]] ; 'input.lvs' -- files with specified cell ; in FHI-AIMS & GPAW allready specified with geometry #
spin          = None         # None=False ; for FHI-AIMS & CP2K: None -- spin-unpolarized/spin-restricted calc. ;  'both' , 'up'='alpha' or 'down" (last 3 are for spin-polarizes or spin-unrestricted calculations #
cp2k_name     = 'CuPc'       # Name used in CP2K calculations or GPAW calc/ #
#
# ***** Informations for x,y,z tip_type = 'fixed' ******
#
x = [  0.0, 20.0, 0.25 ]     # [xmin, xmax, dx] #
y = [  0.0, 15.0, 0.25 ]     # [ymin, ymax, dy] #
z = [ 10.0, 12.0, 0.1  ]     # !!!! z-starts from zero - normally zmin >= 3-4 Ang above heighest atoms !!!! [zmin, zmax, dz] ; for single height scan use : [z, z, 1.0] #
#
# ***** Informations for PP positions, tip_type = 'relaxed' ******
#
Q = 0.00                     # charge (PP-AFM) ; Ocharge PP-AFM complex_tip autumn 2018) ; [e] (monopole), [e*A] (dipole), [e*A^2] (quadrupole) #
K = 0.24                     # x stiffness (PP-AFM master autumn 2018); klat (PP-AFM dev/OpenCl autumn 2018); Oklat (PP-AFM complex_tip autumn 2018) ; [N/m] #
data_format = 'xsf'          # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which PPpos are stored from PP-AFM run #
#
# *****Output options ******
#
png  = True                  # True / False -- plot "png" images (2D constant height) #
WSxM = False                 # True / False -- write ".xyz" WSxM files (2D constant height) #
XSF  = False                 # True / False -- write ".xsf" files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
NPY  = False                 # True / False -- write ".npy" numpy binary files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
#
# ***** Advanced options ******
#
cut_atoms   = None           # None = -1 -- All atoms of the sample contributes to tunelling ; 1 -- only 1st atom of the sample contributes to the tunelling ; 57 -- first 57 atoms of the sample contributes to the tunelling ; ... #
lower_atoms = []             # [] = None -- No atoms has lowered hopping ; be aware python numbering occurs here: [0] - means lowering of the 1st atom; [0,1,2,3] -- lowering of 1st 4 atoms ... #
lower_coefs = []             # [] = None -- No lowering of the hoppings  ; [0.5] -- lowering of the 1st atom hopping to 0.5                           ; [0.5,0.5,0.5,0.5] -- lowering of 1st 4 atoms to 0.5 ... #
#
# ***** More advanced options ******
#
WorkFunction =  5.0          # 5.0 eV is standart #
fermi        = None          # None=0.0 -- no change to the Fermi Level ; -0.1 -- shifts the Fermi Level by 0.1 eV lower ... #
cut_min      = -2.5          # cut out all orbitals lower than  -2.5 eV bellow Rermi (should be: cut_min <= Vmin-2*eta) . taken to the Fermi Level #
cut_max      = +2.5          # cut out all orbitals higher than -2.5 eV above  Fermi (should be: cut_max >= Vmax+2*eta) . taken to the Fermi Level #
files_path   = ''            # where are files fron DFT code ; rather do not use this #
#
#
##########################################################################################################################
#                                                                                                                        #
#                                 DO NOT TOUCH LATER CODE (unless you know what you are doing)                           #
#                                                                                                                        #
##########################################################################################################################

print "Importing libraries"

import os
import sys
sys.path.append(ppstm_path) 

import numpy as np
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
if (XSF or NPY or (tip_tipe == 'relaxed') or (tip_tipe == 'r'_):
    print "For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory "
    import pyProbeParticle.GridUtils as GU

print "Libraries imported"

# --- the grid on which the STM signal is calculated --- #

if (tip_type =='relaxed') or (tip_tipe == 'r'_):
    print "Importing positions of PP from the PP-AFM calculations. Path for the data:"
    path_pos="Q%1.2fK%1.2f" %(Q,K)
    print path_pos
    tip_r, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)
    extent = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[1,1])
    print "DEBUG: extent", extent
    print "PP postions imported"
else:
    print "Priparing the scan grid for fixed scan"
    extent = (x[0],x[1],y[0],y[1])
    tip_r  = RS.mkSpaceGrid(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2])
    lvec   = np.array([[x[0],y[0],z[0]],[x[1]-x[0],0.,0.],[0.,y[1]-y[0],0.],[0.,0.,z[1]-z[0]]])
    print "DEBUG: extent", extent
    print "DEBUG: lvec", lvec
    print "scan grids prepared"

# --- downloading and examples of downloading of the eigen-energies, the LCAO coefficients and geometry (this time for spin-unpolarized calculations):

print "Reading electronic & geometry structure files"

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL')):
    if isinstance(lvs, (list, tuple, np.ndarray)):
	
    else
	cell = np.loadtxt(lvs)
    eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = files_path + 'phik_0001_', geom=path+geometry_file, lvs = cell, fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

elif ((dft_code == 'gpaw') or(dft_code == 'GPAW')):
    eigEn, coefs, Ratin = RS.read_GPAW_all(    name = cp2k_name + '.gpw', fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

elif ((dft_code == 'aims') or(dft_code == 'AIMS') or (dft_code == 'FHI-AIMS')):
    if (spin == None):
	name = 'KS_eigenvectors.band_1.kpt_1.out'
    elif ((spin == 'up')or(spin == 'alpha')or(spin='=both'):
	name = 'KS_eigenvectors_up.band_1.kpt_1.out'
    else:
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
    eigEn, coefs, Ratin = RS.read_AIMS_all(name = name , geom=geometry_file, fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    if (spin == 'both'):
	eigEn1 = eigEn.copy(); coefs1 = coefs.copy(); del eigEn, coefs ;
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
	eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = name , geom=geometry_file, fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)

elif ((dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if (spin == None):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = cp2k_name , fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    elif ((spin == 'up')or(spin == 'alpha')):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = cp2k_name , fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
    elif (spin == 'both'):
	eigEn1, coefs1, Ratin  = RS.read_CP2K_all(name = cp2k_name , fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
	eigEn2, coefs2, Ratin  = RS.read_CP2K_all(name = cp2k_name , fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)
    else:
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = cp2k_name , fermi=fermi, orbs = orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');

print "DEBUG: eigEn.shape ", eigEn.shape
print "DEBUG: coefs.shape ", coefs.shape
print "DEBUG: Ratin.shape ", Ratin.shape

# --- specification on which voltages the STM (dI/dV ...) calculations are performed - two methods - direct specification or sequence of voltages

# Not necessary

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
