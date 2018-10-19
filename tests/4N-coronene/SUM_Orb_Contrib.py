#!/usr/bin/python
#
##########################################################################################################################
#                                                                                                                        #
#                         What follows are options for summing tip-orbitals' contributions                               #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** System information: path (absolute or relative) to your PPSTM code *****
#
ppstm_path = './PPSTM/'
#
# ***** Main informations ******
#
scan_type     = 'didv'       # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'relaxed'    # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
V             = -0.5         # !!!! V = Vmin for SCAN !!!! #
V_max         = +0.5         # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1         # voltage step , dV <= 0.1 V #
eta           =  0.1         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.0         # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb1      = 's'          # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
tip_orb2      = 'pxy'        # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
sample_orbs   = 'sp'         # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
data_format   = 'npy'        # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which PPpos are stored from PP-AFM run #
#
# *****Output options ******
#
tip_orb1_amount = 0.10       #
tip_orb2_amount = 0.20       # 
PNG  = True                  # True / False -- plot "png" images (2D constant height) #
WSxM = True                  # True / False -- write ".xyz" WSxM files (2D constant height) #
XSF  = True                  # True / False -- write ".xsf" files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
NPY  = True                  # True / False -- write ".npy" numpy binary files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
#
# ***** More advanced options ******
#
WorkFunction =  5.0          # 5.0 eV is standart #
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
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
print "For XSF or NPY outputs and inputs you have to have installed PPAFM in your PPSTM directory "
import pyProbeParticle.GridUtils as GU

print "Libraries imported"

# --- Initial check --- #

assert( PNG or WSxM or XSF or NPY ), "No output set to be True; I'm not going to do anything if there is no output. I'm too lazy like a Gartfield. "
didv_b   = False
STM_b    = False
V_scan_b = False

if ( (scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
    didv_b = True
elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
    STM_b = True
else:
    didv_b   = True
    STM_b    = True
    V_scan_b = True

# --- importing STM (dIdV) calculated signal on the grid --- #

if V_scan_b:
    ran = range(V, V_max, dV);
else:
    ran = [V]
if didv_b :
    print "Importing dIdV data"
    print name1
    tip_r, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)
    extent = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[2,1])
    print "DEBUG: extent", extent
    print "PP postions imported"
if STM_b:
    print "Priparing the scan grid for fixed scan"
    extent = (x[0],x[1],y[0],y[1])
    tip_r  = RS.mkSpaceGrid(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2])
    lvec   = np.array([[x[0],y[0],z[0]],[x[1]-x[0],0.,0.],[0.,y[1]-y[0],0.],[0.,0.,z[1]-z[0]]])
    print "DEBUG: extent", extent
    print "DEBUG: lvec", lvec
    print "scan grids prepared"

# --- reading of the eigen-energies, the LCAO coefficients and geometry --- #

print "Reading electronic & geometry structure files"

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL')):
    if isinstance(lvs, (list, tuple, np.ndarray)):
	cell = np.array([[lvs[0][0],lvs[0][1],0.0],[lvs[1][0],lvs[1][1],0.0],[0.0,0.0,99.9]]) if (len(lvs) == 2) else lvs
    elif isinstance(lvs, (str)):
	cell = np.loadtxt(lvs)
    elif ((pbc == (0,0)) or (pbc == (0.,0.))):
	cell = np.array([[0,0,0],[0,0,0],[0,0,0]]);
    else:
	print "PBC required, but lattice vector not specified. What can I do with that? I rather go to eat something."; exit()
    print "DEBUG: cell.shape", cell.shape
    eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = files_path + 'phik_0001_', geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

elif ((dft_code == 'gpaw') or(dft_code == 'GPAW')):
    eigEn, coefs, Ratin = RS.read_GPAW_all(    name = files_path + cp2k_name + '.gpw', fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

elif ((dft_code == 'aims') or(dft_code == 'AIMS') or (dft_code == 'FHI-AIMS')):
    if (spin == None):
	name = 'KS_eigenvectors.band_1.kpt_1.out'
    elif ((spin == 'up')or(spin == 'alpha')or(spin == 'both')):
	name = 'KS_eigenvectors_up.band_1.kpt_1.out'
    else:
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
    eigEn, coefs, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    if (spin == 'both'):
	eigEn1 = eigEn.copy(); coefs1 = coefs.copy(); del eigEn, coefs ;
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
	eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)

elif ((dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if (spin == None):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    elif ((spin == 'up')or(spin == 'alpha')):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
    elif (spin == 'both'):
	eigEn1, coefs1, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
	eigEn2, coefs2, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)
    else:
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');

print "DEBUG: eigEn.shape ", eigEn.shape
print "DEBUG: coefs.shape ", coefs.shape
print "DEBUG: Ratin.shape ", Ratin.shape

print "energies prepared, coeffecients read"

# --- the Main calculations --- #
# 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #

didv_b = False
STM_b  = False

if ( (scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
    didv    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
    didv_b = True
    print "DEBUG: didv.shape ", didv.shape
elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
    nV = abs(V/dV)+1
    print "DEBUG: V, nV:", V, nV
    current = np.array([   PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay) ])
    STM_b = True
    print "DEBUG: current.shape ", current.shape
else:
    current, didv = PS.MSTM( V, V_max, dV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay)
    didv_b = True
    STM_b  = True
    print "DEBUG: didv.shape ", didv.shape
    print "DEBUG: current.shape ", current.shape

# --- plotting part here, plots all calculated signals --- #

Voltages=np.arange(V,V_max+0.001,dV) # this part is important for scans over slabs at different voltages
namez = []
for V in Voltages:
    namez.append(str(round(V,2)))

NoV = len(didv) if didv_b else len(current)
NoH = len(didv[0]) if didv_b else len(current[0])

print "DEBUG: Voltages", Voltages
print "DEBUG: namez", namez
print "DEBUG: NoV", NoV
print "DEBUG: NoH", NoH

if PNG :
    print "We go to plotting "
    for vv in range(NoV):
	for k in range(NoH):
	    #print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
	    name_plot=namez[vv]+';height:'+str(k)+';tip:'+tip_type+';'+tip_orb
	    if didv_b :
		# ploting part here:
		plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
		plt.imshow(didv[vv,k,:,:], origin='image', extent=extent , cmap='gray')
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title("dIdV:"+name_plot)
		plt.savefig( 'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction+vv*WF_decay)+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
		plt.close()
	    if STM_b :
		# ploting part here:
		plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
		plt.imshow(current[vv,k,:,:], origin='image', extent=extent , cmap='gray')
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title("STM:"+name_plot)
		plt.savefig( 'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
		plt.close()
    print "Everything plotted"
if WSxM :
    print "writing WSxM files"
    for vv in range(NoV):
	for k in range(NoH):
	    if didv_b :
		name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction+vv*WF_decay)+"_eta_"+str(eta)+'_%03d.xyz' %k 
		tmp_curr=didv[vv,k,:,:].flatten()
		out_curr=np.zeros((len(tmp_curr),3))
		out_curr[:,0]=tip_r[k,:,:,0].flatten()
		out_curr[:,1]=tip_r[k,:,:,1].flatten()
		out_curr[:,2]=tmp_curr.copy()
		f=open(name_file,'w')
		print >> f, "WSxM file copyright Nanotec Electronica"
		print >> f, "WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al."
		print >> f, "X[A]  Y[A]  Z[A]"
		print >> f, ""
		np.savetxt(f, out_curr)
		f.close()
		#
	    if STM_b :
		name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d.xyz' %k 
		tmp_curr=current[vv,k,:,:].flatten()
		out_curr=np.zeros((len(tmp_curr),3))
		out_curr[:,0]=tip_r[k,:,:,0].flatten()
		out_curr[:,1]=tip_r[k,:,:,1].flatten()
		out_curr[:,2]=tmp_curr.copy()
		f=open(name_file,'w')
		print >> f, "WSxM file copyright Nanotec Electronica"
		print >> f, "WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al."
		print >> f, "X[A]  Y[A]  Z[A]"
		print >> f, ""
		np.savetxt(f, out_curr)
		f.close()
		#
    print "WSxM files written"

if XSF :
    print "writing XSF files"
    for vv in range(NoV):
	if didv_b :
	    name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction+vv*WF_decay)+"_eta_"+str(eta)+'.xsf'
	    GU.saveXSF(name_file, didv[vv], lvec)#, head=XSF_HEAD_DEFAULT )
	if STM_b :
	    name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'.xsf'
	    GU.saveXSF(name_file, current[vv], lvec)#, head=XSF_HEAD_DEFAULT )
    print "XSF files written"

if NPY :
    print "writing npy binary files"
    for vv in range(NoV):
	if didv_b :
	    name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction+vv*WF_decay)+"_eta_"+str(eta)
	    GU.saveNpy(name_file, didv[vv], lvec)#, head=XSF_HEAD_DEFAULT )
	if STM_b :
	    name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)
	    GU.saveNpy(name_file, current[vv], lvec)#, head=XSF_HEAD_DEFAULT )
    print "npy files written"

# --- the end --- #

print 
print
print "Done"
print
