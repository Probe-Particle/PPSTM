#!/usr/bin/python
#
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
ppstm_path = '/storage/praha4-fzu/home/shayaned/bin/edit/PPSTM-OMP/'
#
# ***** Main informations ******
mode          = 'all'      # 'real' , 'imag' or (coming soon) 'all'. LCAO coefs.
kInputName = 'phik_0002_'   # This is for fireball input files only. Select your k-point file. No need to put in the orbital names.
kpointFile = 'temp.kpts'
scan_type     = 'dIdV'     # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'f'    # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
Voltages      = [+0.56413775]
#V             = -0.249012     # !!!! V = Vmin for SCAN !!!! # -0.086174479, -0.155764479
V_max         = +2.0         # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1         # voltage step , dV <= 0.1 V #
eta           =  0.01         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.5         # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb       = 's'          # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
sample_orbs   = 'sp'         # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
dft_code      = 'FIREBALL'   # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
geometry_file = 'input.xyz' # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW #
pbc           = (2,0)        # (0,0) = None = False -- only original geometry ; (0.5,0.5) -- 2x2 cell ; (1,1) -- 3x3 cell (around original) ; (2,2) -- 5x5 cell (around original) ... #
lvs           =  'annulene.lvs'  #[[6.89188341,0.0],[0.0,30]]        # None ; [[ax,ay,0],[bx,by,0]],[0,0,cz]] or [[ax,ay],[bx,by]] ; 'input.lvs' -- files with specified cell ; in FHI-AIMS & GPAW allready specified with geometry #
spin          =  None        # None=False ; for FHI-AIMS & CP2K: None -- spin-unpolarized/spin-restricted calc. ;  'both' , 'up'='alpha' or 'down" (last 3 are for spin-polarizes or spin-unrestricted calculations) #
#cp2k_name     = 'crazy_mol'  # Name used in CP2K calculations or GPAW calc #
#
# ***** Informations for x,y,z tip_type = 'fixed' ******
#
x = [ -15.0, 20.0, 0.1 ]     # [xmin, xmax, dx] #
y = [ -2.0, 14.0, 0.1 ]     # [ymin, ymax, dy] #
z = [  5.0, 5.0, 0.01  ]     # !!!! z-starts from zero - normally zmin >= 3-4 Ang above heighest atoms !!!! [zmin, zmax, dz] ; for single height scan use : [z, z, 1.0] #
#
# ***** Informations for PP positions, tip_type = 'relaxed' ******
#
#Q = 0.00                     # charge (PP-AFM) ; Ocharge PP-AFM complex_tip autumn 2018) ; [e] (monopole), [e*A] (dipole), [e*A^2] (quadrupole) #
#K = 0.50                     # x stiffness (PP-AFM master autumn 2018); klat (PP-AFM dev/OpenCl autumn 2018); Oklat (PP-AFM complex_tip autumn 2018) ; [N/m] #
#data_format = 'npy'          # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which PPpos are stored from PP-AFM run #
#
# *****Output options ******
#
PNG  = True                  # True / False -- plot "png" images (2D constant height) #
WSxM = False                 # True / False -- write ".xyz" WSxM files (2D constant height) #
XSF  = False                 # True / False -- write ".xsf" files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
NPY  = False                 # True / False -- write ".npy" numpy binary files with 3D stucks of data . For this option you have to have "installed" PPAFM in your PPSTM directory #
plot_atoms = True            # True / False -- plot geometry (position of atoms into the PNG images and into the XSF files). You have to have your geometry, which you want to plot in input_plot.xyz. This doesn't change the name of the output files #
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
cut_min      = -2.0          # cut out all orbitals lower than  -2.5 eV below Fermi (should be: cut_min <= Vmin-2*eta) . taken to the Fermi Level #
cut_max      = +2.0          # cut out all orbitals higher than -2.5 eV above  Fermi (should be: cut_max >= Vmax+2*eta) . taken to the Fermi Level #
files_path   = '/storage/praha4-fzu/home/shayaned/fireball/pentalene/'            # where are files fron DFT code ; rather do not use this #
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
import pyPPSTM.basUtils as Bu

import numpy as np
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
from pylab import genfromtxt
import math
if (XSF or NPY or (tip_type == 'relaxed') or (tip_type == 'r' )):
    print "For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory "
    import pyProbeParticle.GridUtils as GU
if (plot_atoms):
    import pyPPSTM.basUtils as Bu
    import pyPPSTM.elements as elements


print "Libraries imported"

# --- Initial check --- #

assert( PNG or WSxM or XSF or NPY ), "No output set to be True; I'm not going to do anything if there is no output. I'm too lazy like a Gartfield. "

# --- specification of tip orbitals --- # 
# 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #

if (tip_orb == 's'):
    tc = [1.,0.,0.,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'pxy'):
    tc = [0.,0.5,0.5,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'spxy'):
    tc = [0.5,0.25,0.25,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == '5spxy'):
    tc = [0.05,0.475,0.475,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == '10spxy'):
    tc = [0.10,0.45,0.45,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'CO'):
    tc = [0.15,0.5,0.5,0.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'pz'):
    tc = [0.,0.,0.,1.,0.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'dz2'):
    tc = [0.,0.,0.,0.,1.,0.,0.] # [s, px, py, pz, dz2, dxz, dyz ] 
elif (tip_orb == 'dxzyz'):
    tc = [0.,0.,0.,0.,0.,0.5,0.5] # [s, px, py, pz, dz2, dxz, dyz ] 
else:
    print "Don't know what kind od tip you mean. I rather going to exit." ; exit()

#print "DEBUG: tc ", tc , " [s, px, py, pz, dz2, dxz, dyz ] "

# --- Reading geometry for plotting (if necessary)  --- #

if (plot_atoms):
    geom_plot, tmp1, tmp2 = Bu.loadAtoms('input.xyz'); del tmp1, tmp2;
    print "DEBUG: geom_plot", geom_plot
else:
    geom_plot = None

# --- the grid on which the STM signal is calculated --- #

if ((tip_type =='relaxed') or (tip_type == 'r')):
    print "Importing positions of PP from the PP-AFM calculations. Path for the data:"
    path_pos="Q%1.2fK%1.2f/" %(Q,K)
    print path_pos
    tip_r, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)
    extent = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[2,1])
    #print "DEBUG: extent", extent
    print "PP postions imported"
    dx=lvec[1,0]/(nDim[2]-1); dy=lvec[2,1]/(nDim[1]-1); dz=lvec[3,2]/(nDim[0]-1);
    tip_r0 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+lvec[1,0],dx,lvec[0,1],lvec[0,1]+lvec[2,1],dy,lvec[0,2],lvec[0,2]+lvec[3,2],dz)
    #print "DEBUG: dx, dy, dz", dx, dy, dz
    #print "DEBUG: tip_r.shape, tip_r0.shape", tip_r.shape, tip_r0.shape
else:
    print "Preparing the scan grid for fixed scan"
    extent = (x[0],x[1],y[0],y[1])
    tip_r  = RS.mkSpaceGrid(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2])
    lvec   = np.array([[x[0],y[0],z[0]],[x[1]-x[0],0.,0.],[0.,y[1]-y[0],0.],[0.,0.,z[1]-z[0]]])
    #print "DEBUG: extent", extent
    #print "DEBUG: lvec", lvec
    print "scan grids prepared"

# --- reading of the eigen-energies, the LCAO coefficients and geometry --- #

print "Reading electronic & geometry structure files"

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL') or (dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if isinstance(lvs, (list, tuple, np.ndarray)):
	cell = np.array([[lvs[0][0],lvs[0][1],0.0],[lvs[1][0],lvs[1][1],0.0],[0.0,0.0,99.9]]) if (len(lvs) == 2) else lvs
    elif isinstance(lvs, (str)):
	cell = np.loadtxt(lvs)
    elif ((pbc == (0,0)) or (pbc == (0.,0.))):
	cell = np.array([[0,0,0],[0,0,0],[0,0,0]]);
    else:
	print "PBC required, but lattice vector not specified. What can I do with that? I rather go to eat something."; exit()
    #print "DEBUG: cell.shape", cell.shape
# SHY ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL')):
#    if (kpointz == 'off'):
#       eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = files_path + kInputName, geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, kpoint = files_path+kpointFile, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#    elif (kpointz == 'on'):
       eigEn, coefsRe, coeffsIm, Ratin, gMode = RS.read_FIREBALL_all(name = files_path + kInputName, geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, kpoint = files_path+kpointFile, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);  
#    elif (mode == 'all'):
#       eigEn, coefs1, Ratin = RS.read_FIREBALL_real(name = files_path + kInputName, geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#       print "REAL COEFS", coefs1, coefs1.shape
#       eigEn, coefs2, Ratin = RS.read_FIREBALL_imag(name = files_path + kInputName, geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
#       print "IMAG COEFS", coefs2, coefs2.shape

#       coefs = np.concatenate((coefs1, coefs2), axis=1)
#       print "ALL COEFS", coefs, coefs.shape
    else:
       print "Please provide a file containing the k-points!"    

elif ((dft_code == 'gpaw') or(dft_code == 'GPAW')):
    eigEn, coefs, Ratin = RS.read_GPAW_all(    name = files_path + cp2k_name + '.gpw', fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

elif ((dft_code == 'aims') or(dft_code == 'AIMS') or (dft_code == 'FHI-AIMS')):
    if ((spin == None) or (spin == False)):
	name = 'KS_eigenvectors.band_1.kpt_1.out'
    elif ((spin == 'up')or(spin == 'alpha')or(spin == 'both')):
	name = 'KS_eigenvectors_up.band_1.kpt_1.out'
    elif ((spin == 'down')or(spin == 'beta')or(spin == 'dn')):
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
    else :
	print "unknown spin, I'm going to sleep. Good Night"; exit()
    eigEn, coefs, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    if (spin == 'both'):
	eigEn1 = eigEn.copy(); coefs1 = coefs.copy(); del eigEn, coefs ;
	name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
	eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)

elif ((dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if ((spin == None)or(spin == False)):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    elif ((spin == 'up')or(spin == 'alpha')):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
    elif (spin == 'both'):
	eigEn1, coefs1, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
	eigEn2, coefs2, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
	eigEn = np.concatenate((eigEn1, eigEn2), axis=0)
	coefs = np.concatenate((coefs1, coefs2), axis=0)
    elif ((spin == 'down')or(spin == 'beta')or(spin == 'dn')):
	eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
    else :
	print "unknown spin, I'm going to sleep. Good Night"; exit()

#print "DEBUG: eigEn.shape ", eigEn.shape
#print "DEBUG: coefs.shape ", coefs.shape
#print "DEBUG: Ratin.shape ", Ratin.shape

print "energies prepared, coeffecients read"
print "eigEn: ", eigEn

# --- the Main calculations --- #
# 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
'''
didv_b = False
STM_b  = False

if ( (scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
    didv    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
    didv_b = True; WF_decay= 0.0;
    #print "DEBUG: didv.shape ", didv.shape
elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
    nV = abs(V/dV)+1
    #print "DEBUG: V, nV:", V, nV
    current = np.array([   PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay) ])
    STM_b = True
    #print "DEBUG: current.shape ", current.shape
else:
    current, didv = PS.MSTM( V, V_max, dV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay)
    didv_b = True
    STM_b  = True
    #print "DEBUG: didv.shape ", didv.shape
    #print "DEBUG: current.shape ", current.shape
'''
# =========== Utils for plotting atoms =========================

def plotAtoms( atoms, atomSize=0.1, edge=True, ec='k', color='w' ):
    plt.fig = plt.gcf()
    es = atoms[0]; xs = atoms[1]; ys = atoms[2]
    for i in range(len(xs)):
	fc = '#%02x%02x%02x' % elements.ELEMENT_DICT[es[i]][7] #; print "DEBUG: fc", fc ; ##fc = '#FFFFFF' ##
	if not edge:
	    ec=fc
	circle=plt.Circle( ( xs[i], ys[i] ), atomSize, fc=fc, ec=ec  )
	plt.fig.gca().add_artist(circle)

def plotGeom( atoms=None, atomSize=0.1 ):
    if atoms is not None:
	plotAtoms( atoms, atomSize=atomSize )

#SHY
mat0 = genfromtxt('input.xyz', skip_header=2)
#data = pd.read_csv('input_plot.xyz', error_bad_lines=False)ss
data, target = np.array_split(np.loadtxt('input.xyz', skiprows=2, dtype=str), [-1], axis=1)

xs = mat0[:,1]
ys = mat0[:,2]
zs = mat0[:,3]
#print 'xs=', xs
#print 'zs=' ,zs

def findBonds(rmax):
   rr=float(rmax)
   bonds = []
   n = len( xs )
   for i in xrange(n):
      for j in xrange(i):
         dx=xs[j]-xs[i]
         dy=ys[j]-ys[i]
         dz=zs[j]-zs[i]
         r=math.sqrt(dx*dx+dy*dy+dz*dz)
         if (r<rr) :
            bonds.append( (i,j) )
#            plt.arrow(xs[i], ys[i], xs[j]-xs[i], ys[j]-ys[i], head_width=0.1, head_length=0.1,  fc='r', ec='C1', lw= 1.0,ls='solid',zorder=2 )
   return bonds
#   return plt.arrow(xs[i], ys[i], xs[j]-xs[i], ys[j]-ys[i], head_width=0.1, head_length=0.1,  fc='r', ec='C1', lw= 1.0,ls='solid',zorder=2 )

bonds = findBonds(1.59)
print "bonds are: " , bonds

def plotBonds(bb):
   for (q1,q2) in bb:
#     i=q1[0]; j=q2[1];
       plt.arrow(xs[q1], ys[q1], xs[q2]-xs[q1], ys[q2]-ys[q1], head_width=0.0, head_length=0.0,  fc='none', ec='C1', lw= 1.0,ls='dotted',zorder=2, alpha= 0.75)
     #  print (q1,q2)


def pltcolor(lst):
    cols=[]
    for l in lst:
        if l=='C':
            cols.append('g')
        elif l=='H':
            cols.append('r')
        else:
            cols.append('b')
    return cols

cols=pltcolor(data[:,0])
# --- plotting part here, plots all calculated signals --- #

#Voltages=np.arange(V,V_max+0.001,dV) # this part is important for scans over slabs at different voltages
namez = []
for V in Voltages:
    namez.append(str(round(V,2)))

didv_b = False
STM_b  = False

if ((scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
   if ((dft_code == 'fireball') or (dft_code == 'Fireball') or (dft_code == 'FIREBALL') and (gMode = 1)):
      didv    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefsRe, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
      didv_b = True; WF_decay= 0.0;
       #print "DEBUG: didv.shape ", didv.shape
   elif ((dft_code == 'fireball') or (dft_code == 'Fireball') or (dft_code == 'FIREBALL') and (gMode = 0)):
      didvR    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefsRe, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
      didvI    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefsIm, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
      didv     = didvR + didvI
      didv_b = True; WF_decay= 0.0;
   else:
      print "No modes selected, really? after all this time?"
       #print "DEBUG: didv.shape ", didv.shape

elif ( (scan_type == 'STM') or (scan_type == 'STM-single') ):
   nV = abs(V/dV)+1
   #print "DEBUG: V, nV:", V, nV
   current = np.array([   PS.STM( V, nV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay) ])
   STM_b = True
    #print "DEBUG: current.shape ", current.shape
else:
   current, didv = PS.MSTM( V, V_max, dV, WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6], WF_decay=WF_decay)
   didv_b = True
   STM_b  = True
  #print "DEBUG: didv.shape ", didv.shape
  #print "DEBUG: current.shape ", current.shape

NoV = len(didv) if didv_b else len(current)
NoH = len(didv[0]) if didv_b else len(current[0])
numVz = len(namez)

print "DEBUG: Voltages", Voltages
print "DEBUG: namez", namez
print "DEBUG: NoV", NoV
print "DEBUG: NoH", NoH
print "DEBUG: didv0", didv[0]
#print didv.shape


#
def splitter(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]


#
if PNG :
    print "We go to plotting "
    for vv in range(NoV):
	for k in range(NoH):
#            print "namez[vv]:", namez[vv]
	    #print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
	    name_plot=namez[vv]+';height: 5.0 $\AA$;tip:'+tip_type+';'+tip_orb+';'+mode
	    if didv_b :
		# ploting part here:
		plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
                plotBonds(bonds)
		plt.imshow(didv[vv,k,:,:], origin='image', extent=extent , cmap='gist_heat')
                tmp_curr=didv[vv,k,:,:].flatten()
                print "flattening done"
		out_curr=np.zeros((len(tmp_curr),3))
                print "writing to file"
		name_file = "./test/TEST"+tip_type+';'+tip_orb+';'+mode+namez[vv]+".dat"
		out_curr[:,2]=tmp_curr.copy()
		f=open(name_file,'w')
                print "file created"
                #print "tip0r", tip_r[k,:,:,0]
                print "tip2r", tip_r[k,:,:,2]
		out_curr[:,0]=tip_r[k,:,:,0].flatten()
		out_curr[:,1]=tip_r[k,:,:,1].flatten()
		print >> f, "WSxM file copyright Nanotec Electronica"
		print >> f, "WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al."
		print >> f, "X[A]  Y[A]  Z[A]"
		print >> f, ""
		np.savetxt(f, out_curr)
		f.close()

		#f=open(name_file,'w')
		#np.savetxt(f, out_curr)
		#f.close()
               # np.savetxt('curTEST'+str(k)+'.txt', didv[vv,k,:,:], fmt='%.18f', delimiter=' ', newline=os.linesep)
		plotGeom(atoms=geom_plot)
		plt.xlabel(r' Tip_x $\AA$')
		plt.ylabel(r' Tip_y $\AA$')
		plt.title("dIdV:"+name_plot)
		plt.savefig( './test/TEST_K02_I_Stip_5.0_TEST4'+'_'+mode+"-"+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
		plt.close()
	    if STM_b :
		# ploting part here:
		plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
		plt.imshow(current[vv,k,:,:], origin='image', extent=extent , cmap='gray')
		plotGeom(atoms=geom_plot)
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
		name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.dat' %k 
		tmp_curr=didv[vv,k,:,:].flatten()
		out_curr=np.zeros((len(tmp_curr),3))
		out_curr[:,0]=tip_r0[k,:,:,0].flatten()
		out_curr[:,1]=tip_r0[k,:,:,1].flatten()
		out_curr[:,2]=tmp_curr.copy()
		f=open(name_file,'w')
		print >> f, "WSxM file copyright Nanotec Electronica"
		print >> f, "WSxM ASCII DAT file; obtained from dIdV code by Krejci et al."
		print >> f, "X[A]  Y[A]  Z[A]"
		print >> f, ""
		np.savetxt(f, out_curr)
		f.close()
		#
	    if STM_b :
		name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'_%03d.xyz' %k 
		tmp_curr=current[vv,k,:,:].flatten()
		out_curr=np.zeros((len(tmp_curr),3))
		out_curr[:,0]=tip_r0[k,:,:,0].flatten()
		out_curr[:,1]=tip_r0[k,:,:,1].flatten()
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
    xsf_head = Bu.At2XSF(geom_plot) if plot_atoms else GU.XSF_HEAD_DEFAULT
    for vv in range(NoV):
	if didv_b :
	    name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'.xsf'
	    GU.saveXSF(name_file, didv[vv], lvec, head=xsf_head )
	if STM_b :
	    name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'.xsf'
	    GU.saveXSF(name_file, current[vv], lvec, head=xsf_head )
    print "XSF files written"

if NPY :
    print "writing npy binary files"
    for vv in range(NoV):
	if didv_b :
	    name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)
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
