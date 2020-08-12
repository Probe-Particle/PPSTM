#!/usr/bin/python
#
##########################################################################################################################
#                                                                                                                        #
#                   What follows are options of the PP-STM code for calculating inner PDOS                               #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** System information: *****
#
ppstm_path = './PPSTM/'      # path (absolute or relative) to your PPSTM code #
#
# ***** Main informations ******
#
V_min         = -2.0         # V_min for plotting
V_max         = +2.0         # V_max for plotting #
eta           =  0.1         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
dV            =  eta/5.0     # voltage/energy step , default : eta/5.0  #
sample_orbs   = 'sp'         # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
dft_code      = 'CP2K'       # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
geometry_file = 'TOAT.xyz'   # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW #
spin          =  None        # None=False ; for FHI-AIMS & CP2K: None -- spin-unpolarized/spin-restricted calc. ;  'both'--  up will be plotted as + and down as - , 'up'='alpha' or 'down"='beta' (last 3 are for spin-polarizes or spin-unrestricted calculations) #
cp2k_name     = 'TOAT'       # Name used in CP2K calculations or GPAW calc #
#
# *****Plotting options ******
#
red_line_atoms    = [0]      # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; [0] -- means 1st atom #
red_line_shell    = 's'      # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
blue_line_atoms   = [1]      # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; [1] -- means 2nd atom #
blue_line_shell   = 'pz'     # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
green_line_atoms  = range(2,8)     # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; range(2,8) = [2,3,4,5,6,7]-- means atoms: 3,4,5,6,7 and 8 #
green_line_shell  = 'p'      # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
orange_line_atoms = [8,9]    # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; [8,9] -- means 9th and 10th atom #
orange_line_shell = 'px'     # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
black_line_atoms  = 'all'    # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; None -- no line #
black_line_shell  = 'all'    # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
yellow_line_atoms = None     # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; None -- no line #
yellow_line_shell = 'all'    # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
gray_line_atoms   = None     # atoms to plot: !!! DO NOT FORGET PYTHON NUMBERING starts from 0 !!!, always use list-like structure; 'all' = -1 -- all atoms; None -- no line #
gray_line_shell   = 'all'    # which chells to plot: 'all' = all sample orbitals, 's', 'p', 'd', 'px', 'py', 'pz', 'dxy', 'dyz', 'dyz', 'dz2', 'dxz', 'dx2y2' #
#
# *****Output options ******
#
PNG  = True                  # True / False -- plot "png" images (2D graph height) #
TXT = False                  # **** Not working yet *** True / False -- write ".txt" files with 1st column energy, other columns PDOS#
#
# ***** Advanced options ******
#
cut_atoms   = None           # Not really neccessary for these calculations, unless really big files and memory; None = -1 -- All atoms of the sample contributes to tunelling ; 1 -- only 1st atom of the sample contributes to the tunelling ; 57 -- first 57 atoms of the sample contributes to the tunelling ; ... #
#
# ***** More advanced options ******
#
fermi        = None          # None=0.0 -- no change to the Fermi Level ; -0.1 -- shifts the Fermi Level by 0.1 eV lower ... #
cut_min      = V_min-3*eta   # cut out all orbitals lower than  -default: V_min-3*eta bellow Fermi (should be: cut_min <= Vmin-2*eta) . taken to the Fermi Level #
cut_max      = V_max+3*eta   # cut out all orbitals higher than -default: V_max+3*eta above  Fermi (should be: cut_max >= Vmax+2*eta) . taken to the Fermi Level #
files_path   = ''            # where are files fron DFT code ; rather do not use this #
lower_atoms  = 'no-d-rescalling' # normally d-orbs are rescalled by factor of 0.2 #
#
#
##########################################################################################################################
#                                                                                                                        #
#                                 DO NOT TOUCH LATER CODE (unless you know what you are doing)                           #
#                                                                                                                        #
##########################################################################################################################

print("Importing libraries")

import os
import sys
sys.path.append(ppstm_path) 

import numpy as np
#import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import pyPPSTM.PreSTMutils       as SU
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
#if (plot_atoms):
#    import pyPPSTM.basUtils as Bu
#    import pyPPSTM.elements as elements


print("Libraries imported")

# --- some function definition --- #

def printf(*args):
    together = ''.join(map(str, args))    # avoid the arg is not str
    #print together
    return together

# --- Initial check --- #

assert( PNG or TXT ), "No output set to be True; I'm not going to do anything if there is no output. I'm too lazy like a Gartfield. "

# --- reading of the eigen-energies, the LCAO coefficients and geometry --- #

print("Reading electronic & geometry structure files")

cell=[[0,0],[0,0]];pbc=(0,0);lower_coefs=[];

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL')):
    eigEn, coefs, Ratin = RS.read_FIREBALL_all(name = files_path + 'phik_0001_', geom=files_path+geometry_file, lvs = cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max,cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);

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
        print("unknown spin, I'm going to sleep. Good Night"); exit()
    eigEn, coefs, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    if (spin == 'both'):
        name = 'KS_eigenvectors_dn.band_1.kpt_1.out'
        eigEn2, coefs2, Ratin = RS.read_AIMS_all(name = files_path + name , geom= files_path + geometry_file, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
        #coefs2 *= -1;

elif ((dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if ((spin == None)or(spin == False)):
        eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs);
    elif ((spin == 'up')or(spin == 'alpha')):
        eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
    elif (spin == 'both'):
        eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='alpha');
        eigEn2, coefs2, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
        #coefs2 *= -1;
    elif ((spin == 'down')or(spin == 'beta')or(spin == 'dn')):
        eigEn, coefs, Ratin  = RS.read_CP2K_all(name = files_path + cp2k_name , lvs=cell, fermi=fermi, orbs = sample_orbs, pbc=pbc, cut_min=cut_min, cut_max=cut_max, cut_at=cut_atoms, lower_atoms=lower_atoms, lower_coefs=lower_coefs, spin='beta');
    else :
        print("unknown spin, I'm going to sleep. Good Night"); exit()

#print "DEBUG: eigEn.shape ", eigEn.shape
#print "DEBUG: coefs.shape ", coefs.shape
#print "DEBUG: Ratin.shape ", Ratin.shape

energies = np.arange(V_min,V_max,dV)

print("energies prepared, coeffecients read")

# --- the PDOS calculations --- #

plot = True if ( (red_line_atoms != None)or(blue_line_atoms != None)or(green_line_atoms != None)or(orange_line_atoms != None)or(black_line_atoms != None)or(yellow_line_atoms != None)or (gray_line_atoms != None) ) else False

print("DEBUG: plot", plot)

assert plot!=False, "No lines to plot"

latomslegend = 4

if (red_line_atoms != None) :
    atoms=red_line_atoms
    PDOS1 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms   , spherical=red_line_shell   )
    red_line_legend     = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
    print("DEBUG: red_line_legend", red_line_legend)
else:	#(red_line_atoms != None) :
    PDOS1 = None

if (blue_line_atoms != None) :
    atoms=blue_line_atoms
    PDOS2 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms  , spherical=blue_line_shell  )
    blue_line_legend    = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(blue_line_atoms != None) :
    PDOS2 = None

if (green_line_atoms != None) :
    atoms=green_line_atoms
    PDOS3 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms , spherical=green_line_shell )
    green_line_legend   = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(green_line_atoms != None) :
    PDOS3 = None

if (orange_line_atoms != None) :
    atoms=orange_line_atoms
    PDOS4 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms, spherical=orange_line_shell)
    orange_line_legend  = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(orange_line_atoms != None) :
    PDOS4 = None

if (black_line_atoms != None) :
    atoms=black_line_atoms
    PDOS5 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms , spherical=black_line_shell )
    black_line_legend   = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(black_line_atoms != None) :
    PDOS5 = None

if (yellow_line_atoms != None) :
    atoms=yellow_line_atoms
    PDOS6 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms, spherical=yellow_line_shell)
    yellow_line_legend  = ','.join(map(str, atoms)) if len(atoms) < latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(yellow_line_atoms != None) :
    PDOS6 = None

if (gray_line_atoms != None) :
    atoms=gray_line_atoms
    PDOS7 = SU.pPDOS(eigEn,coefs, energies, eta=eta, orbs= sample_orbs, atoms=atoms  , spherical=gray_line_shell  )
    gray_line_legend    = ','.join(map(str, atoms)) if len(atoms) > latomslegend else str(atoms[0])+".."+str(atoms[-1])
else:	#(gray_line_atoms != None) :
    PDOS7 = None

if spin=='both' :
    print("DEBUG: printing spin down in both spins")

    PDOS1d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=red_line_atoms   , spherical=red_line_shell   ) if (red_line_atoms != None) else None
    PDOS2d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=blue_line_atoms  , spherical=blue_line_shell  ) if (blue_line_atoms != None) else None
    PDOS3d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=green_line_atoms , spherical=green_line_shell ) if (green_line_atoms != None) else None
    PDOS4d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=orange_line_atoms, spherical=orange_line_shell) if (orange_line_atoms != None) else None
    PDOS5d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=black_line_atoms , spherical=black_line_shell ) if (black_line_atoms != None) else None
    PDOS6d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=yellow_line_atoms, spherical=yellow_line_shell) if (yellow_line_atoms != None) else None
    PDOS7d = -1*SU.pPDOS(eigEn2,coefs2, energies, eta=eta, orbs= sample_orbs, atoms=gray_line_atoms  , spherical=gray_line_shell  ) if (gray_line_atoms != None) else None

else: 
    print("DEBUG: both spins are not used")
    PDOS1d = PDOS2d = PDOS3d = PDOS4d = PDOS5d = PDOS6d = PDOS7d = None

# --- plotting part here, plots all calculated signals --- #

if PDOS1 is not None:
    plt.plot(energies, PDOS1, ls='-',c='r'      ,label= 'PDOS atoms:'+red_line_legend+";"+red_line_shell+"-shell")
if PDOS2 is not None:
    plt.plot(energies, PDOS2, ls='-',c='g'      ,label= 'PDOS atoms:'+blue_line_legend+";"+blue_line_shell+"-shell")
if PDOS3 is not None:
    plt.plot(energies, PDOS3, ls='-',c='b'      ,label= 'PDOS atoms:'+green_line_legend+";"+green_line_shell+"-shell")
if PDOS4 is not None:
    plt.plot(energies, PDOS4, ls='-',c='#FFA500',label= 'PDOS atoms:'+orange_line_legend+";"+orange_line_shell+"-shell")
if PDOS5 is not None:
    plt.plot(energies, PDOS5, ls='-',c='k'      ,label= 'PDOS atoms:'+black_line_legend+";"+black_line_shell+"-shell")
if PDOS6 is not None:
    plt.plot(energies, PDOS6, ls='-',c='y'      ,label= 'PDOS atoms:'+yellow_line_legend+";"+yellow_line_shell+"-shell")
if PDOS7 is not None:
    plt.plot(energies, PDOS7, ls='-',c='g'      ,label= 'PDOS atoms:'+gray_line_legend+";"+gray_line_shell+"-shell")

# spin down if "both"
if PDOS1d is not None:
    plt.plot(energies, PDOS1d, ls='--',c='r')#,label= 'PDOS atoms:'+str(red_line_atoms)+";"+red_line_shell+"shell")
if PDOS2d is not None:
    plt.plot(energies, PDOS2d, ls='--',c='g')#,label= 'PDOS atoms:'+str(blue_line_atoms)+";"+red_line_shell+"shell")
if PDOS3d is not None:
    plt.plot(energies, PDOS3d, ls='--',c='b')#,label= 'PDOS atoms:'+str(green_line_atoms)+";"+red_line_shell+"shell")
if PDOS4d is not None:
    plt.plot(energies, PDOS4d, ls='--',c='#FFA500')#,label= 'PDOS atoms:'+str(orange_line_atoms)+";"+red_line_shell+"shell")
if PDOS5d is not None:
    plt.plot(energies, PDOS5d, ls='--',c='k')#,label= 'PDOS atoms:'+str(red_line_atoms)+";"+red_line_shell+"shell")
if PDOS6d is not None:
    plt.plot(energies, PDOS6d, ls='--',c='y')#,label= 'PDOS atoms:'+str(red_line_atoms)+";"+red_line_shell+"shell")
if PDOS7d is not None:
    plt.plot(energies, PDOS7d, ls='--',c='g')#,label= 'PDOS atoms:'+str(red_line_atoms)+";"+red_line_shell+"shell")

if (plot and PNG) :
    plt.title("PDOS eta: "+str(eta)+' eV')
    plt.xlabel("E-Efermi [eV]")
    plt.ylabel("DOS [arb.un]")
    plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
    plt.savefig('PDOS_eta_'+str(eta)+'eV_orbs_'+sample_orbs+'.png', bbox_inches='tight')


# --- the end --- #

print() 
print()
print("Done")
print()

