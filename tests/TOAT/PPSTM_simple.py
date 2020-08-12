#!/usr/bin/python3
#
##########################################################################################################################
#                                                                                                                        #
#                                  What follows are options of the PP-STM (dIdV) code                                    #
#                                                                                                                        #
##########################################################################################################################
#
# Note : This type of simulations works for solid slabs or molecules on slabs (substrate) ; for freestanding molecule it can give you nonsences
#
# ***** System information: *****
#
ppstm_path = './PPSTM/'      # path (absolute or relative) to your PPSTM code #
#
ncpu       = 1               # number of cpu cores for OMP paralelization: ncpu = 1 -- serial compilation & calculations; iff ncpu > 1, then OMP paralel recompilation is used and C++ calculations are running on more cores #
#
# ***** Main informations ******
#
scan_type     = 'v-scan'     # 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #
tip_type      = 'fixed'      # 'fixed'='f' -- for stiff/metal tip apexes ; 'relaxed'='r' -- for flexible tip apexes (precalculated by PP-AFM) . For this option you have to have "installed" PPAFM in your PPSTM directory #
V             = -0.5         # !!!! V = Vmin for SCAN !!!! #
V_max         = +0.5         # V = V_min >= -2.0 V ; V_max <= 2.0 V (othervise changes in the later code needed) #
dV            =  0.1         # voltage step , dV <= 0.1 V #
eta           =  0.1         # Lorentzian width of states in energy scale: typically 0.1; can be in range of 0.3-0.05 eV in some cases (low amount of layers ...) even up to 1.0 eV #
WF_decay      =  0.5         # 0.0 <= WF_decay <= 1.0 ; How fast WorkFunction tunnelling barrier is changing with Voltage : (WF = WF_0 + V*WF_decay) -- 0.0 no change ; 1.0 - the same change as voltage #
tip_orb       = 's'          # 's' ; 'pxy' -- px & py ; 'spxy' -- 50% s & 50% pxy ; '5spxy' -- 5% s & 95% pxy ; '10spxy' -- 10% s & 90% pxy ; 'CO' -- 13% s & 87% pxy (PRL 119, 166001 (2017)) ; 'pz' ; For sample_orbs = 'sp' , possible 'dz2' and 'dxzyz' -- dxz & dyz #
sample_orbs   = 'sp'         # orbitals of the sample 'sp' (light atoms only, faster) or 'spd' (all atoms) #
dft_code      = 'cp2k'       # 'fireball'='Fireball'='FIREBALL' ; 'aims'='AIMS'='FHI-AIMS' ; 'cp2k'='CP2K' ; 'gpaw'='GPAW' #
geometry_file = 'TOAT.xyz'   # E.G. 'input.xyz' , 'input.bas' , 'geometry.in'; None for GPAW #
pbc           = (0,0)      # (0,0) = None = False -- only original geometry ; (0.5,0.5) -- 2x2 cell ; (1,1) -- 3x3 cell (around original) ; (2,2) -- 5x5 cell (around original) ... #
lvs           =  None        # None ; [[ax,ay,0],[bx,by,0]],[0,0,cz]] or [[ax,ay],[bx,by]] ; 'input.lvs' -- files with specified cell ; in FHI-AIMS & GPAW allready specified with geometry #
spin          =  None        # None=False ; for FHI-AIMS & CP2K: None -- spin-unpolarized/spin-restricted calc. ;  'both' , 'up'='alpha' or 'down" (last 3 are for spin-polarizes or spin-unrestricted calculations) #
cp2k_name     = 'TOAT'       # Name used in CP2K calculations or GPAW calc #
#
# ***** Informations for x,y,z tip_type = 'fixed' ******
#
x = [ -3.0, 15.0, 0.10 ]     # [xmin, xmax, dx] #
y = [ -3.0, 15.0, 0.10 ]     # [ymin, ymax, dy] #
z = [ 12.7, 13.7, 1.0  ]     # !!!! z-starts from zero - normally zmin >= 3-4 Ang above heighest atoms !!!! [zmin, zmax, dz] ; for single height scan use : [z, z, 1.0] #
#
# ***** Informations for PP positions, tip_type = 'relaxed' ******
#
Q = 0.00                     # charge (PP-AFM) ; Ocharge PP-AFM complex_tip autumn 2018) ; [e] (monopole), [e*A] (dipole), [e*A^2] (quadrupole) #
K = 0.24                     # x stiffness (PP-AFM master autumn 2018); klat (PP-AFM dev/OpenCl autumn 2018); Oklat (PP-AFM complex_tip autumn 2018) ; [N/m] #
data_format = 'npy'          # 'xsf'='XSF' ; 'npy'='NPY' ; -- format in which PPpos are stored from PP-AFM run #
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
cut_atoms   = -1             # None = -1 -- All atoms of the sample contributes to tunelling ; 1 -- only 1st atom of the sample contributes to the tunelling ; 57 -- first 57 atoms of the sample contributes to the tunelling ; ... #
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

print("Importing libraries")

import os
import sys
sys.path.append(ppstm_path) 

if (ncpu > 1):
    os.environ['OMP_NUM_THREADS'] = str(ncpu)
    print('OMP_NUM_THREADS:', os.environ['OMP_NUM_THREADS'])

import numpy as np
import pyPPSTM                   as PS
import pyPPSTM.ReadSTM           as RS
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend. ## !!! important for working on clusters !!!!
import matplotlib.pyplot as plt
if (XSF or NPY or (tip_type == 'relaxed') or (tip_type == 'r' )):
    print("For XSF or NPY outputs or tip_type = relaxed you have to have installed PPAFM in your PPSTM directory ")
    import pyProbeParticle.GridUtils as GU
if (plot_atoms):
    import pyPPSTM.basUtils as Bu
    import pyPPSTM.elements as elements


print("Libraries imported")

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
    print("Don't know what kind od tip you mean. I rather going to exit.") ; exit()

#print "DEBUG: tc ", tc , " [s, px, py, pz, dz2, dxz, dyz ] "

# --- Reading geometry for plotting (if necessary)  --- #

if (plot_atoms):
    geom_plot, tmp1, tmp2 = Bu.loadAtoms('input_plot.xyz'); del tmp1, tmp2;
    #print "DEBUG: geom_plot", geom_plot
else:
    geom_plot = None

# --- the grid on which the STM signal is calculated --- #

if ((tip_type =='relaxed') or (tip_type == 'r')):
    print("Importing positions of PP from the PP-AFM calculations. Path for the data:")
    path_pos="Q%1.2fK%1.2f/" %(Q,K)
    print(path_pos)
    tip_r, lvec, nDim = GU.load_vec_field( path_pos+'PPpos' ,data_format=data_format)
    extent = (lvec[0,0],lvec[0,0]+lvec[1,0],lvec[0,1],lvec[0,1]+lvec[2,1])
    #print "DEBUG: extent", extent
    print("PP postions imported")
    dx=lvec[1,0]/(nDim[2]-1); dy=lvec[2,1]/(nDim[1]-1); dz=lvec[3,2]/(nDim[0]-1);
    tip_r0 = RS.mkSpaceGrid(lvec[0,0],lvec[0,0]+lvec[1,0],dx,lvec[0,1],lvec[0,1]+lvec[2,1],dy,lvec[0,2],lvec[0,2]+lvec[3,2],dz)
    #print "DEBUG: dx, dy, dz", dx, dy, dz
    #print "DEBUG: tip_r.shape, tip_r0.shape", tip_r.shape, tip_r0.shape
else:
    print("Priparing the scan grid for fixed scan")
    extent = (x[0],x[1],y[0],y[1])
    tip_r  = RS.mkSpaceGrid(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2])
    lvec   = np.array([[x[0],y[0],z[0]],[x[1]-x[0],0.,0.],[0.,y[1]-y[0],0.],[0.,0.,z[1]-z[0]]])
    #print "DEBUG: extent", extent
    #print "DEBUG: lvec", lvec
    tip_r0 = tip_r
    print("scan grids prepared")

# --- reading of the eigen-energies, the LCAO coefficients and geometry --- #

print("Reading electronic & geometry structure files")

if ((dft_code == 'fireball') or(dft_code == 'Fireball') or (dft_code == 'FIREBALL') or (dft_code == 'cp2k') or(dft_code == 'CP2K')):
    if isinstance(lvs, (list, tuple, np.ndarray)):
        cell = np.array([[lvs[0][0],lvs[0][1],0.0],[lvs[1][0],lvs[1][1],0.0],[0.0,0.0,99.9]]) if (len(lvs) == 2) else lvs
    elif isinstance(lvs, (str)):
        cell = np.loadtxt(lvs)
    elif ((pbc == (0,0)) or (pbc == (0.,0.))):
        cell = np.array([[0,0,0],[0,0,0],[0,0,0]]);
    else:
        print("PBC required, but lattice vector not specified. What can I do with that? I rather go to eat something."); exit()
    #print "DEBUG: cell.shape", cell.shape

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
        print("unknown spin, I'm going to sleep. Good Night"); exit()

#print "DEBUG: eigEn.shape ", eigEn.shape
#print "DEBUG: coefs.shape ", coefs.shape
#print "DEBUG: Ratin.shape ", Ratin.shape

print("energies prepared, coeffecients read")

# --- the Main calculations --- #
# 'didv'='dIdV''='didv-single' -- only dIdV for one voltage = V ; 'v-scan'='V-scan'='Voltage-scan' -- both STM & dIdV scan - V .. Vmax; 'STM'='STM-single' -- STM for one Voltage = V, use V-scan rather #

didv_b = False
STM_b  = False
states_b = False

if ( (scan_type == 'didv') or (scan_type == 'dIdV') or (scan_type == 'didv-single')):
    didv    = np.array([   PS.dIdV( V,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
    didv_b = True; WF_decay= 0.0;
    #print "DEBUG: didv.shape ", didv.shape
elif ( (scan_type == 'states') or (scan_type == 'STATES') ):
    states = np.sort(eigEn); mask = states >= V; states = states[mask]; del mask;
    mask = states <= V_max; states = states[mask]; del mask;
    fst = True
    print("DEBUG: states:", states)
    for isi in states:
        if fst:
            didv    = np.array([   PS.dIdV( isi,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ])
            fst = False
        else :
            didv =np.append(didv, [PS.dIdV( isi,    WorkFunction, eta, eigEn, tip_r, Ratin, coefs, orbs=sample_orbs, s=tc[0], px =tc[1], py=tc[2], pz=tc[3], dz2=tc[4], dxz=tc[5], dyz=tc[6] ) ], axis =0)
        #print "DEBUG: didv.shape ", didv.shape
        didv_b = True; states_b = True; WF_decay= 0.0;

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


# --- plotting part here, plots all calculated signals --- #

Voltages= np.arange(V,V_max+0.001,dV) if not states_b else states # this part is important for scans over slabs at different voltages
round_index = 2 if not states_b else 5
print("Voltages", Voltages)
namez = []
for V in Voltages:
    namez.append(str(round(V,round_index)))

NoV = len(didv) if didv_b else len(current)
NoH = len(didv[0]) if didv_b else len(current[0])

#print "DEBUG: Voltages", Voltages
#print "DEBUG: namez", namez
#print "DEBUG: NoV", NoV
#print "DEBUG: NoH", NoH

if PNG :
    print("We go to plotting ")
    for vv in range(NoV):
        for k in range(NoH):
            #print "DEBUG: long name:::", namez[vv],';height:%03d;tip:'  %k,tip_type,';',tip_orb
            name_plot=namez[vv]+';height:'+str(k)+';tip:'+tip_type+';'+tip_orb
            if didv_b :
                # ploting part here:
                plt.figure( figsize=(0.5 * lvec[1,0] , 0.5 * lvec[2,1] ) )
                plt.imshow(didv[vv,k,:,:], origin='image', extent=extent , cmap='gray')
                plotGeom(atoms=geom_plot)
                plt.xlabel(r' Tip_x $\AA$')
                plt.ylabel(r' Tip_y $\AA$')
                plt.title("dIdV:"+name_plot)
                plt.savefig( 'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.png' %k , bbox_inches='tight' )
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
    print("Everything plotted")
if WSxM :
    print("writing WSxM files")
    for vv in range(NoV):
        for k in range(NoH):
            if didv_b :
                name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'_%03d.xyz' %k 
                tmp_curr=didv[vv,k,:,:].flatten()
                out_curr=np.zeros((len(tmp_curr),3))
                out_curr[:,0]=tip_r0[k,:,:,0].flatten()
                out_curr[:,1]=tip_r0[k,:,:,1].flatten()
                out_curr[:,2]=tmp_curr.copy()
                f=open(name_file,'w')
                print("WSxM file copyright Nanotec Electronica", file=f)
                print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                print("X[A]  Y[A]  Z[A]", file=f)
                print("", file=f)
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
                print("WSxM file copyright Nanotec Electronica", file=f)
                print("WSxM ASCII XYZ file; obtained from dIdV code by Krejci et al.", file=f)
                print("X[A]  Y[A]  Z[A]", file=f)
                print("", file=f)
                np.savetxt(f, out_curr)
                f.close()
        #
    print("WSxM files written")

if XSF :
    print("writing XSF files")
    xsf_head = Bu.At2XSF(geom_plot) if plot_atoms else GU.XSF_HEAD_DEFAULT
    for vv in range(NoV):
        if didv_b :
            name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)+'.xsf'
            GU.saveXSF(name_file, didv[vv], lvec, head=xsf_head )
        if STM_b :
            name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)+'.xsf'
            GU.saveXSF(name_file, current[vv], lvec, head=xsf_head )
    print("XSF files written")

if NPY :
    print("writing npy binary files")
    for vv in range(NoV):
        if didv_b :
            name_file =  'didv_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction-Voltages[vv]*WF_decay)+"_eta_"+str(eta)
            GU.saveNpy(name_file, didv[vv], lvec)#, head=XSF_HEAD_DEFAULT )
        if STM_b :
            name_file =  'STM_'+namez[vv]+"_tip_"+tip_type+"-"+tip_orb+"_WF_"+str(WorkFunction)+"_WF_decay_"+str(round(WF_decay,1))+"_eta_"+str(eta)
            GU.saveNpy(name_file, current[vv], lvec)#, head=XSF_HEAD_DEFAULT )
    print("npy files written")

# --- the end --- #

print() 
print()
print("Done")
print()

