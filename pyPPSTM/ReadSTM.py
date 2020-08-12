#!/usr/bin/python

import os
import sys
import numpy as np
from . import basUtils as bU
from . import elements

from   ctypes import c_int, c_double, c_char_p
import ctypes
from . import cpp_utils

import time

# this library has functions for reading STM coefficients and make a grid for non-relaxed 3D scan

# global variables:

cut_at_ =-1
pbc_ = (0,0)

lower_atoms_ = []
lower_coefs_ = []
cut_min_ = -15.0
cut_max_ = 5.0

n_min_ = 0
n_max_ = 0
Ynum_ = 4
num_at_ = -1

dOrbRes_ = True

# ==============================
# ============================== C++ compilations
# ==============================

LIB_PATH = os.path.dirname( os.path.realpath(__file__) )
print(" ProbeParticle Library DIR = ", LIB_PATH)

cpp_name='IO'
#cpp_utils.compile_lib( cpp_name  )
make_name='MIO' if sys.platform=='darwin' else 'IO'
print("DEBUG: make_name", make_name)
cpp_utils.make(make_name)
lib    = ctypes.CDLL(  cpp_utils.CPP_PATH + "/" + cpp_name + cpp_utils.lib_ext )     # load dynamic librady object using ctypes 

# define used numpy array types for interfacing with C++

array1i = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
array4d = np.ctypeslib.ndpointer(dtype=np.double, ndim=4, flags='CONTIGUOUS')

# ========
# ======== Python warper function for C++ functions
# ========

#************* sp(d) now as well *************
# int read_AIMS_coefs          (char *fname, double* coefs, int* period, int nMOmax, int nMOmin, int nAtoms, int nPerAtoms ){
lib.read_AIMS_coefs.argtypes = [ c_char_p,   array3d,       array1i,     c_int,      c_int,      c_int,      c_int ]
lib.read_AIMS_coefs.restype  = c_int
def read_AIMS_coefs(fname, at_nums ):
    #eigs = ReadSTM.getAimsEigenE(fname)
    #nMO = len(eigs)
    periods = np.array([ elements.ELEMENTS[iZ][2] for iZ in at_nums ], dtype=np.int32)
    #for iZ,per in zip(at_nums,periods): print iZ,per # for DEBUG
    coefs = np.zeros( (n_max_ - n_min_ ,num_at_ , Ynum_) ); # print "DEBUG: coefs.shape", coefs.shape
    lib.read_AIMS_coefs( str.encode(fname),coefs, periods, n_max_, n_min_ , num_at_ , Ynum_ );
    return coefs.copy()#, eigs.copy()

# ==============================
# ============================== Pure python functions
# ==============================


def mkSpaceGrid(xmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz):
    '''
    mkSpaceGridsxmin,xmax,dx,ymin,ymax,dy,zmin,zmax,dz):
    Give rectangular grid along the main cartesian axes for non-relaxed dI/dV or STM - 4D grid of xyz coordinates.
        '''
    h = np.mgrid[xmin:xmax+0.0001:dx,ymin:ymax+0.0001:dy,zmin:zmax+0.0001:dz]
    f = np.transpose(h)
    sh = f.shape
    print("Grid has dimensios: ", sh)
    return f;	#, sh;

# preparing procedures:

def initial_check(orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
    '''
    do some initial checks of incoming parameters (orbs, imaginary) and most of the global parameters for reading procedures
    '''
    assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
    assert (imaginary == False), "sorry imaginary version is under development" 	
    print("reading FHI-AIMS LCAO coefficients for basis: ",orbs)	
    global cut_at_ ; cut_at_ = -1 if cut_at == None else cut_at
    global pbc_    ; pbc_ = pbc
    global cut_min_ ; cut_min_ = cut_min
    global cut_max_ ; cut_max_ = cut_max
    global Ynum_    ; Ynum_ = 4 if (orbs =='sp') else 9
    global lower_atoms_ ; lower_atoms_ = lower_atoms if lower_atoms != 'no-d-rescalling' else []
    #print  "DEBUG: lower_atoms_", lower_atoms_
    global dOrbRes_     ; dOrbRes_     = True        if lower_atoms != 'no-d-rescalling' else False
    #print  "DEBUG: dOrbRes_", dOrbRes_
    global lower_coefs_ ; lower_coefs_ = lower_coefs

# geometries underprocedures: 

def cut_atoms(atoms):
    '''
    Cut unwanted atoms !!! from the end of the geometry !!! important atoms should be first
    '''
    assert (cut_at_ <= len(atoms[1])), "wrong cut for atoms"
    if not ((cut_at_ == -1)or(cut_at_ == len(atoms[1]))):
        atoms2 = [atoms[0][:cut_at_],atoms[1][:cut_at_],atoms[2][:cut_at_],atoms[3][:cut_at_]]
    else:
        atoms2 = atoms
    global num_at_ ; num_at_= len(atoms2[1])
    return atoms2;

def for_PBC(atoms,lvs):
    '''
    Apply PBC onto the geometry
    '''
    if (pbc_ != ((0,0)or(0.,0.))):
        assert (lvs != (None or []) ), "Lattice vectors (cell) not specified"
        print("Applying PBC")
        if (pbc_ == (0.5,0.5)):
            atoms = bU.multCell( atoms, lvs, m=(2,2,1) )
            Rs = np.array([atoms[1],atoms[2],atoms[3]])
        else:
            atoms = bU.multCell( atoms, lvs, m=( (int(2*pbc_[0])+1),(int(2*pbc_[1])+1),1 ) )
            Rs = np.array([atoms[1],atoms[2],atoms[3]]); 
            Rs[0] -= int(pbc_[0])*lvs[0,0]+int(pbc_[1])*lvs[1,0]
            Rs[1] -= int(pbc_[0])*lvs[0,1]+int(pbc_[1])*lvs[1,1]
        print(" Number of atoms after PBC: ", len(Rs[0]))
    else:
        Rs = np.array([atoms[1],atoms[2],atoms[3]])
    Ratin    = np.transpose(Rs).copy()
    return Ratin

# procedures for preparing geometries for STM:

def get_FIREBALL_geom(geom='answer.bas', lvs=None, sl=False):
    '''
    Prepares geometry from the FIREBALL files format, sl = means automatically skip line in the xyz file
    '''
    print(" # ============ define atoms ============")
    atoms, nDim, tmp = bU.loadAtoms(geom, sl=sl)
    #print "DEBUG: atoms", atoms
    del nDim, tmp
    atoms = cut_atoms(atoms)
    print(" Number of atoms: ", num_at_)
    Ratin = for_PBC(atoms,lvs)
    print("atomic geometry read")
    return Ratin ;

def get_AIMS_geom(geom='geometry.in'):
    '''
    Prepares geometry from the FHI-AIMS files format
    '''
    print(" # ============ define atoms ============")
    atoms, nDim, lvs = bU.loadGeometryIN(geom)
    lvs = np.array(lvs)
    del nDim
    atoms = cut_atoms(atoms)
    at_num = []
    for i in atoms[0]:
        at_num.append(elements.ELEMENT_DICT[i][0])
    print(" Number of atoms: ", num_at_)
    Ratin = for_PBC(atoms,lvs)
    print("atomic geometry read")
    return Ratin, np.array(at_num);

def get_GPAW_geom(geom=None):
    '''
    Prepares geometry from the ASE atoms binary
    '''
    print(" # ============ define atoms ============")
    from ase import Atoms
    tmp = geom.get_positions()
    atoms = [geom.get_atomic_numbers(), tmp[:,0], tmp[:,1], tmp[:,2] ]
    lvs = geom.get_cell()
    del tmp
    atoms = cut_atoms(atoms)
    print(" Number of atoms: ", num_at_)
    Ratin = for_PBC(atoms,lvs)
    print("atomic geometry read")
    return Ratin ;

# procedures for sorting eigenstates:

def to_fermi(eig, fermi, orig_fermi=0.0):
    '''
    Shift the fermi level & shift the eigenenergy to the Fermi-Level
    '''
    fermi = orig_fermi if (fermi == None) else fermi + orig_fermi
    print("The Fermi Level: ", fermi, " eV; in FHI-AIMS is the Fermi automatically 0.")
    eig -= fermi
    return eig;

def cut_eigenenergies(eig):
    '''
    Removes eigenstates (molecular orbitals) that are far from the energy important for scanning
    '''	
    j = 1
    global n_min_ , n_max_
    for i in eig:
        n_min_ = j if (i < cut_min_ ) else n_min_
        n_max_ = j if (i < cut_max_ ) else n_max_
        j += 1
    assert (n_min_ < n_max_), "no orbitals left for dI/dV"
    return eig[n_min_:n_max_];

# procedure for handling the coefficients:

def lower_Allorb(coef):
    '''
    Lowering hoppings for some atoms predefined by user
    '''
    if (lower_atoms_ != []):
        print('lowering atoms hoppings for atoms:', lower_atoms_)
        i_coef = 0;
        for j in lower_atoms_:
            coef[:,j,:] *= lower_coefs_[i_coef]
            i_coef +=1
    return coef;

def lower_Dorb(coef):
    '''
    Lowering hoppings for all d-orbitals (automatically); It has physical reason, but simple rescalling is nasty
    '''
    d_rescale=0.2
    if (Ynum_ > 4) and dOrbRes_: #(orbs=='spd') & and d-orbital-rescalling:
        print("!!! Be aware d-orbs are now rescaled by factor of" ,d_rescale) 
        print("This is due to a faster decay of d-orbs in the original basis sets, but simple rescaling is nasty !!!")
        coef[:,:,4:] *= d_rescale
    coeff = coef.flatten()
    return coeff.reshape((n_max_-n_min_,num_at_*Ynum_));

def remove_coeffs(coeffs):
    '''
    Removing the LCAO coefficients for cutted atoms
    '''
    if (cut_at_ != -1):
        coeffs=np.delete(coeffs,list(range(cut_at_*Ynum_,num_at_*Ynum_)),1)
    return coeffs;
    
def pbc_coef(coeffs):
    '''
    Applying PBC to the LCAO Coefficients
    '''
    if ((pbc_ != (0,0))or(pbc_ != (0.0,0.0))) :
        print("applying pbc")
        coeff =np.repeat(coeffs,int(pbc_[0]*2+1)*int(pbc_[1]*2+1),0).flatten()
        global num_at_;	num_at_ *=int(pbc_[0]*2+1)*int(pbc_[1]*2+1)
        coeffs = coeff.reshape((n_max_-n_min_,num_at_*Ynum_))
    return coeffs;

def	handle_coef(coef):
    '''
    Do all the necessary procedures - rescalling (user & d-orbs), cutting and applying PBC
    '''
    coef = lower_Allorb(coef)
    coeffs = lower_Dorb(coef)
    coeffs = remove_coeffs(coeffs)
    return pbc_coef(coeffs);

# procedures for preparing everything for STM:

def	read_AIMS_all(name = 'KS_eigenvectors.band_1.kpt_1.out', geom='geometry.in', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
    '''
    read_AIMS_all(name = 'KS_eigenvectors.band_1.kpt_1.out', geom='geometry.in', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
    read eigen energies, coffecients (0=Fermi Level) from the 'name' file and geometry  from the 'geom' file.
    orbs - 'sp' read only sp structure of valence orbitals or 'spd' orbitals of the sample
    Fermi - set to zero by AIMS itself
    pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
    imaginary = False (other options for future k-points dependency
    cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
    cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
    lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
    header - newer version of aims gives one aditional line with AIMS-UUID to the output files
    '''
    initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs)
    # obtaining the geometry :
    Ratin, at_num = get_AIMS_geom(geom=geom)
    #print "at_num:",at_num

    # getting eigen-energies:
    filein = open(name )
    skip_header = 2
    for i in range(20):
        tmp=filein.readline().split()
        skip_header += 1
        if (len(tmp)>1):
            if (tmp[1]=='Basis'):
                break
    tmp=filein.readline()
    pre_eig = filein.readline().split()
    filein.close()
    pre_eig=np.delete(pre_eig,[0,1,2],0)
    n_bands = len(pre_eig)
    eig = np.zeros(n_bands)
    for i in range(n_bands):
        eig[i] = float(pre_eig[i])
    del pre_eig, tmp;
    eig = to_fermi(eig, fermi, orig_fermi=0.0)
    eig = cut_eigenenergies(eig)
    print("eigenenergies read")
    # ****** TO BE REMOVED ********
    """
    # old slow python procedure
    # finding position of the LCAO coeficients in the AIMS output file & its phase - sign
    tmp = np.genfromtxt(name,skip_header=skip_header, usecols=(1,2,3,4,5),dtype=None)
    orb_pos=np.zeros((num_at_,Ynum_), dtype=np.int)
    orb_sign=np.zeros((num_at_,Ynum_), dtype=np.int)
    orb_pos += -1
    el = elements.ELEMENTS
    for j in range(num_at_):
        Z = at_num[j];
        per = el[Z][2]
        temp=int((np.mod(per,2)-0.5)*2)	# phase of radial function in long distance for l=0: if n even - +1, if odd - -1
        if (orbs == 'sp'):
            orb_sign[j]=[temp,-1*temp,-1*temp,temp]		# {1, 1, 1, -1};(*Dont change, means - +s, +py +pz -px*) but l=1 has opposite phase than l=0 ==>  sign[s]*{1, -1, -1, 1};
        else: # (orbs == 'spd'):
            orb_sign[j]=[temp,-1*temp,-1*temp,temp,-1*temp,-1*temp,-1*temp,temp,-1*temp]		# {1, 1, 1, -1, 1, 1, 1, 1, -1, 1};(*Dont change, means - +s, +py +pz -px +dxy +dyz +dz2 -dxz +dx2y2)
            # but l=1 has opposite phase than l=0 and l=2 is n-1 - the same phase as l=1 ==>  sign[s]*{1, -1, -1, 1, -1, -1, -1, 1, -1};
    for i in range(len(tmp)):
        for j in range(num_at_):
            Z = at_num[j];
            per = el[Z][2]
            if ((tmp[i][0]==j+1)and(tmp[i][1]=='atomic')):
                if (tmp[i][2]==per):
                    if 	(tmp[i][3]=='s'):
                        orb_pos[j,0]=i
                    elif (tmp[i][3]=='p'):
                        if  (tmp[i][4]==-1):
                            orb_pos[j,1]=i
                        elif (tmp[i][4]==0):
                            orb_pos[j,2]=i
                        elif (tmp[i][4]==1):
                            orb_pos[j,3]=i
                elif ((tmp[i][2]==per-1)and(orbs=='spd')and(per>3)):
                    if (tmp[i][3]=='d'):
                        if   (tmp[i][4]==-2):
                            orb_pos[j,4]=i
                        elif (tmp[i][4]==-1):
                            orb_pos[j,5]=i
                        elif (tmp[i][4]==0):
                            orb_pos[j,6]=i
                        elif (tmp[i][4]==1):
                            orb_pos[j,7]=i
                        elif (tmp[i][4]==2):
                            orb_pos[j,8]=i
    # Reading the coefficients and assigning proper sign, just for wanted eigen-energies
    print "The main reading procedure, it can take some time, numpy reading txt can be slow."
    del tmp; del temp;
    tmp = np.genfromtxt(name,skip_header=skip_header, usecols=tuple(xrange(6, n_bands*2+6, 2))) #tmp = np.genfromtxt(name,skip_header=5)#, usecols=(6,))
    tmp = tmp[:,n_min_:n_max_]
    coef = np.zeros((n_max_-n_min_,num_at_,Ynum_))
    for j in range(num_at_):
        for l in range(Ynum_):
            if (orb_pos[j,l]!=-1):
                coef[:,j,l] = tmp[orb_pos[j,l],:]
                coef[:,j,l] *= orb_sign[j,l]
    del tmp;
    """
    # lowering over atoms and applying PBC
    coef = read_AIMS_coefs(name, at_num ) # Maximal error between C++ and python reading ~1.5E-08; reading procedure ~100x faster
    coeffs = handle_coef(coef)
    print("All coefficients read")
    return eig.copy(), coeffs.copy(), Ratin.copy();

def	read_GPAW_all(name = 'OUTPUT.gpw', fermi = None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[] ):
    '''
    read_GPAW_all(name = 'OUTPUT.gpw', fermi = None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lower_atoms=[], lower_coefs=[]):
    This procedure nead to import ASE and GPAW
    read eigen energies, coffecients, Fermi Level and geometry  from the GPAW  *.gpw file.
    If fermi = None then Fermi comes from the GPAW calculation
    orbs - only 'sp' works 	can read only sp structure of valence orbitals (hydrogens_has to be at the end !!!!)
    pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
    imaginary = False (other options for future k-points dependency
    cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
    cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
    lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
    '''
    initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs)
    # obtaining the geometry :
    from ase import Atoms
    from gpaw import GPAW
    calc = GPAW(name)
    slab = calc.get_atoms()
    Ratin = get_GPAW_geom(geom=slab)

    # getting eigen-energies
    n_bands = calc.get_number_of_bands()
    eig = calc.get_eigenvalues(kpt=0, spin=0, broadcast=True)
    at_num = slab.get_atomic_numbers()
    eig = to_fermi(eig, fermi, orig_fermi=calc.get_fermi_level())
    eig = cut_eigenenergies(eig)
    print("eigen-energies read")
    # obtaining the LCAO coefficients (automatically removed unwanted states - molecular orbitals - and atoms)
    coef = np.zeros((n_max_-n_min_,num_at_,Ynum_))
    if (orbs=='spd'):
        print("!!! WARNING: d-orbitals should be in principle working, but coefficients can be wrong, according to my experiences !!!")
        print("DEBUG: going to crazy procedure, which finds, where the d-orbs starts")
        print("from gpaw.utilities.dos import print_projectors; print_projectors('X')")
        print("this prints you where the d-orb should start")
        from gpaw.setup_data import SetupData
        chem_sym=slab.get_chemical_symbols()
        d_orb=np.zeros((num_at_));
        for i in range(num_at_):
            if at_num[i]>2:
                setup =  SetupData(chem_sym[i],'LDA','paw');l_j = setup.l_j;tmp=l_j[:l_j.index(2)];a=[1,3];oo=0;
                for j in range(len(tmp)):
                    oo +=a[tmp[j]];
                d_orb[i]=oo;
    for i in range(n_min_,n_max_):
        h=0
        for j in range(num_at_):
            ii = i-n_min_
            coef[ii,j,0] = calc.wfs.kpt_u[0].C_nM[i,h]
            if (at_num[j]>2):
                coef[ii,j,1] = calc.wfs.kpt_u[0].C_nM[i,h+1]
                coef[ii,j,2] = calc.wfs.kpt_u[0].C_nM[i,h+2]
                coef[ii,j,3] = calc.wfs.kpt_u[0].C_nM[i,h+3]
                if ((orbs=='spd')and(d_orb[j]>1)):
                    coef[ii,j,4] = calc.wfs.kpt_u[0].C_nM[i,h+d_orb[j]]
                    coef[ii,j,5] = calc.wfs.kpt_u[0].C_nM[i,h+d_orb[j]+1]
                    coef[ii,j,6] = calc.wfs.kpt_u[0].C_nM[i,h+d_orb[j]+2]
                    coef[ii,j,7] = calc.wfs.kpt_u[0].C_nM[i,h+d_orb[j]+3]
                    coef[ii,j,8] = calc.wfs.kpt_u[0].C_nM[i,h+d_orb[j]+4]
            h += calc.wfs.setups[j].nao
    #from gpaw.utilities.dos import print_projectors; print_projectors('Cu')
    #print "DEBUG: Cu coeffs:"
    #for i in range(n_min,n_max):
    #	for j in range(15):
    #		print j, calc.wfs.kpt_u[0].C_nM[i,j]
    #	print "DEBUG: coef[sth,0,:]" , coef[i-n_min,0,:] 
    # lowering tunneling for predefined atoms
    # lowering over atoms and applying PBC
    coeffs = handle_coef(coef)
    print("All coefficients read")
    return eig.copy(), coeffs.copy(), Ratin.copy();

def	read_FIREBALL_all(name = 'phi_' , geom='answer.bas', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[]):
    '''
    read_FIREBALL_all(name = 'phi_' , geom='answer.bas', fermi=None, orbs = 'sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[]):
    This procedure uses only local libraries;
    read coffecients and eigen numbers from Fireball made (iwrtcoefs = -2) files phik_0001_s.dat, phik_0001_py.dat ....
    fermi - If None the Fermi Level from the Fireball calculations (in case of molecule and visualising some molecular orbitals it can be move to their energy by putting there real value)
    orbs = 'sp' read only sp structure of valence orbitals or 'spd' orbitals of the sample
    pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell etc.
    imaginary = False (other options for future k-points dependency
    cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
    cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
    lvs = None no lattice vector (cell); for PBC 3x3 array containing the cell vectors has to be put here
    lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
    note: sometimes oxygens have to have hoppings lowered by 0.5 this is under investigation
    '''
    initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs)
    # obtaining the geometry :
    Ratin = get_FIREBALL_geom(geom=geom, lvs=lvs)

    # getting eigen-energies
    filein = open(name+'s.dat' )
    pre_eig = filein.readline().split()
    filein.close()
    #print "DEBUG: num_at_  and pre_eig[0] " , num_at_ ,pre_eig[0], type(num_at_), type(pre_eig[0][0]), "int(pre_eig[0])", int(pre_eig[0]) , "(num_at_<=int(pre_eig[0][0]))", (num_at_<=int(pre_eig[0]))
    assert (num_at_<=int(pre_eig[0])),"coefficients for lower amount of atoms, that atoms in geometry file";
    n_bands= int(pre_eig[1]);
    eig = np.loadtxt(name+'s.dat',skiprows=1, usecols=(0,))
    assert (len(eig)==n_bands), "number of bands wrongly specified"
    eig = to_fermi(eig, fermi, orig_fermi=float(pre_eig[2]))
    del pre_eig;
    eig = cut_eigenenergies(eig)
    print("eigen-energies read")

    print(" loading the LCAO coefficients")
    coef = np.zeros((n_bands,num_at_,Ynum_))
    if (num_at_ > 1):
        coef[:,:,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,:,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,:,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,:,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        if (orbs =='spd'):
            coef[:,:,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,:,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,:,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,:,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,:,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
    else:
        coef[:,0,0] = np.loadtxt(name+'s.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,0,1] = np.loadtxt(name+'py.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,0,2] = np.loadtxt(name+'pz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        coef[:,0,3] = np.loadtxt(name+'px.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
        if (orbs =='spd'):
            coef[:,0,4] = np.loadtxt(name+'dxy.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,0,5] = np.loadtxt(name+'dyz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,0,6] = np.loadtxt(name+'dz2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,0,7] = np.loadtxt(name+'dxz.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )
            coef[:,0,8] = np.loadtxt(name+'dx2y2.dat',skiprows=1,usecols=tuple(range(1, num_at_*2+1, 2)) )

    # removing states (molecular orbitals) that are not wanted
    coef=coef[n_min_:n_max_,:,:]
    # lowering over atoms and applying PBC
    coeffs = handle_coef(coef)
    print("All coefficients read")
    return eig.copy(), coeffs.copy(), Ratin.copy();


def read_CP2K_all(name, fermi=None, orbs='sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[], spin="closed_shell"):
    '''
    read_CP2K_all(name, fermi=None, orbs='sp', pbc=(1,1), imaginary = False, cut_min=-15.0, cut_max=5.0, cut_at=-1, lvs = None, lower_atoms=[], lower_coefs=[], spin="closed_shell"):
    This procedure uses only local libraries;
    read coffecients and eigen numbers from CP2k made (iwrtcoefs = -2) file name-cartesian-mos-1_0.MOLog
    fermi - If None the Fermi Level from the Fireball calculations (in case of molecule and visualising some molecular orbitals it can be move to their energy by putting there real value)
    orbs = 'sp' read only sp structure of valence orbitals or 'spd' orbitals of the sample
    pbc (1,1) - means 3 times 3 cell around the original, (0,0) cluster, (0.5,0.5) 2x2 cell (not arround original) etc.
    imaginary = False (other options for future k-points dependency
    cut_min = -15.0, cut_max = 5.0 - cut off states(=mol  orbitals) bellow cut_min and above cut_max; energy in eV
    cut_at = -1 .. all atoms; eg. cut_at = 15 --> only first fifteen atoms for the current calculations (mostly the 1st layer is the important one)
    lvs = None no lattice vector (cell); for PBC 3x3 array containing the cell vectors has to be put here
    lower_atotms=[], lower_coefs=[] ... do nothing; lower_atoms=[0,1,2,3], lower_coefs=[0.5,0.5,0.5,0.5] lower coefficients (=hoppings) for the first four atoms by 0.5
    note: sometimes oxygens have to have hoppings lowered by 0.5 this is under investigation
    spin="closed_shell" or "alpha" or "beta" for spin-unrestricted calculations
    '''
    initial_check(orbs=orbs, pbc=pbc, imaginary=imaginary, cut_min=cut_min, cut_max=cut_max, cut_at=cut_at, lower_atoms=lower_atoms, lower_coefs=lower_coefs)
    # read geometry
    Ratin = get_FIREBALL_geom(geom=name+".xyz", lvs=lvs, sl=True) # sl - to always skip the 2nd line from the xyz file
    #import ase.io # this is ommited, now we are using only inner functions "
    #geom = ase.io.read(name+".xyz")
    #Ratin = get_GPAW_geom(geom=geom)
    #at_num = geom.get_atomic_numbers()
    labels, eig, occs, evecs, fermi_energy = read_cp2k_MO_file(name+"-cartesian-mos-1_0.MOLog", spin=spin)
    lumo = np.argmax(occs==0.0)
    homo = lumo -1

    # select relevant MOs
    eig = to_fermi(eig, fermi, orig_fermi=fermi_energy)
    eig = cut_eigenenergies(eig)

    # copy coefficients of relevant MOs
    coef = np.zeros((n_max_-n_min_,num_at_,Ynum_))
    for i in range(n_min_,n_max_):
        ii = i-n_min_
        for j, label in enumerate(labels):
            iatom = int(label[1]) - 1
            #print "DEBUG: j, label, iatom", j, label, iatom
            if (iatom >= num_at_):
                #print "DEBUG iatom>= num_at_", iatom, num_at_
                break
            func = label[3]
            if func.endswith("s"):
                coef[ii,iatom,0] += evecs[j,i]
            # beware: unusual order of directions
            elif func.endswith("py"):
                coef[ii,iatom,1] += evecs[j,i]
            elif func.endswith("pz"):
                coef[ii,iatom,2] += evecs[j,i]
            elif func.endswith("px"):
                coef[ii,iatom,3] += evecs[j,i]
            elif func.endswith("dxy"):
                if (orbs =='spd'):
                    coef[ii,iatom,4] += evecs[j,i]
            elif func.endswith("dyz"):
                if (orbs =='spd'):
                    coef[ii,iatom,5] += evecs[j,i]
            elif func.endswith("dz2"):
                if (orbs =='spd'):
                    coef[ii,iatom,6] += evecs[j,i]
            elif func.endswith("dxz"):
                if (orbs =='spd'):
                    coef[ii,iatom,7] += evecs[j,i]
            elif func.endswith("dx2"):
                if (orbs =='spd'):
                    coef[ii,iatom,8] += evecs[j,i]
            elif func.endswith("dy2"):
                if (orbs =='spd'):
                    coef[ii,iatom,8] -= evecs[j,i]  # coef[,,8] = dx2 - dy2
            else:
                raise Exception

    # lowering tunneling for predefined atoms
    # lowering over atoms and applying PBC
    coeffs = handle_coef(coef)
    print("All coefficients read")
    return eig.copy(), coeffs.copy(), Ratin.copy();

#===============================================================================
def read_cp2k_MO_file(fn, spin):
    '''reads files with basis decomposition of states==Molecular orbitals - it can read closed-shell systems, alpha or beta spin separately, meaning you need to run it twice to get both spins,'''
    print(("Reading CP2K MOs from:"+fn))

    # read all lines into memory
    f = open(fn)
    lines = []
    for l in f.readlines():
        l = l.strip()
        if(len(l)==0): continue
        lines.append(l)
    f.close()

    # detect dimensions
    #parts = lines[-3].split()		# works for files with HOMO-LUMO gap only
    parts=None; tmp_i = None; 
    for i in np.arange(1,5):
        if lines[-1*i][0] == "Fermi" or lines[-1*i][0] == "F":
            parts= lines[-1*i-1].split(); tmp_i = -1*i-1; break;
    assert parts != None, "We did't find the -- Fermi -- line - meaning the end of the file"
    nbasis = int(parts[0])
    natoms = int(parts[1])
    nmos = int(lines[-nbasis-2+tmp_i].split()[-1])
    nlines_per_spin = (nbasis+3) * ((nmos+3)//4) + 0 # 2 originally but at the end of the file, there can be different endings - I try to run everything flexible
    print(("Found %d MOs spanned by %d basis functions centered on %d atoms."%(nmos, nbasis, natoms)))
    del tmp_i;

    # handle spin
    if spin == "alpha":
        first_line = 0
        assert lines[first_line].strip() == "ALPHA MO EIGENVALUES, MO OCCUPATION NUMBERS, AND CARTESIAN MO EIGENVECTORS", ("Wrong_line:", lines[first_line].strip())
    elif spin == "beta":
        for i in  np.arange(0,5):
            first_line = nlines_per_spin + i
            if lines[first_line].strip() == "BETA MO EIGENVALUES, MO OCCUPATION NUMBERS, AND CARTESIAN MO EIGENVECTORS":
                break;
        assert lines[first_line].strip() == "BETA MO EIGENVALUES, MO OCCUPATION NUMBERS, AND CARTESIAN MO EIGENVECTORS"
    elif spin == "closed_shell":
        first_line = 0
        assert lines[first_line].strip() == "MO EIGENVALUES, MO OCCUPATION NUMBERS, AND CARTESIAN MO EIGENVECTORS"
    else:
        raise Exception

    # read fermi energy
    last_line = first_line + nlines_per_spin
    #assert lines[last_line].startswith("HOMO-LUMO gap:")
    fermi_line = -1
    for i in range(-3,6):
        if lines[last_line+i].startswith("Fermi energy:"):
            fermi_line = last_line+i; break;
    assert fermi_line >= 0,("fermi_line number <= 0- we din't find fermi line in the MOLog file:", fermi_line)
    fermi_energy = 27.211385 * float(lines[fermi_line].split()[2])

    # unfold table
    idx = []
    evals = []
    occs = []
    evecs = [list() for i in range(nbasis)]
    labels = [l.split()[:4] for l in lines[first_line+4:first_line+nbasis+4]]
    for i in range((nmos+3)//4): # round up
        a = first_line + i*(nbasis+3) + 1
        idx.extend(lines[a].split())
        evals.extend(lines[a+1].split())
        occs.extend(lines[a+2].split())
        for j in range(nbasis):
            parts = lines[a+3+j].split()
            assert parts[:4] == labels[j]
            evecs[j].extend(parts[4:])

    # convert to numpy arrays
    assert np.all(np.array(idx, int) == np.arange(nmos)+1)
    evals = np.array(evals, float)
    occs = np.array(occs, float)
    evecs = np.array(evecs, float)
    assert evals.shape == (nmos,)
    assert occs.shape == (nmos,)
    assert evecs.shape == (nbasis, nmos)

    # convert hartree to eV
    evals = 27.211385 * evals

    # done
    return labels, evals, occs, evecs, fermi_energy

############## END OF LIBRARY ##################################
