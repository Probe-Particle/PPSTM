#!/usr/bin/python

import os
import sys
import numpy as np
from   scipy.ndimage import uniform_filter
from   ctypes import c_int, c_double, c_char_p
import ctypes

from . import basUtils as bU
from . import elements

from . import cpp_utils

#important constants:

hbar       = 6.58211951440e-16 # [eV.s]
aumass     = 1.66053904020e-27 # [kg] 
eVA2_to_Nm = 16.0217662        # [eV/A^2] / [N/m] 
G2Amp      = 7.7480917346E-05  # rescaling into Amper


# ==============================
# ============================== Pure python functions
# ==============================

LIB_PATH = os.path.dirname( os.path.realpath(__file__) )
print(" ProbeParticle Library DIR = ", LIB_PATH)

def standart_check(orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxy=0.0, dxz=0.0, dyz=0.0, dz2=0.0):
    assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
    #assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
    assert ((s > 0.0)or(px > 0.0)or(py > 0.0)or(pz > 0.0)or(dz2 > 0.0)or(dxy > 0.0)or(dxz > 0.0)or(dyz > 0.0)), "all tip orbitals are zero"
    assert ((s >= 0.0)and(px >= 0.0)and(py >= 0.0)and(pz >= 0.0)and(dz2 >= 0.0)and(dxy >= 0.0)and(dxz >= 0.0)and(dyz >= 0.0)), "you cannot have negative current"
    tip = np.zeros((9))
    tip[0] = s
    tip[1] = py
    tip[2] = pz
    tip[3] = px
    tip[4] = dxy
    tip[5] = dyz
    tip[6] = dz2
    tip[7] = dxz
    orb_t = 4 if (orbs=='sp') else 9
    return tip, orb_t;


def dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxy=0.0, dxz=0.0, dyz=0.0, dz2=0.0):
    '''
    dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0):
    V - voltage = (energy vs. the Fermi Level in eV);
    WF - workfunction (normally ~5 eV gives reasonable results),
    eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
    eig - eigenenergies of sample states (=molecular orbitals)
    R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
    coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
    orbs = 'sp' orbitals of the sample (spd don't work at the moment
    s and/or px and/or py and/or pz orbitals at the PP
    unification of all the predefined dI/dV procedures from C++, you can choose, whatever PP orbital you want
    '''
    tip, orb_t = standart_check(orbs=orbs, s=s, px=px, py=py, pz=pz, dxy=dxy, dxz=dxz, dyz=dyz, dz2=dz2)
    cur = dIdV_sp_sp( V, WF, eta, eig, R, Rat, coes, tip, orb_t )
    return cur;

def dIdV_tilt( V, WF, eta ,eig, R, R0, Rat, coes, orbs='sp', pz=0.0, pxy =0.0, dz2=0.0, dxyz=0.0, len_R=4.0, al=1.0):
    '''
    dIdV_tilt( V, WF, eta ,eig, R, R0, Rat, coes, orbs='sp', pz=1.0, pxy =0.0, len_R=4.0, al=1.0):
    V - voltage = (energy vs. the Fermi Level in eV);
    WF - workfunction (normally ~5 eV gives reasonable results),
    eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
    eig - eigenenergies of sample states (=molecular orbitals)
    R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
    coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
    orbs = 'sp' orbitals of the sample (spd don't work at the moment
    pz and/or pxy tilting orbitals on the PP
    len_R - length between the tip and the PP
    al = 1.0 rescaling of the actuall tilting
    '''
    #assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
    assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
    assert ((pz > 0.0)or(pxy > 0.0)or(dz2 > 0.0)or(dxyz > 0.0)), "all tip orbitals are zero"
    assert ((pz >= 0.0)and(pxy >= 0.0)and(dz2 >= 0.0)and(dxyz >= 0.0)), "you cannot have negative current"
    tip = np.zeros((4))
    tip[0] = pz
    tip[1] = pxy
    tip[2] = dz2
    tip[3] = dxyz
    print("tip coefs:", tip)
    orb_t = 4
    cur = dIdV_sp_sp_tilt( V, WF, eta, eig, R, R0, Rat, coes, tip, len_R, al, orb_t)
    return cur;

def STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxy=0.0, dxz=0.0, dyz=0.0, dz2=0.0, WF_decay=1.0):
    '''
    STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, WF_decay=1.0):
    summing more dI/dV via rectangle integration, be aware Work Function is changing with Voltage!
    '''
    assert (float(V) != 0.0),"you cannot have zero Voltage"
    print("STM simulation via more dI/dV calculations")
    print("Be aware Work Function is changing with Voltage by a factor: ",  WF_decay)
    ii=1;
    for v_ in np.linspace(0.,V,nV):
        print("Start to calculate voltage step %d of %d in total." %(ii, nV))
        ii +=1
        assert (WF-v_*WF_decay> 0.1), "Non-physical Work Function or Voltage, together WF <= 0.1 eV	"
        print("WF for this step is: " , WF-v_*WF_decay, " eV")
        i_ = dIdV( v_, WF-v_*WF_decay, eta ,eig, R, Rat, coes, orbs=orbs, s=s, px =px, py=py, pz=pz, dxy=dxy, dxz=dxz, dyz=dyz, dz2=dz2)
        #print "maximal dI/dV: " , max(i_)
        if (v_ == 0):
            cur = i_
        else:
            cur += i_
    cur *= abs(V)*G2Amp
    print("All dI/dV steps done, current rescalled into Ampers")
    return cur;

def MSTM( Vmin, Vmax, dV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, dxy=0.0, dxz=0.0, dyz=0.0, dz2=0.0, WF_decay=1.0):
    '''
    MSTM( Vmin, Vmax,  dV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, WF_decay=1.0):
    summing more dI/dV via rectangle integration, be aware Work Function is changing with Voltage!
    '''
    #assert (float(V) != 0.0),"you cannot have zero Voltage" # You has to have 0 voltage in this scan
    print("Multiple STM simulation via more dI/dV calculations, store - all STMs and dIdVs")
    print("Be aware Work Function is changing with Voltage by a factor: ",  WF_decay)
    Vmin = Vmin if Vmin <= 0.0 else 0.0;
    Vmax = Vmax if Vmax >= 0.0 else 0.0;
    ii=0;
    v0 = -1
    for v_ in np.arange(Vmin,Vmax+0.001,dV):
        print("Start to calculate voltage step %d of %d in total." %(ii, len(np.arange(Vmin,Vmax+0.001,dV))))
        assert (WF-v_*WF_decay> 0.1), "Non-physical Work Function or Voltage, together WF <= 0.1 eV	"
        print("WF for this step is: " , WF-v_*WF_decay, " eV")
        i_ = dIdV( v_, WF-v_*WF_decay, eta ,eig, R, Rat, coes, orbs=orbs, s=s, px =px, py=py, pz=pz, dxy=dxy, dxz=dxz, dyz=dyz, dz2=dz2)
        dIdVs = np.array([i_]) if v_ == Vmin else np.append(dIdVs, np.array([i_]),axis=0)
        if -0.005 < v_ < 0.005:
            v0 = ii
        ii +=1
    print("dIdVs calculated, voltage 0 has number %d (python logic)" %(v0))
    assert (v0 != -1), "no zero voltage in the scan"
    ii=0;
    print("Now integrating STM current")
    for v_ in np.arange(Vmin,Vmax+0.001,dV):
        print("Sample bias Voltage:" , v_)
        cur = dIdVs[v0]
        if (v_ < -0.005):
            #print "DEBUG: Voltages for STM without zero:", np.arange(ii,v0,1) ;
            for i in np.arange(ii,v0,1):
                cur +=dIdVs[int(i)] 
        if (v_ > 0.005):
            #print "DEBUG: Voltages for STM without zero:", np.arange(v0+1,ii+0.001,1) ;
            for i in np.arange(v0+1,ii+0.001,1):
                cur +=dIdVs[int(i)] 
        cur *= abs(0.001)*G2Amp if -0.005 < v_ < 0.005 else abs(v_)*G2Amp;
        STMs = np.array([cur]) if v_ == Vmin else np.append(STMs, np.array([cur]) , axis=0)
        del cur; ii += 1;
    print("All dI/dV steps done, current rescalled into Ampers")
    #print "DEBUG: size of STMs and dIdVs arrays:", STMs.shape , dIdVs.shape
    return STMs, dIdVs;

def IETS_simple( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.05):
    '''
    IETS_simple( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.05)
    V - voltage = (energy vs. the Fermi Level in eV);
    WF - workfunction (normally ~5 eV gives reasonable results),
    eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
    eig - eigenenergies of sample states (=molecular orbitals)
    R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
    coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
    orbs = 'sp' orbitals of the sample (spd don't work at the moment
    s and/or px and/or py and/or pz orbitals at the PP
    unification of all the predefined dI/dV procedures from C++, you can choose, whatever PP orbital you want
    Amp=0.05 amplitude of vibrations in x and y directions (the same as in PPAFM)
    IETS = (dI/dV)/dx+(dI/dV)/dy
    '''
    tip, orb_t = standart_check(orbs=orbs, s=s, px=px, py=py, pz=pz, dxz=dxz, dyz=dyz, dz2=dz2)
    print("!!! Not full IETS calculations !!!: It is calculating only the electronic part of the full IETS.")
    print("You entered very approximative IETS calculations that consist of dT/dx calc. in different positoins of PP")
    print("Vibration (x,y) Amplitude - for numerical derivation - is:",Amp)
    print("approx IETS = (dI/dV)/dx+(dI/dV)/dy")
    cur1 = IETS_sp_sp( V, WF, eta, eig, R, Rat, coes, tip, Amp, orb_t)
    print("IETS done")
    return cur1;

def IETS_complex( V, WF, eta ,eig, R, eigenEner, eigenVec1, eigenVec2, eigenVec3, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.05, M=16):
    '''
    IETS_complex( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.05)
    V - voltage = (energy vs. the Fermi Level in eV);
    WF - workfunction (normally ~5 eV gives reasonable results),
    eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
    eig - eigenenergies of sample states (=molecular orbitals)
    R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
    coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
    orbs = 'sp' orbitals of the sample (spd don't work at the moment
    s and/or px and/or py and/or pz orbitals at the PP
    Amp=0.05 amplitude of vibrations in x and y directions (the same as in PPAFM)
    M = 16 effectiva mass of the vibrating molecule/atom ; in Atomic Units
    denomin = 1/w_vib^2
    IETS = 1/w_vib1^2 * d(dI/dV)/dvib1+ 1/w_vib2^2 *d(dI/dV)/dvib2 ;
    last part: 1/w_vib3^2 *d(dI/dV)/dvib3 is not connected with Frustrated translation
    '''
    tip, orb_t = standart_check(orbs=orbs, s=s, px=px, py=py, pz=pz, dxz=dxz, dyz=dyz, dz2=dz2)
    print("You entered complex IETS calculations that consist of IETS calculations in different positoins of PP")
    print("Vibration (x,y) Amplitude - for numerical derivation - is:",Amp)
    print("IETS = 1/w_vib1^2 * d(dI/dV)/dvib1+ 1/w_vib2^2 *d(dI/dV)/dvib2")
    sh1 = np.array(R.shape)
    sh2 = np.array(eigenEner.shape)
    #print "dimensions of arrays:", sh1, sh2 # DEBUG only
    assert ((sh1[0]==sh2[0])and(sh1[1]==sh2[1])and(sh1[2]==sh2[2])and(sh1[3]==sh2[3])) , "different shape of R (tip positions ) & eigenEner(gies)"
    tmp = eigenEner.copy()
    #due to negative eigen-energies in special points:
    for i in range(len(eigenEner)):
        for l in [0]:    #range(len(eigvalK[0,0,0])):
            eigenEner[i,:,:,l]=uniform_filter(tmp[i,:,:,l], size=3, mode='nearest')
    del tmp
    # the end
    Evib = hbar * np.sqrt( ( eVA2_to_Nm * eigenEner )/( M * aumass ) )
    denomin = 1/(Evib*Evib) # 1/(Evib[:,:,:0]*Evib[:,:,:0]) + ...
    print("Calculating IETS along the softest & middle-soft vibration")
    iets, iets_stm  = IETScomplex( V, WF, eta, eig, R, eigenVec1*Amp, eigenVec2*Amp, denomin, Rat, coes, tip, orb_t )
    print("IETS done, gives you back denominators - 1/w1^2, 1/w2^2  & full IETS signal")
    #return denomin[:,:,:,0]+denomin[:,:,:,1]+denomin[:,:,:,2],iets;
    return denomin[:,:,:,0]+denomin[:,:,:,1], iets_stm, iets;

def before_C( eig, R, Rat, coes, orb_t):
    NoAt = len(Rat)
    NoOrb = len(eig)
    sh = R.shape
    cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
    Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
    assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
    if (len(coes) != 0):
        #DEBUG:#
        #print "NoAt",NoAt
        #print "NoOrb", NoOrb
        #print "orb_t", orb_t
        #print "len(coes)", len(coes)
        #print "len(coes[0])", len(coes[0])
        assert (NoOrb == len(coes)*len(coes[0])/(orb_t*NoAt)), "Different eigennumbers, than basis"	
    print("We're going to C++")
    return NoAt, NoOrb, Npoints, cur_1d, sh;


# ==============================
# ============================== interface to C++ core 
# ==============================

# ============================== interface to C++ core 

cpp_name='ProbeSTM_spd'
#cpp_utils.compile_lib( cpp_name  )
make_name='MSTM' if sys.platform=='darwin' else 'STM'

try:
    ncpu = int(os.environ['OMP_NUM_THREADS'])
except:
    ncpu = 1;
    print("DEBUG: OMP_NUM_THREADS not defined - serial calculations")

#print "DEBUG: ncpu:", ncpu
make_name_end ='PAR' if ncpu > 1.01 else ''
del ncpu;
make_name += make_name_end
print("DEBUG: make_name", make_name);
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
# void proc_dIdVspdspd(   int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur)
lib.proc_dIdVspdspd.argtypes = [ c_int, c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d, array1d ]
lib.proc_dIdVspdspd.restype  = None
def dIdV_sp_sp( V, WF, eta ,eig, R, Rat, coes, tip_coes, orb_t):
    print("Entering the dI/dV ( sp(d)-sp(d) ) procedure")
    NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes, orb_t)
    lib.proc_dIdVspdspd( orb_t, NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, tip_coes, cur_1d)
    print("We're back in Python")
    return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

# void proc_dIdVspsp_tilt( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double len_R, double al, double* eig, double* R_, double* R0_, double* Rat_, double* coesin, double* tip_coes, double* cur)
lib.proc_dIdVspsp_tilt.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, c_double, c_double, array1d, array4d, array4d, array2d, array2d, array1d, array1d ]
lib.proc_dIdVspsp_tilt.restype  = None
def dIdV_sp_sp_tilt( V, WF, eta ,eig, R, R0, Rat, coes, tip_coes, len_R, al, orb_t):
    print("Entering the dI/dV (sp-sp) procedure with tilting orbitals")
    NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes, orb_t)
    lib.proc_dIdVspsp_tilt( NoAt, NoOrb, Npoints, V, WF, eta, len_R, al, eig, R.copy(), R0.copy(), Rat, coes, tip_coes, cur_1d)
    print("We're back in Python")
    return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

# void proc_dIdVspspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double Amp, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur)
lib.proc_IETSspspd.argtypes = [ c_int, c_int, c_int, c_int, c_double, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d, array1d ]
lib.proc_IETSspspd.restype  = None
def IETS_sp_sp( V, WF, eta ,eig, R, Rat, coes, tip_coes, Amp, orb_t):
    print("Entering the IETS (sp-sp(d)) procedure")
    NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes, orb_t)
    lib.proc_IETSspspd( orb_t, NoAt, NoOrb, Npoints, V, WF, eta, Amp, eig, R.copy(), Rat, coes, tip_coes, cur_1d)
    print("We're back in Python")
    return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

#  void proc_IETScomplex( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* eigenVec1_, double* eigenVec2_, double* denomin_,
#                             double* Rat_, double* coesin, double* tip_coes, double* cur, double* cur2)
lib.proc_IETScomplex.argtypes = [ c_int, c_int, c_int, c_int, c_double, c_double, c_double, array1d, array4d, array4d, array4d, array4d,
                                       array2d, array2d, array1d, array1d, array1d ]
lib.proc_IETScomplex.restype  = None
def IETScomplex( V, WF, eta ,eig, R, eigenVec1, eigenVec2, denomin, Rat, coes, tip_coes, orb_t):
    print("Entering the complex IETS (sp-sp(d)) procedure")
    NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes, orb_t)
    cur_2d = cur_1d.copy()
    lib.proc_IETScomplex( orb_t, NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), eigenVec1.copy(), eigenVec2.copy(), denomin.copy(), Rat, coes, tip_coes, cur_1d, cur_2d)
    print("We're back in Python")
    return cur_1d.reshape((sh[0],sh[1],sh[2])).copy(), cur_2d.reshape((sh[0],sh[1],sh[2])).copy();


############## END OF LIBRARY ##################################
