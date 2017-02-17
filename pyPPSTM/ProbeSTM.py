#!/usr/bin/python

import os
import numpy as np
from   ctypes import c_int, c_double, c_char_p
import ctypes
import basUtils as bU
import elements

import cpp_utils


# ==============================
# ============================== Pure python functions
# ==============================

LIB_PATH = os.path.dirname( os.path.realpath(__file__) )
print " ProbeParticle Library DIR = ", LIB_PATH

def standart_check(orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0):
	assert ((orbs == 'sp')or(orbs == 'spd')), "sorry I can't do different orbitals" 
	#assert (orbs == 'sp'), "sorry I can't do different orbitals" 	
	assert ((s > 0.0)or(px > 0.0)or(py > 0.0)or(pz > 0.0)or(dz2 > 0.0)or(dxz > 0.0)or(dyz > 0.0)), "all tip orbitals are zero"
	assert ((s >= 0.0)and(px >= 0.0)and(py >= 0.0)and(pz >= 0.0)and(dz2 >= 0.0)and(dxz >= 0.0)and(dyz >= 0.0)), "you cannot have negative current"
	tip = np.zeros((9))
	tip[0] = s
	tip[1] = py
	tip[2] = pz
	tip[3] = px
	tip[5] = dyz
	tip[6] = dz2
	tip[7] = dxz
	orb_t = 4 if (orbs=='sp') else 9
	return tip, orb_t;


def dIdV( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0):
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
	tip, orb_t = standart_check(orbs=orbs, s=s, px=px, py=py, pz=pz, dxz=dxz, dyz=dyz, dz2=dz2)
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
	print "tip coefs:", tip
	cur = dIdV_sp_sp_tilt( V, WF, eta, eig, R, R0, Rat, coes, tip, len_R, al)
	return cur;

def STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, WF_decay=1.0):
	'''
	STM( V, nV, WF, eta ,eig, R, Rat, coes, orbs='sp', s=1.0, px =0.0, py=0.0, pz=0.0, WF_decay=1.0):
	summing more dI/dV via rectangle integration, be aware Work Function is changing with Voltage!
	'''
	assert (float(V) != 0.0),"you cannot have zero Voltage"
	print "STM simulation via more dI/dV calculations"
	print "Be aware Work Function is changing with Voltage by a factor: ",  WF_decay
	ii=1;
	for v_ in np.linspace(0.,V,nV):
		print "Start to calculate voltage step %d of %d in total." %(ii, nV)
		ii +=1
		assert (WF-v_*WF_decay> 0.1), "Non-physical Work Function or Voltage, together WF <= 0.1 eV	"
		print "WF for this step is: " , WF-v_*WF_decay, " eV"
		i_ = dIdV( v_, WF-v_*WF_decay, eta ,eig, R, Rat, coes, orbs=orbs, s=s, px =px, py=py, pz=pz, dxz=dxz, dyz=dyz, dz2=dz2)
		#print "maximal dI/dV: " , max(i_)
		if (v_ == 0):
			cur = i_
		else:
			cur += i_
	cur *= abs(V)*7.7480917346E-05 # rescaling into Amper
	print "All dI/dV steps done, current rescalled into Ampers"
	return cur;

def IETS_simple( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.02):
	'''
	IETS_simple( V, WF, eta ,eig, R, Rat, coes, orbs='sp', s=0.0, px =0.0, py=0.0, pz=0.0, dxz=0.0, dyz=0.0, dz2=0.0, Amp=0.02)
	V - voltage = (energy vs. the Fermi Level in eV);
	WF - workfunction (normally ~5 eV gives reasonable results),
	eta - energy smearing (0.5-0.30 eV) deppending on system (for single orbital very low number
	eig - eigenenergies of sample states (=molecular orbitals)
	R input of points in whish you calculate dI/dV (relaxed via PP afm, or nonrelaxed via mkSpaceGrid)
	coes -- LCAO coefficients from read_fire_coes (Fireball, maybe FHI-AIMS & mathematica) or read_GPAW_all
	orbs = 'sp' orbitals of the sample (spd don't work at the moment
	s and/or px and/or py and/or pz orbitals at the PP
	unification of all the predefined dI/dV procedures from C++, you can choose, whatever PP orbital you want
	Amp=0.02 amplitude of vibrations in x and y directions
	IETS = (dI/dV)/dx+(dI/dV)/dy
	'''
	tip, orb_t = standart_check(orbs=orbs, s=s, px=px, py=py, pz=pz, dxz=dxz, dyz=dyz, dz2=dz2)
	print "Not working yet"
	print "You entered very simple IETS calculations that consist of IETS calculations in different positoins of PP"
	print "Vibration (x,y) Amplitude is:",Amp
	cur1 = IETS_sp_sp( V, WF, eta, eig, R, Rat, coes, tip, Amp, orb_t)
	print "IETS done"
	return cur1;

def before_C( eig, R, Rat, coes):
	NoAt = len(Rat)
	NoOrb = len(eig)
	sh = R.shape
	cur_1d = np.zeros((sh[0]*sh[1]*sh[2]))
	Npoints = sh[0]*sh[1]*sh[2]	#len(R)/3
	assert (NoOrb == len(coes)), "Different eigennumbers, than basis"
	if (len(coes) != 0):
		assert (NoOrb == len(coes)*len(coes[0])/(4*NoAt)), "Different eigennumbers, than basis"	
	print "We're going to C++"
	return NoAt, NoOrb, Npoints, cur_1d, sh;


# ==============================
# ============================== interface to C++ core 
# ==============================

# ============================== interface to C++ core 

cpp_name='ProbeSTM_spd'
#cpp_utils.compile_lib( cpp_name  )
cpp_utils.make("STM")
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
	print "Entering the dI/dV ( sp(d)-sp(d) ) procedure"
	NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes)
	lib.proc_dIdVspdspd( orb_t, NoAt, NoOrb, Npoints, V, WF, eta, eig, R.copy(), Rat, coes, tip_coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

# void proc_dIdVspsp_tilt( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double len_R, double al, double* eig, double* R_, double* R0_, double* Rat_, double* coesin, double* tip_coes, double* cur)
lib.proc_dIdVspsp_tilt.argtypes = [ c_int, c_int, c_int, c_double, c_double, c_double, c_double, c_double, array1d, array4d, array4d, array2d, array2d, array1d, array1d ]
lib.proc_dIdVspsp_tilt.restype  = None
def dIdV_sp_sp_tilt( V, WF, eta ,eig, R, R0, Rat, coes, tip_coes, len_R, al):
	print "Entering the dI/dV (sp-sp) procedure with tilting orbitals"
	NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes)
	lib.proc_dIdVspsp_tilt( NoAt, NoOrb, Npoints, V, WF, eta, len_R, al, eig, R.copy(), R0.copy(), Rat, coes, tip_coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();

# void proc_dIdVspspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double Amp, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur)
lib.proc_IETSspspd.argtypes = [ c_int, c_int, c_int, c_int, c_double, c_double, c_double, c_double, array1d, array4d, array2d, array2d, array1d, array1d ]
lib.proc_IETSspspd.restype  = None
def IETS_sp_sp( V, WF, eta ,eig, R, Rat, coes, tip_coes, Amp, orb_t):
	print "Entering the IETS (sp-sp(d) procedure"
	NoAt, NoOrb, Npoints, cur_1d, sh = before_C( eig, R, Rat, coes)
	lib.proc_IETSspspd( orb_t, NoAt, NoOrb, Npoints, V, WF, eta, Amp, eig, R.copy(), Rat, coes, tip_coes, cur_1d)
	print "We're back in Python"
	return cur_1d.reshape((sh[0],sh[1],sh[2])).copy();


############## END OF LIBRARY ##################################
