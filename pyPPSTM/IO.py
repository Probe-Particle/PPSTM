#!/usr/bin/python

import os
import numpy as np
from scipy.ndimage import uniform_filter
from   ctypes import c_int, c_double, c_char_p
import ctypes

import basUtils as bU
import elements

import cpp_utils
import ReadSTM

#important constants:

hbar       = 6.58211951440e-16 # [eV.s]
aumass     = 1.66053904020e-27 # [kg] 
eVA2_to_Nm = 16.0217662        # [eV/A^2] / [N/m] 
G2Amp      = 7.7480917346E-05  # rescaling into Amper


# ==============================
# ============================== Pure python functions
# ==============================

LIB_PATH = os.path.dirname( os.path.realpath(__file__) )
print " ProbeParticle Library DIR = ", LIB_PATH


cpp_name='IO'
#cpp_utils.compile_lib( cpp_name  )
cpp_utils.make("IO")
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
# int read_AIMS_coefs(char *fname, double* coefs, int nMO, int nAtom, int nPerAtoms ){
lib.read_AIMS_coefs.argtypes = [ c_char_p, array3d, array1i, c_int, c_int, c_int ]
lib.read_AIMS_coefs.restype  = c_int
def read_AIMS_coefs(fname, at_nums, nAtom, nPerAtom=9 ):
    eigs = ReadSTM.getAimsEigenE(fname)
    nMO = len(eigs)
    periods = np.array([ elements.ELEMENTS[iZ][2] for iZ in at_nums ], dtype=np.int32)
    for iZ,per in zip(at_nums,periods): print iZ,per
    coefs = np.zeros( (nMO,nAtom,nPerAtom) )
    lib.read_AIMS_coefs( fname,coefs, periods, nMO, nAtom, nPerAtom );
    return coefs.copy(), eigs.copy()

