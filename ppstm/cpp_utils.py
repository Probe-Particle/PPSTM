import logging
import os, sys
import ctypes
import numpy as np

logger = logging.getLogger(__name__)

if "PPSTM_RECOMPILE" in os.environ and os.environ["PPSTM_RECOMPILE"] != "":
    _recompile = True
else:
    _recompile = False


# Shared libraries are called .dll on Windows and .so on Linux/MAC by convention
# Note: Windows were not tested yet
system = sys.platform
if system == "win32":
    _lib_ext = "_lib.dll"
else:
    _lib_ext = "_lib.so"

def work_dir( v__file__ ): 
    return os.path.dirname( os.path.realpath( v__file__ ) )

PACKAGE_PATH = work_dir( __file__ )
CPP_PATH     = os.path.normpath( PACKAGE_PATH + '../../cpp/' )

logger.debug(" PACKAGE_PATH = ", PACKAGE_PATH)
logger.debug(" CPP_PATH     = ", CPP_PATH)

def ctypes_make( make_name: str , lib_name:str ):
    '''
    Compile the given part of C++ code using make_name is _recompile=True or if the library, created from lib_name does not exist.
    E.g. make_name = STM & lib_name ProbeSTM_spd ; This ends up with following command  make STM -> creates ProbeSTM_spd.lib.so
    Returns Ctypes dynamic library object necessary for communication between python and C++
    '''
    lib_path    =   CPP_PATH + "/" + lib_name + _lib_ext
    if _recompile or not os.path.exists(lib_path) : # checks if the libraries exist
        what = 'M' + make_name if system =='darwin' else make_name # different recompilation for mac !
        logger.debug("make command:",what)
        current_directory = os.getcwd()
        os.chdir ( CPP_PATH          )
        os.system( "make "+what       )
        os.chdir ( current_directory )
    return ctypes.CDLL( lib_path )      # load dynamic librady object using ctypes 

# define used numpy array types for interfacing with C++

array1i = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array1d = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
array4d = np.ctypeslib.ndpointer(dtype=np.double, ndim=4, flags='CONTIGUOUS')