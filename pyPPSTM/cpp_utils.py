
import os, sys

_recompile = True # global settings for recompilation
lib_ext   ='_lib.so'

def work_dir( v__file__ ): 
    return os.path.dirname( os.path.realpath( v__file__ ) )

PACKAGE_PATH = work_dir( __file__ )
CPP_PATH     = os.path.normpath( PACKAGE_PATH + '../../cpp/' )

print(" PACKAGE_PATH = ", PACKAGE_PATH)
print(" CPP_PATH     = ", CPP_PATH)

def make( what="" ):
    if _recompile:
        what = 'M' + what if sys.platform=='darwin' else what # different recompilation for mac !
        print ("DEBUG: make command:",what)
        current_directory = os.getcwd()
        os.chdir ( CPP_PATH          )
        os.system( "make "+what       )
        os.chdir ( current_directory )

def makeclean( ):
    CWD=os.getcwd()
    os.chdir( CPP_PATH )
    os.system("make clean")
    os.chdir(CWD)

