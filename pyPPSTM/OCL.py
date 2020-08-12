#!/usr/bin/python

import os
import pyopencl as cl
import numpy    as np 

# initialize OpenCL

def initCl():
    PACKAGE_PATH = os.path.dirname( os.path.realpath( __file__ ) ); print(PACKAGE_PATH)
    #CL_PATH  = os.path.normpath( PACKAGE_PATH + '../../cl/' )
    CL_PATH  = os.path.normpath( PACKAGE_PATH + '/../cl' )
    #CL_PATH = PACKAGE_PATH+"/cl/"
    print(CL_PATH)
    plats   = cl.get_platforms()
    ctx     = cl.Context(properties=[(cl.context_properties.PLATFORM, plats[0])], devices=None)       
    queue   = cl.CommandQueue(ctx)
    f       = open(CL_PATH+"/STM.cl", 'r')
    fstr    = "".join(f.readlines())
    program = cl.Program(ctx, fstr).build()
    return ctx,queue,program

ctx,queue,program = initCl()

def initArgs(atoms, CAOs, Spectral, rTips ):
    '''
    int nAtoms, int nMOs, 
    __global float4*  atoms,  // [nAtoms]
    __global float4*  CAOs,   // [nMOs*nAtoms] 
    __global float2*  DOSs,   // [nMOs] occupation 
    __global float4*  rTips,  // [global_size]
    __global float *  Iout,   // [global_size] output current
    
   '''
    nDim = rTips.shape
    ntot = (nDim[0]*nDim[1]*nDim[2],)
    print("initArgs rTips ", rTips.shape, ntot)
    nAtoms      = np.int32( len(atoms) )
    nMOs        = np.int32( len(CAOs) )  
    print("initArgs nAtoms, nMOs", nAtoms, nMOs)
    mf          = cl.mem_flags
    cl_Gout     = cl.Buffer(ctx, mf.WRITE_ONLY                   , rTips.nbytes/4     )
    cl_atoms    = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=atoms      )
    cl_CAOs     = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=CAOs       )
    cl_Spectral = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=Spectral   )
    cl_rTips    = cl.Buffer(ctx, mf.READ_ONLY  | mf.COPY_HOST_PTR, hostbuf=rTips      )
    kargs = ( nAtoms, nMOs, cl_atoms, cl_CAOs, cl_Spectral, cl_rTips, cl_Gout )
    return kargs 
        
def run( kargs, nDim, local_size=(32,) ):
    print("run opencl kernel ...")
    global_size = (nDim[0]*nDim[1]*nDim[2],)
    assert ( global_size[0]%local_size[0]==0 ), "number of grid points %i must be divisible by local_group_size %i" %(global_size[0],local_size[0]);
    Gout         = np.zeros( nDim, dtype=np.float32 )
    print("FE.shape",      Gout.shape)
    print("global_size: ", global_size)
    print("local_size:  ", local_size)
    program.Conductance_s_sp( queue, global_size, local_size, *(kargs))
    cl.enqueue_copy    ( queue, Gout, kargs[6] );
    queue.finish()
    print("... opencl kernel DONE")
    return Gout

def getPos(lvec, nDim=None, step=(0.1,0.1,0.1) ):
    if nDim is None:
        nDim = (    int(np.linalg.norm(lvec[3,:])/step[2]),
                    int(np.linalg.norm(lvec[2,:])/step[1]),
                    int(np.linalg.norm(lvec[1,:])/step[0]))
    dCell = np.array( ( lvec[1,:]/nDim[2], lvec[2,:]/nDim[1], lvec[3,:]/nDim[0] ) ) 
    ABC   = np.mgrid[0:nDim[0],0:nDim[1],0:nDim[2]]
    print("nDim",nDim)
    print("ABC[0].shape ", ABC[0].shape)
    X = lvec[0,0] + ABC[2]*dCell[0,0] + ABC[1]*dCell[1,0] + ABC[0]*dCell[2,0]
    Y = lvec[0,1] + ABC[2]*dCell[0,1] + ABC[1]*dCell[1,1] + ABC[0]*dCell[2,1] 
    Z = lvec[0,2] + ABC[2]*dCell[0,2] + ABC[1]*dCell[1,2] + ABC[0]*dCell[2,2] 
    return X, Y, Z
    
def XYZ2float4(X,Y,Z):
    nDim = X.shape
    XYZW = np.zeros( (nDim[0],nDim[1],nDim[2],4), dtype=np.float32)
    XYZW[:,:,:,0] = X
    XYZW[:,:,:,1] = Y
    XYZW[:,:,:,2] = Z
    return XYZW
    
def getPos_f4( lvec, nDim=None, step=(0.1,0.1,0.1) ):
    X,Y,Z = getPos(lvec, nDim=nDim, step=step )
    return  XYZ2float4(X,Y,Z)  
    
def xyzq2float4(xyzs,qs):
    atoms_      = np.zeros( (len(qs),4), dtype=np.float32)
    atoms_[:,:3] = xyzs[:,:]
    atoms_[:, 3] = qs[:]      
    return atoms_
    
def getSpectral( eigenvals, Wf = 1.0, w=0.2 ):
    w2 = w*w
    lorentz     = np.ones( len(eigenvals), dtype=np.float32 )
    lorentz[:] /= ( w**2 + (eigenvals-Wf)**2 )
    return lorentz
    
def CAOsp2f4(CAOs,nAtoms):
    CAOs = CAOs.reshape((-1,nAtoms,4))
    print("CAOs.shape ", CAOs.shape)
    CAOs_      = np.zeros( (len(CAOs),nAtoms,4), dtype=np.float32)
    CAOs_[:,:,0] = CAOs[:,:,0] # s
    CAOs_[:,:,1] = CAOs[:,:,3] # px
    CAOs_[:,:,2] = CAOs[:,:,1] # py  --- because fireball  
    CAOs_[:,:,3] = CAOs[:,:,2] # pz    
    return CAOs_
    
    
    

