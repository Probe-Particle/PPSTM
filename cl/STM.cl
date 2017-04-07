
#define SQRT3    1.7320508f
#define R_SAFE   1e-4

float4 chen_s_sp( float3 dR, float invR ){  // NOTE :  s-Tip vs. sp sample
    return (float4)( 1.0f, dR*(SQRT3*invR) );
}

__kernel void Conductance_s_sp(
    int nAtoms, int nMOs, 
    __global float4*  atoms,    // [nAtoms]
    __global float4*  CAOs,     // [nMOs*nAtoms] 
    __global float *  Spectral, // [nMOs] occupation 
    __global float4*  rTips,    // [global_size]
    __global float *  Gout      // [global_size] output current
){
    __local float4 LATOMs[32];
    __local float4 LCs   [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 tipPos = rTips[iG].xyz;
    
    float G = 0;
    for (int iMO=0; iMO<nMOs; iMO++ ){
        int    iMO0 = iMO*nAtoms;
        float  Ti   = 0;
        for (int iAtom0=0; iAtom0<nAtoms; iAtom0+= nL ){
            int iAtom = iAtom0 + iL;   if(iAtom>=nAtoms) break;
            LATOMs[iL] = atoms[ iAtom        ];
            LCs   [iL] =  CAOs[ iAtom + iMO0 ];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int j=0; j<nL; j++){  
                //  it would be nice to precompute Angular and Radial for all MOs 
                //  but we probably do not have enough local memory
                //  but we cannot switch the loops (iMO,iAtom) since hopping of MOs should be combined incoherently
                float4 atom     = LATOMs[j];
                float3 dR       = tipPos - atom.xyz;
                float  r        = sqrt( dot(dR,dR) );
                float  Radial   = exp(atom.w*r);
                float4 Angular  = chen_s_sp( dR, 1.0f/(r+R_SAFE) );
                Ti             += Radial*dot( LCs[j], Angular );  
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        G += Spectral[iMO] * Ti*Ti;
    }
    Gout[iG] = G;
}



