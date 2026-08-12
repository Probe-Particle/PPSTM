#define GO      39.47841760435743f //(2*pi)^2
#define EV      0.036749034f
#define AB      1.889725989f
#define REV_DPI 0.1591549431f  // 1/(2 Pi)
#define FOUR_PI 12.566371f
#define N_P     1.7320508f  // Sqrt(3)
#define N_D     3.8729833f  // Sqrt(15)
#define N_D2    1.1180300f  // Sqrt(5)*0.5
#define I_3     0.3333333f  // 1/3

inline float sqr(const float x){return x*x;}

inline float trr(const float x){return x*x*x;}

inline float qur(const float x){return x*x*x*x;}

inline float s(constant float* coe, const float rev_rr, const float3 dr, const int const_orb){
    float f = coe[0];					                       // s  orb. of sample
    f += coe[1]*dr.y*rev_rr*N_P;		                       // py orb. of sample
    f += coe[2]*dr.z*rev_rr*N_P;		                       // pz orb. of sample
    f += coe[3]*dr.x*rev_rr*N_P;		                       // px orb. of sample
    if (const_orb == 9) {
        f += coe[4]*dr.x*dr.y*sqr(rev_rr)*N_D;                     //dxy orb. of sample
        f += coe[5]*dr.y*dr.z*sqr(rev_rr)*N_D;                     //dyz orb. of sample
        f += coe[6]*( 3*sqr(dr.z)*sqr(rev_rr) - 1 )*N_D2;          //dz2 orb. of sample
        f += coe[7]*dr.x*dr.z*sqr(rev_rr)*N_D;                     //dxz orb. of sample
        f += coe[8]*0.5f*sqr(rev_rr)*( sqr(dr.x) - sqr(dr.y) )*N_D; //dx2-y2 orb. of sample
    }
    return f;
}

inline float px(constant float* coe, const float rev_rr, const float3 dr, const int const_orb, const float decay){
    float f = coe[0]*dr.x*decay;													// s  orb. of sample
    f += coe[1]*N_P*dr.x*dr.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
    f += coe[2]*N_P*dr.x*dr.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
    f += coe[3]*N_P*( -1 + decay*rev_rr*sqr(dr.x) + sqr(rev_rr)*sqr(dr.x) );		// px orb. of sample
    if (const_orb == 9) {
        f += coe[4]*N_D*dr.y*rev_rr*( 2*sqr(dr.x)*sqr(rev_rr) + decay*sqr(dr.x)*rev_rr - 1 );                           //dxy orb. of sample
        f += coe[5]*N_D*dr.x*dr.y*dr.z*sqr(rev_rr)*( 2*rev_rr + decay );                                                //dyz orb. of sample
        f += coe[6]*N_D2*dr.x*( 6*sqr(dr.z)*trr(rev_rr) + decay*(3*sqr(dr.z)*sqr(rev_rr)-1) );                          //dz2 orb. of sample
        f += coe[7]*N_D*dr.z*rev_rr*( 2*sqr(dr.x)*sqr(rev_rr) + decay*sqr(dr.x)*rev_rr - 1 );                           //dxz orb. of sample
        f += coe[8]*N_D*rev_rr*dr.x*( (sqr(dr.x)-sqr(dr.y))*sqr(rev_rr) + 0.5f*decay*(sqr(dr.x)-sqr(dr.y))*rev_rr - 1 ); //dx2-y2 orb. of sample
    }
    return f*rev_rr;
}

inline float py(constant float* coe, const float rev_rr, const float3 dr, const int const_orb, const float decay){
    float f = coe[0]*dr.y*decay;													// s  orb. of sample
    f += coe[1]*N_P*( -1 + decay*rev_rr*sqr(dr.y) + sqr(rev_rr)*sqr(dr.y) );		// py orb. of sample
    f += coe[2]*N_P*dr.y*dr.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
    f += coe[3]*N_P*dr.y*dr.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
    if (const_orb == 9) {
        f += coe[4]*N_D*dr.x*rev_rr*( 2*sqr(dr.y)*sqr(rev_rr) + decay*sqr(dr.y)*rev_rr - 1 );                           //dxy orb. of sample
        f += coe[5]*N_D*dr.z*rev_rr*( 2*sqr(dr.y)*sqr(rev_rr) + decay*sqr(dr.y)*rev_rr - 1 );                           //dyz orb. of sample
        f += coe[6]*N_D2*dr.y*( 6*sqr(dr.z)*trr(rev_rr) + decay*(3*sqr(dr.z)*sqr(rev_rr)-1) );                          //dz2 orb. of sample
        f += coe[7]*N_D*dr.x*dr.y*dr.z*sqr(rev_rr)*( 2*rev_rr + decay );                                                //dxz orb. of sample
        f += coe[8]*N_D*rev_rr*dr.y*( (sqr(dr.x)-sqr(dr.y))*sqr(rev_rr) + 0.5f*decay*(sqr(dr.x)-sqr(dr.y))*rev_rr + 1 ); //dx2-y2 orb. of sample
    }
    return f*rev_rr;
}

inline float pz(constant float* coe, const float rev_rr, const float3 dr, const int const_orb, const float decay){
    float f = coe[0]*dr.z*decay;													// s  orb. of sample
    f += coe[1]*N_P*dr.z*dr.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
    f += coe[2]*N_P*( -1 + decay*rev_rr*sqr(dr.z) + sqr(rev_rr)*sqr(dr.z) );		// pz orb. of sample
    f += coe[3]*N_P*dr.z*dr.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
    if (const_orb == 9) {
        f += coe[4]*N_D*dr.x*dr.y*dr.z*sqr(rev_rr)*( 2*rev_rr + decay );                                  //dxy orb. of sample
        f += coe[5]*N_D*dr.y*rev_rr*( 2*sqr(dr.z)*sqr(rev_rr) + decay*sqr(dr.z)*rev_rr - 1 );             //dyz orb. of sample
        f += coe[6]*N_D2*( (6*trr(dr.z)*trr(rev_rr)-6*dr.z*rev_rr) + decay*(3*sqr(dr.z)*sqr(rev_rr)-1) ); //dz2 orb. of sample
        f += coe[7]*N_D*dr.x*rev_rr*( 2*sqr(dr.z)*sqr(rev_rr) + decay*sqr(dr.z)*rev_rr - 1 );             //dxz orb. of sample
        f += coe[8]*N_D*sqr(rev_rr)*dr.z*(sqr(dr.x)-sqr(dr.y))*( rev_rr + 0.5f*decay  );                   //dx2-y2 orb. of sample
    }
    return f*rev_rr;
}

inline float dxz_sp(constant float* coe, const float rev_rr, const float3 dr, const float decay){
    float f = coe[0]*dr.x*dr.z*decay*( rev_rr + decay) ;																			 // s  orb. of sample
    f += coe[1]*N_P*dr.x*dr.y*dr.z*rev_rr*( 3*sqr(rev_rr) + 3*decay*rev_rr + sqr(decay) );											 // py orb. of sample
    f += coe[2]*N_P*dr.x*( 3*sqr(dr.z)*trr(rev_rr) + 3*decay*sqr(dr.z)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dr.z)*rev_rr - decay);  // pz orb. of sample
    f += coe[3]*N_P*dr.z*( 3*sqr(dr.x)*trr(rev_rr) + 3*decay*sqr(dr.x)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dr.x)*rev_rr - decay);  // px orb. of sample
    return f*sqr(rev_rr);
}

inline float dyz_sp(constant float* coe, const float rev_rr, const float3 dr, const float decay){
    float f = coe[0]*dr.y*dr.z*decay*( rev_rr + decay);																			     // s  orb. of sample
    f += coe[1]*N_P*dr.z*( 3*sqr(dr.y)*trr(rev_rr) + 3*decay*sqr(dr.y)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dr.y)*rev_rr - decay);  // px orb. of sample
    f += coe[2]*N_P*dr.y*( 3*sqr(dr.z)*trr(rev_rr) + 3*decay*sqr(dr.z)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dr.z)*rev_rr - decay);  // pz orb. of sample
    f += coe[3]*N_P*dr.x*dr.y*dr.z*rev_rr*( 3*sqr(rev_rr) + 3*decay*rev_rr + sqr(decay));											 // px orb. of sample
    return f*sqr(rev_rr);
}

inline float dz2_sp(constant float* coe, const float rev_rr, const float3 dr, const float decay){
    float f = coe[0]*( -I_3*sqr(decay) + decay*sqr(dr.z)*trr(rev_rr) + sqr(decay)*sqr(dr.z)*sqr(rev_rr) - decay*rev_rr);  // s  orb. of sample
    f += coe[1]*N_P*dr.y*rev_rr*( 3*sqr(dr.z)*qur(rev_rr) + 3*decay*sqr(dr.z)*trr(rev_rr) -   sqr(rev_rr) + sqr(decay)*sqr(dr.z)*sqr(rev_rr) -   decay*rev_rr - I_3*sqr(decay) );	// py orb. of sample
    f += coe[2]*N_P*dr.z*rev_rr*( 3*sqr(dr.z)*qur(rev_rr) + 3*decay*sqr(dr.z)*trr(rev_rr) - 3*sqr(rev_rr) + sqr(decay)*sqr(dr.z)*sqr(rev_rr) - 3*decay*rev_rr - I_3*sqr(decay) );	// pz orb. of sample
    f += coe[3]*N_P*dr.x*rev_rr*( 3*sqr(dr.z)*qur(rev_rr) + 3*decay*sqr(dr.z)*trr(rev_rr) -   sqr(rev_rr) + sqr(decay)*sqr(dr.z)*sqr(rev_rr) -   decay*rev_rr - I_3*sqr(decay) );	// px orb. of sample
    return f;
}

inline float lor(const float V, const float eig, const float eta){
    return REV_DPI*eta/((V-eig)*(V-eig)+(0.25f*eta*eta));
}

float didv_spsp_vec(
    const float3 r,
    const int no_at,
    const int no_orb,
    const int const_orb,
    const float v,
    const float wf,
    const float eta,
    constant float * eig,
    constant float * ratin,
    constant float * coesin,
    const uint tip_orb
    ){
    const float decay = sqrt( fabs(2*wf)*EV);
    const float norm = FOUR_PI*decay;
    float f = 0.0f;
    for (int iorb=0; iorb<no_orb; iorb++){
        float gatomsp = 0.0f;
        for(int iat=0; iat<no_at; iat++){
            const float3 dri = (r - vload3(iat, ratin)) * AB;
            const float rri = sqrt(dot(dri,dri));
            const float rev_rr = 1/rri;
            const float radial = exp(-(rri*decay));
            constant float* coe = coesin+(iorb*no_at*const_orb)+(const_orb*iat);
            float g_eig_atom_unscaled = 0.0f;
            if (tip_orb == 0) {
                g_eig_atom_unscaled = s(coe, rev_rr, dri, const_orb);
            } else if (tip_orb == 1) {
                g_eig_atom_unscaled = py(coe, rev_rr, dri, const_orb, decay);
            } else if (tip_orb == 2) {
                g_eig_atom_unscaled = pz(coe, rev_rr, dri, const_orb, decay);
            } else if (tip_orb == 3) {
                g_eig_atom_unscaled = px(coe, rev_rr, dri, const_orb, decay);
            } else if (tip_orb == 5) {
                g_eig_atom_unscaled = dyz_sp(coe, rev_rr, dri, decay);
            } else if (tip_orb == 6) {
                g_eig_atom_unscaled = dz2_sp(coe, rev_rr, dri, decay);
            } else if (tip_orb == 7) {
                g_eig_atom_unscaled = dxz_sp(coe, rev_rr, dri, decay);
            } else {
                printf("tip_orb=%.d is not currently supported!\n", tip_orb);
                return 0.f;
            }
            gatomsp += radial * g_eig_atom_unscaled;
        }
        f += lor(v,eig[iorb],eta)*sqr(gatomsp);
    }
    return f*GO*norm;
}

void _proc_didv_spd_spd(
    const int s,
    const int const_orb,
    const int no_at,
    const int no_orb,
    const int n_points,
    const float v,
    const float wf,
    const float eta,
    constant float* eig,
    constant float* r,
    constant float* ratin,
    constant float* coesin,
    constant float* tip_coes,
    global float* cur
    ){
    cur[s] = 0.f;
    for (int i = 0; i < 4; i++){
        if (tip_coes[i] > 0.0f) {
            cur[s] += tip_coes[i]*didv_spsp_vec( vload3(s, r), no_at, no_orb, const_orb, v, wf, eta, eig, ratin, coesin, i);
        }
    }
    if (tip_coes[4] > 0.0f) {
        printf("tip_coes[4]=%.1f is not currently supported!\n", tip_coes[4]);
        return;
    }
    for (int i = 5; i < 8; i++){
        if (tip_coes[i] > 0.0f) {
            if (const_orb == 4){  // sp orbitals of sample
                cur[s] += tip_coes[i]*didv_spsp_vec( vload3(s, r), no_at, no_orb, const_orb, v, wf, eta, eig, ratin, coesin, i);
            } else {
                printf("tip_coes[%d]=%.1f is not currently supported for const_orb=%d!\n", i, tip_coes[i], const_orb);
                return;
            }
        }
    }
    if (tip_coes[8] > 0.0f) {
        printf("tip_coes[8]=%.1f is not currently supported!\n", tip_coes[8]);
        return;
    }
}

kernel void proc_didv_spd_spd(
    const int const_orb,
    const int no_at,
    const int no_orb,
    const int n_points,
    const float v,
    const float wf,
    const float eta,
    constant float* eig,
    constant float* r,
    constant float* ratin,
    constant float* coesin,
    constant float* tip_coes,
    global float* cur
    ){
    int s = get_global_id(0);
    _proc_didv_spd_spd(s, const_orb, no_at, no_orb, n_points, v, wf, eta, eig, r, ratin, coesin, tip_coes, cur);
}

kernel void proc_didv_spd_spd_sequential(
    const int const_orb,
    const int no_at,
    const int no_orb,
    const int n_points,
    const float v,
    const float wf,
    const float eta,
    constant float* eig,
    constant float* r,
    constant float* ratin,
    constant float* coesin,
    constant float* tip_coes,
    global float* cur
    ){
    for (int s=0; s<n_points; s++){
        _proc_didv_spd_spd(s, const_orb, no_at, no_orb, n_points, v, wf, eta, eig, r, ratin, coesin, tip_coes, cur);
    }
}