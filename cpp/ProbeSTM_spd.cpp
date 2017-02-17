
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Vec3.cpp"

// ================= CONSTANTS

const double Go = 39.47841760435743; //(2*pi)^2
const double eV = 0.036749034;
const double aB = 1.889725989;
const double rev_dpi = 0.1591549431; // 1/(2 Pi)
const double four_pi = 12.566371;    // pi*4
const double N_p = 1.7320508;        // Sqrt(3)
const double N_d = 3.8729833;        // Sqrt(15)
const double N_d2 = 4.4721360;       // Sqrt(5)*0.5
const double I_3 = 0.3333333;        // 1/3

// ================= GLOBAL VARIABLES


static double decay = 1;
//static double decay2 = 1;
static double Norm = 1;
Vec3d * dr; //=new Vec3d[NoAt];
double* radial; //= new double[NoAt];
double* rev_rr; //= new double[NoAt];

// ================= INNER FUNCTIONS

 inline double sqr(double x){
 return x*x;
 }

 inline double trr(double x){
 return x*x*x;
 }

 inline double qur(double x){
 return x*x*x*x;
 }

 inline double Lor(double V, double eig, double eta){
 double f = rev_dpi*eta/((V-eig)*(V-eig)+(0.25*eta*eta));
 return f;
 }

// ***************** Single sp functions for different sp orb.
// ================= orbital conductances (for s tip only now)
                          
 inline double ssp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0];					// s  orb. of sample
 f += coe[1]*dR.y*rev_rr*N_p;		// py orb. of sample
 f += coe[2]*dR.z*rev_rr*N_p;		// pz orb. of sample
 f += coe[3]*dR.x*rev_rr*N_p;		// px orb. of sample
 return f;
 }

 inline double sspd(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0];					                        // s  orb. of sample
 f += coe[1]*dR.y*rev_rr*N_p;		                        // py orb. of sample
 f += coe[2]*dR.z*rev_rr*N_p;		                        // pz orb. of sample
 f += coe[3]*dR.x*rev_rr*N_p;		                        // px orb. of sample
 f += coe[4]*dR.x*dR.y*sqr(rev_rr)*N_d;                     //dxy orb. of sample
 f += coe[5]*dR.y*dR.z*sqr(rev_rr)*N_d;                     //dyz orb. of sample
 f += coe[6]*( 3*sqr(dR.z)*sqr(rev_rr) - 1 )*N_d2;          //dz2 orb. of sample
 f += coe[7]*dR.x*dR.z*sqr(rev_rr)*N_d;                     //dxz orb. of sample
 f += coe[8]*0.5*sqr(rev_rr)*( sqr(dR.x) - sqr(dR.y) )*N_d; //dx2-y2 orb. of sample
 return f;
 }

// ================= orbital conductances (for py tip only now)

 inline double pysp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.y*decay;													// s  orb. of sample
 f += coe[1]*N_p*( -1 + decay*rev_rr*sqr(dR.y) + sqr(rev_rr)*sqr(dR.y) );		// py orb. of sample
 f += coe[2]*N_p*dR.y*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*dR.y*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 return f*rev_rr;
 }

 inline double pyspd(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.y*decay;													// s  orb. of sample
 f += coe[1]*N_p*( -1 + decay*rev_rr*sqr(dR.y) + sqr(rev_rr)*sqr(dR.y) );		// py orb. of sample
 f += coe[2]*N_p*dR.y*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*dR.y*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 f += coe[4]*N_d*dR.x*rev_rr*( 2*sqr(dR.y)*sqr(rev_rr) + decay*sqr(dR.y)*rev_rr - 1 );                           //dxy orb. of sample
 f += coe[5]*N_d*dR.z*rev_rr*( 2*sqr(dR.y)*sqr(rev_rr) + decay*sqr(dR.y)*rev_rr - 1 );                           //dyz orb. of sample
 f += coe[6]*N_d2*dR.y*( 6*sqr(dR.z)*trr(rev_rr) + decay*(3*sqr(dR.z)*sqr(rev_rr)-1) );                          //dz2 orb. of sample
 f += coe[7]*N_d*dR.x*dR.y*dR.z*sqr(rev_rr)*( 2*rev_rr + decay );                                                //dxz orb. of sample
 f += coe[8]*N_d*rev_rr*dR.y*( (sqr(dR.x)-sqr(dR.y))*sqr(rev_rr) + 0.5*decay*(sqr(dR.x)-sqr(dR.y))*rev_rr + 1 ); //dx2-y2 orb. of sample
 return f*rev_rr;
 }

// ================= orbital conductances (for pz tip only now)

 inline double pzsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.z*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.z*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*( -1 + decay*rev_rr*sqr(dR.z) + sqr(rev_rr)*sqr(dR.z) );		// pz orb. of sample
 f += coe[3]*N_p*dR.z*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 return f*rev_rr;
 }

 inline double pzspd(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.z*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.z*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*( -1 + decay*rev_rr*sqr(dR.z) + sqr(rev_rr)*sqr(dR.z) );		// pz orb. of sample
 f += coe[3]*N_p*dR.z*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 f += coe[4]*N_d*dR.x*dR.y*dR.z*sqr(rev_rr)*( 2*rev_rr + decay )   ;                               //dxy orb. of sample
 f += coe[5]*N_d*dR.y*rev_rr*( 2*sqr(dR.z)*sqr(rev_rr) + decay*sqr(dR.z)*rev_rr - 1 );             //dyz orb. of sample
 f += coe[6]*N_d2*( (6*trr(dR.z)*trr(rev_rr)-6*dR.z*rev_rr) + decay*(3*sqr(dR.z)*sqr(rev_rr)-1) ); //dz2 orb. of sample
 f += coe[7]*N_d*dR.x*rev_rr*( 2*sqr(dR.z)*sqr(rev_rr) + decay*sqr(dR.z)*rev_rr - 1 );             //dxz orb. of sample
 f += coe[8]*N_d*sqr(rev_rr)*dR.z*(sqr(dR.x)-sqr(dR.y))*( rev_rr + 0.5*decay  );                   //dx2-y2 orb. of sample
 return f*rev_rr;
 }

// ================= orbital conductances (for px tip only now)

 inline double pxsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.x*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.x*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*dR.x*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*( -1 + decay*rev_rr*sqr(dR.x) + sqr(rev_rr)*sqr(dR.x) );		// px orb. of sample
 return f*rev_rr;
 }

 inline double pxspd(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.x*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.x*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*dR.x*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*( -1 + decay*rev_rr*sqr(dR.x) + sqr(rev_rr)*sqr(dR.x) );		// px orb. of sample
 f += coe[4]*N_d*dR.y*rev_rr*( 2*sqr(dR.x)*sqr(rev_rr) + decay*sqr(dR.x)*rev_rr - 1 );                           //dxy orb. of sample
 f += coe[5]*N_d*dR.x*dR.y*dR.z*sqr(rev_rr)*( 2*rev_rr + decay );                                                //dyz orb. of sample
 f += coe[6]*N_d2*dR.x*( 6*sqr(dR.z)*trr(rev_rr) + decay*(3*sqr(dR.z)*sqr(rev_rr)-1) );                          //dz2 orb. of sample
 f += coe[7]*N_d*dR.z*rev_rr*( 2*sqr(dR.x)*sqr(rev_rr) + decay*sqr(dR.x)*rev_rr - 1 );                           //dxz orb. of sample
 f += coe[8]*N_d*rev_rr*dR.x*( (sqr(dR.x)-sqr(dR.y))*sqr(rev_rr) + 0.5*decay*(sqr(dR.x)-sqr(dR.y))*rev_rr - 1 ); //dx2-y2 orb. of sample
 return f*rev_rr;
 }

// ================= orbital conductances (for dxz tip only now)

 inline double dxzsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.x*dR.y*decay*( rev_rr + decay) ;																			// s  orb. of sample
 f += coe[1]*N_p*dR.x*dR.y*dR.z*rev_rr*( 3*sqr(rev_rr) + 3*decay*rev_rr + sqr(decay) );											// py orb. of sample
 f += coe[2]*N_p*dR.x*( 3*sqr(dR.z)*trr(rev_rr) + 3*decay*sqr(dR.z)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dR.z)*rev_rr - decay );	// pz orb. of sample
 f += coe[3]*N_p*dR.z*( 3*sqr(dR.x)*trr(rev_rr) + 3*decay*sqr(dR.x)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dR.x)*rev_rr - decay );	// px orb. of sample
 return f*sqr(rev_rr);
 }

// ================= orbital conductances (for dyz tip only now)

 inline double dyzsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.y*dR.z*decay*( rev_rr + decay) ;																			// s  orb. of sample
 f += coe[1]*N_p*dR.z*( 3*sqr(dR.y)*trr(rev_rr) + 3*decay*sqr(dR.y)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dR.y)*rev_rr - decay );	// px orb. of sample
 f += coe[2]*N_p*dR.y*( 3*sqr(dR.z)*trr(rev_rr) + 3*decay*sqr(dR.z)*sqr(rev_rr) - rev_rr + sqr(decay)*sqr(dR.z)*rev_rr - decay );	// pz orb. of sample
 f += coe[3]*N_p*dR.x*dR.y*dR.z*rev_rr*( 3*sqr(rev_rr) + 3*decay*rev_rr + sqr(decay) );											// px orb. of sample
 return f*sqr(rev_rr);
 }

// ================= orbital conductances (for dyz tip only now)

 inline double dz2sp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*( -I_3*sqr(decay) + decay*sqr(dR.z)*trr(rev_rr) + sqr(decay)*sqr(dR.z)*sqr(rev_rr) - decay*rev_rr) ;														// s  orb. of sample
 f += coe[1]*N_p*dR.y*rev_rr*( 3*sqr(dR.z)*qur(rev_rr) + 3*decay*sqr(dR.z)*trr(rev_rr) -   sqr(rev_rr) + sqr(decay)*sqr(dR.z)*sqr(rev_rr) -   decay*rev_rr - I_3*sqr(decay) );	// py orb. of sample
 f += coe[2]*N_p*dR.z*rev_rr*( 3*sqr(dR.z)*qur(rev_rr) + 3*decay*sqr(dR.z)*trr(rev_rr) - 3*sqr(rev_rr) + sqr(decay)*sqr(dR.z)*sqr(rev_rr) - 3*decay*rev_rr - I_3*sqr(decay) );	// pz orb. of sample
 f += coe[3]*N_p*dR.x*rev_rr*( 3*sqr(dR.z)*qur(rev_rr) + 3*decay*sqr(dR.z)*trr(rev_rr) -   sqr(rev_rr) + sqr(decay)*sqr(dR.z)*sqr(rev_rr) -   decay*rev_rr - I_3*sqr(decay) );	// px orb. of sample
 return f;
 }

//******************** SINGLE ORBITAL conductance:
 // ================== sqrt(G) - single mol. orb. (sp) vs. sp-tip over atoms
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double Gatomsp(int NoAt, int const_orb, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*FUNC( coes+(const_orb*iat), rev_rr[iat], dr[iat] );
	}
 return f;
 }

 // ================== single point dI/dV calculation sp - sp
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double dIdVspsp_vec( const Vec3d& r, int NoAt, int NoOrb, int const_orb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*sqr( Gatomsp<FUNC>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) ) );
	}
 f *= Go*Norm;
 return f;
 }


 // ================== single point dI/dV calculation pz tilt. - sp
 double dIdVpzsp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*sqr(   (len_R-al*ri.z)/len_R * Gatomsp<pzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                - (al*ri.x)/len_R * Gatomsp<pxsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                - (al*ri.y)/len_R * Gatomsp<pysp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )       );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation pxy tilt. - sp
 double dIdVpxysp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	double a = Gatomsp<pzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double b = Gatomsp<pxsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double c = Gatomsp<pysp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	f += Lor(V,eig[i],eta)* ( sqr( (al*ri.x)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*ri.x)))/len_R * b )
                             +sqr( (al*ri.y)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*ri.y)))/len_R * c ) );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation pxy tilt. - sp
 double dIdVdz2sp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*sqr(   (len_R-al*2*ri.z)/len_R * Gatomsp<dz2sp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                + (al*2*ri.x)/len_R * Gatomsp<dxzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                + (al*2*ri.y)/len_R * Gatomsp<dyzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )  );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation pxy tilt. - sp
 double dIdVdxyzsp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	double a = Gatomsp<dz2sp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double b = Gatomsp<dxzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double c = Gatomsp<dyzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	f += Lor(V,eig[i],eta)* ( sqr( -(al*2*ri.x)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*2*ri.x)))/len_R * b )
                             +sqr( -(al*2*ri.y)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*2*ri.y)))/len_R * c ) );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point IETS calculation sp - sp
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double IETSspsp_vec( const Vec3d& r, int NoAt, int NoOrb, int const_orb, double V, double eta, double Amp, double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 Vec3d drx1[NoAt];
 Vec3d drx2[NoAt];
 Vec3d dry1[NoAt];
 Vec3d dry2[NoAt];
 double radialx1[NoAt];
 double radialx2[NoAt];
 double radialy1[NoAt];
 double radialy2[NoAt];
 double rev_rrx1[NoAt];
 double rev_rrx2[NoAt];
 double rev_rry1[NoAt];
 double rev_rry2[NoAt];
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	drx1[iat]= dri;
	drx1[iat].x +=Amp*0.5;
	drx2[iat]= dri;
	drx2[iat].x +=-Amp*0.5;
	dry1[iat]= dri;
	dry1[iat].y +=Amp*0.5;
	dry2[iat]= dri;
	dry2[iat].y +=-Amp*0.5;
	//double rri = dri.norm();
	//radial[iat] = exp(-(rri*decay));  
	//rev_rr[iat] = 1/rri;
	double rrix1 = drx1[iat].norm();
	radialx1[iat] = exp(-(rrix1*decay));  
	rev_rrx1[iat] = 1/rrix1;
	double rrix2 = drx2[iat].norm();
	radialx2[iat] = exp(-(rrix2*decay));  
	rev_rrx2[iat] = 1/rrix2;
	double rriy1 = dry1[iat].norm();
	radialy1[iat] = exp(-(rriy1*decay));  
	rev_rry1[iat] = 1/rriy1;
	double rriy2 = dry2[iat].norm();
	radialy2[iat] = exp(-(rriy2*decay));  
	rev_rry2[iat] = 1/rriy2;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*(sqr( ( Gatomsp<FUNC>(NoAt,const_orb,drx1,rev_rrx1,radialx1, coesin+(i*NoAt*const_orb) ) - Gatomsp<FUNC>(NoAt,const_orb,drx2,rev_rrx2,radialx2, coesin+(i*NoAt*const_orb) ) )/Amp ) +
					 		sqr( ( Gatomsp<FUNC>(NoAt,const_orb,dry1,rev_rry1,radialy1, coesin+(i*NoAt*const_orb) ) - Gatomsp<FUNC>(NoAt,const_orb,dry2,rev_rry2,radialy2, coesin+(i*NoAt*const_orb) ) )/Amp )   );

	}
 f *= Go*Norm;
 return f;
 }


// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -sp vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVspdspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	if ( const_orb == 4  ){ //sp orbitals of sample
		//int const_orb=4;
		if (tip_coes[0] > 0){
			printf("calculating s orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*dIdVspsp_vec<ssp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating py orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*dIdVspsp_vec<pysp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating pz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*dIdVspsp_vec<pzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating px orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*dIdVspsp_vec<pxsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
		if (tip_coes[5] > 0){
			printf("calculating dyz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[5]*dIdVspsp_vec<dyzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[6] > 0){
			printf("calculating dz2 orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[6]*dIdVspsp_vec<dz2sp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[7] > 0){
			printf("calculating dxz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[7]*dIdVspsp_vec<dxzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[8] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
   }
	else 	if ( const_orb == 9  ){ //spd orbitals of sample
		//int const_orb=9;
		if (tip_coes[0] > 0){
			printf("calculating s orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*dIdVspsp_vec<sspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating py orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*dIdVspsp_vec<pyspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating pz orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*dIdVspsp_vec<pzspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating px orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*dIdVspsp_vec<pxspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[5] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[6] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[7] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[8] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
	}
	delete dr; delete radial; delete rev_rr;
 }

 // ================== procedure dIdV -sp vs. sp sample with tilting for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVspsp_tilt( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double len_R, double al, double* eig, double* R_, double* R0_, double* Rat_, double* coesin, double* tip_coes, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
    Vec3d * R0 = ( Vec3d * ) R0_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	int const_orb=4;
	if (tip_coes[0] > 0){
		printf("calculating tilting pz orb. on a tip sp\n");
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[0]*dIdVpzsp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
	}
	if (tip_coes[1] > 0){
		printf("calculating tilting pxy orb. on a tip sp\n");
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[1]*dIdVpxysp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	if (tip_coes[2] > 0){
		printf("calculating tilting d2z orb. on a tip sp\n");
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[2]*dIdVdz2sp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	if (tip_coes[3] > 0){
		printf("calculating tilting dxyz orb. on a tip sp\n");
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[3]*dIdVdxyzsp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	delete dr; delete radial; delete rev_rr;
 }

 // ================== procedure IETS -sp vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_IETSspspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double Amp, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	if ( const_orb == 4  ){ //sp orbitals of sample
		//int const_orb=4;
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*IETSspsp_vec<ssp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*IETSspsp_vec<pysp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*IETSspsp_vec<pzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*IETSspsp_vec<pxsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
		if (tip_coes[5] > 0){
			printf("calculating IETS with dyz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[5]*IETSspsp_vec<dyzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[6] > 0){
			printf("calculating IETS with dz2 orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[6]*IETSspsp_vec<dz2sp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[7] > 0){
			printf("calculating IETS with dxz orb. on a tip sp\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[7]*IETSspsp_vec<dxzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[8] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
   }
	else 	if ( const_orb == 9  ){ //spd orbitals of sample
		//int const_orb=9;
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*IETSspsp_vec<sspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*IETSspsp_vec<pyspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*IETSspsp_vec<pzspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip spd\n");
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*IETSspsp_vec<pxspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[5] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[6] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[7] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
		if (tip_coes[8] > 0){
			printf("Are you CRAZY, do you know how LONG the formulas would be!!!\n");
		}
	}
	delete dr; delete radial; delete rev_rr;
 }


}
