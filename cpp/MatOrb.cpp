
#ifndef MatOrb_h 
#define MatOrb_h 

//begin

// ================= CONSTANTS

const double Go = 39.47841760435743; //(2*pi)^2
const double eV = 0.036749034;
const double aB = 1.889725989;
const double rev_dpi = 0.1591549431; // 1/(2 Pi)
const double four_pi = 12.566371;    // pi*4
const double N_p = 1.7320508;        // Sqrt(3)
const double N_d = 3.8729833;        // Sqrt(15)
const double N_d2 = 1.1180300;       // Sqrt(5)*0.5
const double I_3 = 0.3333333;        // 1/3

// ================= GLOBAL VARIABLES

static double decay = 1;
//static double decay2 = 1;
static double Norm = 1;

// ================= OTHER FUNCTIONS

 void set_globals(double WF_){
 decay = sqrt( abs(WF_+WF_)*eV);
 Norm = four_pi*decay;
 }

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


//end

#endif
