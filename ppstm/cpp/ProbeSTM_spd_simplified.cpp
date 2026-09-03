
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "Vec3.cpp"
#include "MatOrb.cpp"

// ================= GLOBAL VARIABLES & CONSTANTS

// global variables & constants are defined in MatOrb

//******************** SINGLE ORBITAL conductance:
 // ================== sqrt(G) - single mol. orb. (sp[d]) vs. sp[d]-tip over atoms
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double Gatomsp(int NoAt, int const_orb, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*FUNC( coes+(const_orb*iat), rev_rr[iat], dr[iat] );
	}
 return f;
 }

 // ================== single point dI/dV calculation sp(d) - sp(d)
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double dIdVspsp_vec( const Vec3d& r, int NoAt, int NoOrb, int const_orb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 Vec3d dr[NoAt];
 double radial[NoAt];
 double rev_rr[NoAt];
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


// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -sp vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVspdspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur){
	set_globals(WF);
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
 //printf("inside a function sp\n");
	if ( const_orb == 4  ){ //sp orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating s orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*dIdVspsp_vec<ssp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating py orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*dIdVspsp_vec<pysp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating pz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*dIdVspsp_vec<pzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating px orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*dIdVspsp_vec<pxsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			//printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
			printf("calculating dxy orb. on a tip sp  -- needs some tets\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[4]*dIdVdxysp_vec(      R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}

		}
		if (tip_coes[5] > 0){
			printf("calculating dyz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[5]*dIdVspsp_vec<dyzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[6] > 0){
			printf("calculating dz2 orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[6]*dIdVspsp_vec<dz2sp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[7] > 0){
			printf("calculating dxz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[7]*dIdVspsp_vec<dxzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[8] > 0){
			printf("What a pitty, I cannot find no formulas for dx2-y2 orb. Shame on the programer!\n");
		}
   }
	else 	if ( const_orb == 9  ){ //spd orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating s orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*dIdVspsp_vec<sspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating py orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*dIdVspsp_vec<pyspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating pz orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*dIdVspsp_vec<pzspd>( R[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating px orb. on a tip spd\n");
			#pragma omp parallel for
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
	//delete dr; delete radial; delete rev_rr;
 }

}
