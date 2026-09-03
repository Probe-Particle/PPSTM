
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


 // ================== single point dIdV calculation sp - sp
 double dIdVdxysp_vec( const Vec3d& r, int NoAt, int NoOrb, int const_orb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 // more testing would be good :-) //
 double f = 0.0;
 double Amp=0.05;
 Vec3d drx1[NoAt];
 Vec3d drx2[NoAt];
 double radialx1[NoAt];
 double radialx2[NoAt];
 double rev_rrx1[NoAt];
 double rev_rrx2[NoAt];
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
	double rrix1 = drx1[iat].norm();
	radialx1[iat] = exp(-(rrix1*decay));  
	rev_rrx1[iat] = 1/rrix1;
	double rrix2 = drx2[iat].norm();
	radialx2[iat] = exp(-(rrix2*decay));  
	rev_rrx2[iat] = 1/rrix2;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*(sqr( ( Gatomsp<pysp>(NoAt,const_orb,drx1,rev_rrx1,radialx1, coesin+(i*NoAt*const_orb) ) - Gatomsp<pysp>(NoAt,const_orb,drx2,rev_rrx2,radialx2, coesin+(i*NoAt*const_orb) ) )/Amp ));

	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation pz tilt. - sp(d)
 double dIdVpzsp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
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
	f += Lor(V,eig[i],eta)*sqr(   (len_R-al*ri.z)/len_R * Gatomsp<pzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                - (al*ri.x)/len_R * Gatomsp<pxsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                - (al*ri.y)/len_R * Gatomsp<pysp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )       );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation pxy tilt. - sp(d)
 double dIdVpxysp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
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
	double a = Gatomsp<pzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double b = Gatomsp<pxsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	double c = Gatomsp<pysp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) );
	f += Lor(V,eig[i],eta)* ( sqr( (al*ri.x)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*ri.x)))/len_R * b )
                             +sqr( (al*ri.y)/len_R * a + (2*len_R-sqrt(sqr(len_R)-sqr(al*ri.y)))/len_R * c ) );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation dz2 tilt. - sp
 double dIdVdz2sp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
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
	f += Lor(V,eig[i],eta)*sqr(   (len_R-al*2*ri.z)/len_R * Gatomsp<dz2sp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                + (al*2*ri.x)/len_R * Gatomsp<dxzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )
                                + (al*2*ri.y)/len_R * Gatomsp<dyzsp>(NoAt,const_orb,dr,rev_rr,radial, coesin+(i*NoAt*const_orb) )  );
	}
 f *= Go*Norm;
 return f;
 }

 // ================== single point dI/dV calculation dxz&dyz tilt. - sp
 double dIdVdxyzsp_vec( const Vec3d& r, const Vec3d& ri, double len_R, double al, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin){
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

 // ================== single point IETS calculation sp - sp
 template < double FUNC(double* coe, double rev_rr, const Vec3d& dR) >
 double IETSComplex_vec( const Vec3d& r, const Vec3d& r1, const Vec3d& r2, const Vec3d& w2, int NoAt, int NoOrb, int const_orb, double V, double eta, double* eig, Vec3d * Ratin, double* coesin, double& g){
 //printf("inside a function \n");
 double f = 0.0;
 double a = 0.0;
 double b = 0.0;
 g = 0.0;
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
 double Amp=2*r1.norm();
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
	drx1[iat]= dri+r1;
	drx2[iat]= dri-r1;
	dry1[iat]= dri+r2;
	dry2[iat]= dri-r2;
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
	a = sqr( ( Gatomsp<FUNC>(NoAt,const_orb,drx1,rev_rrx1,radialx1, coesin+(i*NoAt*const_orb) ) - Gatomsp<FUNC>(NoAt,const_orb,drx2,rev_rrx2,radialx2, coesin+(i*NoAt*const_orb) ) )/Amp );
	b = sqr( ( Gatomsp<FUNC>(NoAt,const_orb,dry1,rev_rry1,radialy1, coesin+(i*NoAt*const_orb) ) - Gatomsp<FUNC>(NoAt,const_orb,dry2,rev_rry2,radialy2, coesin+(i*NoAt*const_orb) ) )/Amp );
	f += Lor(V,eig[i],eta)*( w2.x * a + w2.y * b );
	g += Lor(V,eig[i],eta)*( a + b );
	}
 f *= Go*Norm;
 g *= Go*Norm;
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

 // ================== procedure dIdV -sp vs. sp sample with tilting for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVspsp_tilt( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double len_R, double al, double* eig, double* R_, double* R0_, double* Rat_, double* coesin, double* tip_coes, double* cur){
	set_globals(WF);
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * R0 = ( Vec3d * ) R0_;
	Vec3d * Ratin = (Vec3d *) Rat_;
 //printf("inside a function sp\n");
	int const_orb=4;
	if (tip_coes[0] > 0){
		printf("calculating tilting pz orb. on a tip sp\n");
		#pragma omp parallel for
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[0]*dIdVpzsp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
	}
	if (tip_coes[1] > 0){
		printf("calculating tilting pxy orb. on a tip sp\n");
		#pragma omp parallel for
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[1]*dIdVpxysp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	if (tip_coes[2] > 0){
		printf("calculating tilting d2z orb. on a tip sp\n");
		#pragma omp parallel for
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[2]*dIdVdz2sp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	if (tip_coes[3] > 0){
		printf("calculating tilting dxyz orb. on a tip sp\n");
		#pragma omp parallel for
		for (int s=0; s<Npoints; s++){
		 cur[s] += tip_coes[3]*dIdVdxyzsp_vec( R[s], R[s]-R0[s], len_R, al, NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin);
		}
    }
	//delete dr; delete radial; delete rev_rr;
 }

 // ================== procedure IETS -sp vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_IETSspspd( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double Amp, double* eig, double* R_, double* Rat_, double* coesin, double* tip_coes, double* cur){
	set_globals(WF);
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
 //printf("inside a function sp\n");
	if ( const_orb == 4  ){ //sp orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*IETSspsp_vec<ssp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*IETSspsp_vec<pysp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*IETSspsp_vec<pzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[3]*IETSspsp_vec<pxsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[4] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
		if (tip_coes[5] > 0){
			printf("calculating IETS with dyz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[5]*IETSspsp_vec<dyzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[6] > 0){
			printf("calculating IETS with dz2 orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[6]*IETSspsp_vec<dz2sp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[7] > 0){
			printf("calculating IETS with dxz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[7]*IETSspsp_vec<dxzsp>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[8] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
   }
	else 	if ( const_orb == 9  ){ //spd orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[0]*IETSspsp_vec<sspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[1]*IETSspsp_vec<pyspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){
			 cur[s] += tip_coes[2]*IETSspsp_vec<pzspd>( R[s], NoAt, NoOrb, const_orb, V, eta, Amp, eig, Ratin, coesin);
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip spd\n");
			#pragma omp parallel for
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
	//delete dr; delete radial; delete rev_rr;
 }

 // ================== procedure IETS -sp vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_IETScomplex( int const_orb, int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* eigenVec1_, double* eigenVec2_, double* denomin_,
                             double* Rat_, double* coesin, double* tip_coes, double* cur, double* cur2){
	set_globals(WF);
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	Vec3d * eigenVec1 = ( Vec3d * ) eigenVec1_;
	Vec3d * eigenVec2 = ( Vec3d * ) eigenVec2_;
	Vec3d * denomin   = ( Vec3d * ) denomin_;
 //printf("inside a function sp\n");
	if ( const_orb == 4  ){ //sp orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[0]*IETSComplex_vec<ssp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[0]*gg;
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[1]*IETSComplex_vec<pysp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[1]*gg;
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[2]*IETSComplex_vec<pzsp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[2]*gg;
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[3]*IETSComplex_vec<pxsp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[3]*gg;
			}
		}
		if (tip_coes[4] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
		if (tip_coes[5] > 0){
			printf("calculating IETS with dyz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[5]*IETSComplex_vec<dyzsp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[5]*gg;
			}
		}
		if (tip_coes[6] > 0){
			printf("calculating IETS with dz2 orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[6]*IETSComplex_vec<dz2sp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[6]*gg;
			}
		}
		if (tip_coes[7] > 0){
			printf("calculating IETS with dxz orb. on a tip sp\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[7]*IETSComplex_vec<dxzsp>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[7]*gg;
			}
		}
		if (tip_coes[8] > 0){
			printf("What a pitty, I cannot find no formulas for dxy orb. Shame on the programer!\n");
		}
   }
	else 	if ( const_orb == 9  ){ //spd orbitals of sample
		if (tip_coes[0] > 0){
			printf("calculating IETS with s orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[0]*IETSComplex_vec<sspd>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[0]*gg;
			}
		}
		if (tip_coes[1] > 0){
			printf("calculating IETS with py orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[1]*IETSComplex_vec<pyspd>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[1]*gg;
			}
		}
		if (tip_coes[2] > 0){
			printf("calculating IETS with pz orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[2]*IETSComplex_vec<pzspd>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[2]*gg;
			}
		}
		if (tip_coes[3] > 0){
			printf("calculating IETS with px orb. on a tip spd\n");
			#pragma omp parallel for
			for (int s=0; s<Npoints; s++){ double gg = 0;
			 cur[s] += tip_coes[3]*IETSComplex_vec<pxspd>( R[s], eigenVec1[s], eigenVec2[s], denomin[s], NoAt, NoOrb, const_orb, V, eta, eig, Ratin, coesin, gg);
			 cur2[s] += tip_coes[3]*gg;
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
