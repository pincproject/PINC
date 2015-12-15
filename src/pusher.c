/**
 * @file		pusher.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle pusher
 * @date		02.12.15
 *
 * Particle pusher schemes
 */

#include "pinc.h"


/******************************************************************************
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/


/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/


void puMoveNoBnd(Population *pop){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];
		long int pStart = iStart*nDims;
		long int pStop = iStop*nDims;

		for(long int p=pStart;p<pStop;p++){
			pos[p] += vel[p];
		}
	}
}
// void puMovePeriodic();
// void puMoveOpen();

//
// for(int d=0;d<nDims;d++){
// 	const double *posloc = pos;
// 	double *velloc = vel+d;
// 	for(long int i=0;i<nParticles;i++){
// 		*velloc += interpA(E,posloc,nGPointsProd);
// 		posloc += nDims;
// 		velloc += nDims;
// 	}
// 	E += nGPointsProd[nDims];
// }

// void puAccLeapfrog3D0();
// void puAccLeapfrog3D1(Population *pop){
//
// 	int nSpecies = pop->nSpecies;
// 	double *pos = pop->pos;
// 	double *vel = pop->vel;
//
// 	for(int s=0;s<nSpecies;s++){
//
// 		long int iStart = pop->iStart[s];
// 		long int iStop = pop->iStop[s];
//
// 		for(long int i=iStart;i<iStop;i++){
// 			double acc = puInterp3D1();
// 		}
// 	}
//
// }
// void puAccLeapfrog3D2();
// void puAccLeapfrogND0();
// void puAccLeapfrogND1();
// void puAccLeapfrogND2();
//
// void puAccBoris3D0();
// void puAccBoris3D1();
// void puAccBoris3D2();
// void puAccBorisND0();
// void puAccBorisND1();
// void puAccBorisND2();
//
// inline void puDistrib3D0();
// inline void puDistrib3D1();
// inline void puDistrib3D2();
// inline void puDistribND0();
// inline void puDistribND1();
// inline void puDistribND2();

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

// static inline double puInterp3D0();
static inline double puInterp3D1(const double *val, const double *pos, const long int *nGPointsProd){

	// Integer parts of position
	int j = (int) pos[0];
	int k = (int) pos[1];
	int l = (int) pos[2];

	// Decimal (cell-referenced) parts of position and their complement
	double x = pos[0]-j;
	double y = pos[1]-k;
	double z = pos[2]-l;
	double xcomp = 1-x;
	double ycomp = 1-y;
	double zcomp = 1-z;

	// Index of neighbouring nodes
	long int p 		= j + k*nGPointsProd[1] + l*nGPointsProd[2];
	long int pj 	= p + 1;
	long int pk 	= p + nGPointsProd[1];
	long int pjk 	= pk + 1;
	long int pl 	= p + nGPointsProd[2];
	long int pjl 	= pl + 1;
	long int pkl 	= pl + nGPointsProd[1];
	long int pjkl 	= pkl + 1;

	// Linear interpolation
	double result =	 zcomp*( ycomp*(xcomp*val[p   ]+x*val[pj  ])
							+y    *(xcomp*val[pk  ]+x*val[pjk ]) )
					+z    *( ycomp*(xcomp*val[pl  ]+x*val[pjl ])
							+y    *(xcomp*val[pkl ]+x*val[pjkl]) );

	return result;

}
// static inline double puInterp3D2();
// static inline double puInterpND0();
// static inline double puInterpND1(double *pos, const GridQuantity *E, const int *nGPointsProd){
//
// 	long int p = 0;
// 	for(int d=0;d<nDims;d++){
// 		int integer = (int)pos[d];
// 		decimal[d] = pos[d]-integer[d];
// 		complement[d] = 1-decimal[d];
// 		p += nGPointsProd[d] * integer;
// 	}
//
// 	return puInterpND1Inner(val,&nGPointsProd[nDims-1],p,&decimal[nDims-1],&complement[nDims-1]);
//
// }
// static double puInterpND1Inner(const double *val, const long int *mul, long int p, double *decimal, double *complement){
//
// 	double result;
// 	if(*mul==1){
// 		result  = *complement*val[p];		// stay
// 		result += *decimal   *val[p+1];		// incr.
// 	} else {
// 		result  = *complement*inner2(val,mul-1,p     ,decimal-1,complement-1);	// stay
// 		result += *decimal   *inner2(val,mul-1,p+*mul,decimal-1,complement-1);	// incr.
// 	}
// 	return result;
//
// }
// static inline double puInterpND2();
