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

static inline void puInterp3D1(double *result, const double *pos, const double *val, const long int *sizeProd);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

void puMove(Population *pop){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p++){
			pos[p] += vel[p];
		}
	}
}

void puBndPeriodic(Population *pop, const Grid *grid){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	int *size = grid->size;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){
			for(int d=0;d<nDims;d++){
				if(pos[i*nDims+d]>size[d+1])	pos[i*nDims+d] -= size[d+1];
				if(pos[i*nDims+d]<0)			pos[i*nDims+d] += size[d+1];
			}
		}
	}
}
// void puBndPeriodicDD();
// void puBndOpen();
// void puBndOpenDD();

// void puAcc3D0();
void puAcc3D1(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			puInterp3D1(&vel[p],&pos[p],val,sizeProd);
		}

		// Specie-specific re-normalization
		gMulDouble(E,pop->renormE[s]);
	}
}
// void puAcc3D2();
// void puAccND0();
// void puAccND1();
// void puAccND2();
//
// void puAccBoris3D0();
// void puAccBoris3D1();
// void puAccBoris3D2();
// void puAccBorisND0();
// void puAccBorisND1();
// void puAccBorisND2();
//
// inline void puDistr3D0();
void puDistr3D1(const Population *pop, Grid *rho){

	gZero(rho);
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;

	int nSpecies = pop->nSpecies;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[3*i];

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
			long int p 		= j + k*sizeProd[2] + l*sizeProd[3];
			long int pj 	= p + 1; //sizeProd[1];
			long int pk 	= p + sizeProd[2];
			long int pjk 	= pk + 1; //sizeProd[1];
			long int pl 	= p + sizeProd[3];
			long int pjl 	= pl + 1; //sizeProd[1];
			long int pkl 	= pl + sizeProd[2];
			long int pjkl 	= pkl + 1; //sizeProd[1];

			val[p] 		+= xcomp*ycomp*zcomp;
			val[pj]		+= x    *ycomp*zcomp;
			val[pk]		+= xcomp*y    *zcomp;
			val[pjk]	+= x    *y    *zcomp;
			val[pl]     += xcomp*ycomp*z    ;
			val[pjl]	+= x    *ycomp*z    ;
			val[pjkl]	+= x    *y    *z    ;

		}

	}

}
// inline void puDistr3D2();
// inline void puDistrND0();
// inline void puDistrND1();
// inline void puDistrND2();

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

// static inline double puInterp3D0();
static inline void puInterp3D1(double *result, const double *pos, const double *val, const long int *sizeProd){

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
	long int p 		= j*3 + k*sizeProd[2] + l*sizeProd[3];
	long int pj 	= p + 3; //sizeProd[1];
	long int pk 	= p + sizeProd[2];
	long int pjk 	= pk + 3; //sizeProd[1];
	long int pl 	= p + sizeProd[3];
	long int pjl 	= pl + 3; //sizeProd[1];
	long int pkl 	= pl + sizeProd[2];
	long int pjkl 	= pkl + 3; //sizeProd[1];

	// Linear interpolation
	for(int v=0;v<3;v++)
		result[v] +=	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
								+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
						+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
								+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

// static inline double puInterp3D2();
// static inline double puInterpND0();
// static inline double puInterpND1();
// static inline double puInterpND2();
