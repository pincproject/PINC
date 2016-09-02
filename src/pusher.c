/**
 * @file		pusher.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle pusher
 * @date		02.12.15
 *
 * Particle pusher schemes
 */

#include "core.h"
#include "pusher.h"
#include <math.h>

/******************************************************************************
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/

/** @name Interpolators (used in accelerators)
 * @brief	Interpolates field on grid to position of particle
 * @param[out]		result		Vector value at position
 * @param			pos			Position of particle
 * @param			val			Grid values (e.g. E->val)
 * @param			sizeProd	sizeProd of grid (e.g. E->sizeProd)
 * @param			nDims		Number of dimensions (if not fixed)
 * @param[in,out]	integer		Integer part of particle position
 * @param[in,out]	decimal		Decimal part of particle position
 * @param[in,out]	complement	One minus decimal part of particle position
 * @param			mul			Multiple in order to increment one in a certain direction (from sizeProd, used in recursive algorithm)
 * @param			lastMul		Last mul to use in recursive algorithm (equals sizeProd[1])
 * @param			factor		Multiplicative factor propagated forward in recursion
 * @param			p			Index of lower corner node
 * @return	void
 *
 * Please, just trust me, they work! Go have a look at the functions in
 * pusher.h.
 */
///@{
static inline void puInterp3D1(	double *result, const double *pos,
								const double *val, const long int *sizeProd);
static inline void puInterpND0(	double *result, const double *pos,
								const double *val, const long int *sizeProd,
								int nDims);
static inline void puInterpND1(	double *result, const double *pos,
								const double *val, const long int *sizeProd,
								int nDims, int *integer, double *decimal,
								double *complement);
static void puInterpND1Inner(	double *result, const double *val, long int p,
								const long int *mul, long int lastMul,
								int nDims, double *decimal, double *complement,
								double factor);
static void puDistrND1Inner(	double *val, long int p, const long int *mul,
								long int lastMul, double *decimal,
								double *complement, double factor);
///@}
/**
 * @brief	Adds cross product of a and b to res
 * @param	a		Vector (of length 3)
 * @param	b		Vector (of length 3)
 * @param	res		Resulting vector (of length 3)
 * @return	void
 */
static inline void addCross(const double *a, const double *b, double *res);

/**
 * @brief	Sanity check of accelerator and distributor functions
 * @param	ini		Input file
 * @param	name	Name of function to check for (for use in errors)
 * @param	dim		Dimensionality function works for (0 for all)
 * @param	order	Order of interpolation used
 * @return	void
 *
 * To be used in _set() functions to test validity of ini-file (since this is
 * the same for all accelerator/distributor functions, only depending on
 * order and dimensionality)
 */
static void puSanity(dictionary *ini, const char* name, int dim, int order);

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

funPtr puAcc3D1_set(dictionary *ini){
	puSanity(ini,"puAcc3D1",3,1);
	return puAcc3D1;
}
void puAcc3D1(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			double dv[3];
			puInterp3D1(dv,&pos[p],val,sizeProd);
			for(int d=0;d<nDims;d++) vel[p+d] += dv[d];
		}

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}
}

funPtr puAcc3D1KE_set(dictionary *ini){
	puSanity(ini,"puAcc3D1KE",3,1);
	return puAcc3D1KE;
}
void puAcc3D1KE(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double *kinEnergy = pop->kinEnergy;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		kinEnergy[s]=0;

		for(long int p=pStart;p<pStop;p+=nDims){
			double dv[3];
			puInterp3D1(dv,&pos[p],val,sizeProd);
			double velSquared=0;
			for(int d=0;d<nDims;d++){
				velSquared += vel[p+d]*(vel[p+d]+dv[d]);
				vel[p+d] += dv[d];
			}
			kinEnergy[s]+=velSquared;
		}

		kinEnergy[s]*=mass[s];

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}
}
funPtr puAccND1KE_set(dictionary *ini){
	puSanity(ini,"puAccND1KE",0,1);
	return puAccND1KE;
}
void puAccND1KE(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double *kinEnergy = pop->kinEnergy;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	double *dv = malloc(nDims*sizeof(*dv));
	int *integer = malloc(nDims*sizeof(*integer));
	double *decimal = malloc(nDims*sizeof(*decimal));
	double *complement = malloc(nDims*sizeof(*complement));

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		kinEnergy[s]=0;

		for(long int p=pStart;p<pStop;p+=nDims){

			puInterpND1(dv,&pos[p],val,sizeProd,nDims,integer,decimal,complement);
			double velSquared=0;
			for(int d=0;d<nDims;d++){
				velSquared += vel[p+d]*(vel[p+d]+dv[d]);
				vel[p+d] += dv[d];
			}
			kinEnergy[s]+=velSquared;
		}

		kinEnergy[s]*=mass[s];

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}

	free(dv);
	free(integer);
	free(decimal);
	free(complement);
}

funPtr puAccND1_set(dictionary *ini){
	puSanity(ini,"puAccND1",0,1);
	return puAccND1;
}
void puAccND1(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	double *dv = malloc(nDims*sizeof(*dv));
	int *integer = malloc(nDims*sizeof(*integer));
	double *decimal = malloc(nDims*sizeof(*decimal));
	double *complement = malloc(nDims*sizeof(*complement));

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;


		for(long int p=pStart;p<pStop;p+=nDims){

			puInterpND1(dv,&pos[p],val,sizeProd,nDims,integer,decimal,complement);
			for(int d=0;d<nDims;d++){
				vel[p+d] += dv[d];
			}
		}


		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}

	free(dv);
	free(integer);
	free(decimal);
	free(complement);
}

funPtr puAccND0KE_set(dictionary *ini){
	puSanity(ini,"puAccND0KE",0,0);
	return puAccND0KE;
}
void puAccND0KE(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double *kinEnergy = pop->kinEnergy;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	double *dv = malloc(nDims*sizeof(*dv));

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		kinEnergy[s]=0;

		for(long int p=pStart;p<pStop;p+=nDims){

			puInterpND0(dv,&pos[p],val,sizeProd,nDims);
			double velSquared=0;
			for(int d=0;d<nDims;d++){
				velSquared += vel[p+d]*(vel[p+d]+dv[d]);
				vel[p+d] += dv[d];
			}
			kinEnergy[s]+=velSquared;
		}

		kinEnergy[s]*=mass[s];

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}

	free(dv);
}

funPtr puAccND0_set(dictionary *ini){
	puSanity(ini,"puAccND0",0,0);
	return puAccND0KE;
}
void puAccND0(Population *pop, Grid *E){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	double *dv = malloc(nDims*sizeof(*dv));

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){

			puInterpND0(dv,&pos[p],val,sizeProd,nDims);
			for(int d=0;d<nDims;d++){
				vel[p+d] += dv[d];
			}
		}


		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}

	free(dv);
}


void puBoris3D1(Population *pop, Grid *E, const double *T, const double *S){

	int nSpecies = pop->nSpecies;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			double dv[3], vPrime[3];
			puInterp3D1(dv,&pos[p],val,sizeProd);

			// Add half the acceleration (becomes v minus in B&L notation)
			for(int d=0;d<nDims;d++) vel[p+d] += 0.5*dv[d];

			// Rotate
			memcpy(vPrime,vel,3*sizeof(*vPrime));
			addCross(vel,&T[3*s],vPrime); // vPrime is now v prime
			addCross(vPrime,&S[3*s],vel); // vel is now v plus (B&L)

			// Compute energy here in KE-version

			// Add half the acceleration
			for(int d=0;d<nDims;d++) vel[p+d] += 0.5*dv[d];
		}

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}
}

void puBoris3D1KE(Population *pop, Grid *E, const double *T, const double *S){

	int nSpecies = pop->nSpecies;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;

	double *kinEnergy = pop->kinEnergy;
	double *mass = pop->mass;

	long int *sizeProd = E->sizeProd;
	double *val = E->val;

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		kinEnergy[s]=0;

		for(long int p=pStart;p<pStop;p+=nDims){
			double dv[3], vPrime[3];
			puInterp3D1(dv,&pos[p],val,sizeProd);

			// Add half the acceleration (becomes v minus in B&L notation)
			for(int d=0;d<nDims;d++) vel[p+d] += 0.5*dv[d];

			// Rotate
			memcpy(vPrime,vel,3*sizeof(*vPrime));
			addCross(vel,&T[3*s],vPrime); // vPrime is now v prime
			addCross(vPrime,&S[3*s],vel); // vel is now v plus (B&L)

			// Compute energy
			double velSquared = 0;
			for(int d=0;d<nDims;d++){
				velSquared += pow(vel[p+d],2);
			}
			kinEnergy[s]+=velSquared;

			// Add half the acceleration
			for(int d=0;d<nDims;d++) vel[p+d] += 0.5*dv[d];
		}

		kinEnergy[s]*=mass[s];

		// Specie-specific re-normalization
		gMul(E,pop->renormE[s]);
	}

}

void puGet3DRotationParameters(dictionary *ini, double *T, double *S){

	int nDims = iniGetInt(ini,"grid:nDims");
	int nSpecies = iniGetInt(ini,"grid:nSpecies");
	double *BExt = iniGetDoubleArr(ini,"fields:BExt",nDims);
	double *charge = iniGetDoubleArr(ini,"population:charge",nSpecies);
	double *mass = iniGetDoubleArr(ini,"population:mass",nSpecies);
	double halfTimeStep = 0.5*iniGetDouble(ini,"time:timeStep");

	for(int s=0;s<nSpecies;s++){
		double factor = halfTimeStep*charge[s]/mass[s];
		double denom = 1;
		for(int p=0;p<3;p++){
			T[3*s+p] = factor*BExt[p];
			denom += pow(T[3*s+p],2);
		}
		double mul = 2.0/denom;
		for(int p=0;p<3;p++){
			S[3*s+p] = mul*T[3*s+p];
		}
	}
}


funPtr puDistr3D1_set(dictionary *ini){
	puSanity(ini,"puDistr3D1",3,1);
	return puDistr3D1;
}
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

			// if(pjkl>=sizeProd[4])
			// 	msg(STATUS,"Particle %i at (%f,%f,%f) out-of-bounds, tried to access node %li",i,pos[0],pos[1],pos[2],pjkl);

			val[p] 		+= xcomp*ycomp*zcomp;
			val[pj]		+= x    *ycomp*zcomp;
			val[pk]		+= xcomp*y    *zcomp;
			val[pjk]	+= x    *y    *zcomp;
			val[pl]     += xcomp*ycomp*z    ;
			val[pjl]	+= x    *ycomp*z    ;
			val[pkl]	+= xcomp*y    *z    ;
			val[pjkl]	+= x    *y    *z    ;

		}

		gMul(rho,pop->renormRho[s]);

	}

}

funPtr puDistrND1_set(dictionary *ini){
	puSanity(ini,"puDistrND1",0,1);
	return puDistrND1;
}
void puDistrND1(const Population *pop, Grid *rho){

	gZero(rho);

	int nDims = pop->nDims;
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;

	int nSpecies = pop->nSpecies;

	int *integer = malloc(nDims*sizeof(*integer));
	double *decimal = malloc(nDims*sizeof(*decimal));
	double *complement = malloc(nDims*sizeof(*complement));

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[nDims*i];

			long int p = 0;

			for(int d=0;d<nDims;d++){
				integer[d] = (int) pos[d];
				decimal[d] = pos[d] - integer[d];
				complement[d] = 1 - decimal[d];

				p += integer[d]*sizeProd[d+1];
			}

			puDistrND1Inner(val,p,&sizeProd[nDims],sizeProd[1],&decimal[nDims-1],&complement[nDims-1],1);

		}

		gMul(rho,pop->renormRho[s]);

	}

	free(integer);
	free(decimal);
	free(complement);
}

static void puDistrND1Inner(	double *val, long int p, const long int *mul,
								long int lastMul, double *decimal,
								double *complement, double factor){

	if(*mul==lastMul){
		val[p     ] += *complement*factor;
		val[p+*mul] += *decimal*factor;
	} else {
		puDistrND1Inner(val,p     ,mul-1,lastMul,decimal-1,complement-1,*complement*factor);
		puDistrND1Inner(val,p+*mul,mul-1,lastMul,decimal-1,complement-1,*decimal   *factor);
	}

}

funPtr puDistrND0_set(dictionary *ini){
	puSanity(ini,"puDistrND0",0,0);
	return puDistrND0;
}
void puDistrND0(const Population *pop, Grid *rho){

	gZero(rho);

	int nDims = pop->nDims;
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;

	int nSpecies = pop->nSpecies;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[nDims*i];

			long int p = 0;

			for(int d=0;d<nDims;d++){
				int integer = (int)(pos[d]+0.5);
				p += integer*sizeProd[d+1];
			}
			val[p]++;

		}

		gMul(rho,pop->renormRho[s]);

	}
}

/******************************************************************************
 * MIGRATION FUNCTIONS (TO BE MOVED TO SEPARATE MODULE)
 *****************************************************************************/

// DEPRECATED (ID-technique doesn't work)
void puBndIdMigrants3D(Population *pop, MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	double *pos = pop->pos;
	double *thresholds = mpiInfo->thresholds;
	int neighborhoodCenter = mpiInfo->neighborhoodCenter;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;

	// By using the dummy to hold data we won't lose track of the beginning of
	// the arrays by, say, migrants[neigh]++.
	long int **migrants = mpiInfo->migrantsDummy;
	for(int neigh=0;neigh<nNeighbors;neigh++){
		migrants[neigh] = mpiInfo->migrants[neigh];
	}
	alSetAll(nEmigrants,nSpecies*nNeighbors,0);

	double lx = thresholds[0];
	double ly = thresholds[1];
	double lz = thresholds[2];
	double ux = thresholds[3];
	double uy = thresholds[4];
	double uz = thresholds[5];

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*3;
		long int pStop = pop->iStop[s]*3;

		for(long int p=pStart;p<pStop;p+=3){
			double x = pos[p];
			double y = pos[p+1];
			double z = pos[p+2];
			int nx = - (x<lx) + (x>=ux);
			int ny = - (y<ly) + (y>=uy);
			int nz = - (z<lz) + (z>=uz);
			int n = 13 + nx + 3*ny + 9*nz;

			if(n!=neighborhoodCenter){
				*(migrants[n]++) = p;
				nEmigrants[n*nSpecies+s]++;
			}
		}
	}

}

// DEPRECATED (ID-technique doesn't work)
void puBndIdMigrantsND(Population *pop, MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *thresholds = mpiInfo->thresholds;
	int neighborhoodCenter = mpiInfo->neighborhoodCenter;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;

	// By using the dummy to hold data we won't lose track of the beginning of
	// the arrays by, say, migrants[neigh]++.
	long int **migrants = mpiInfo->migrantsDummy;
	for(int neigh=0;neigh<nNeighbors;neigh++){
		migrants[neigh] = mpiInfo->migrants[neigh];
	}
	alSetAll(nEmigrants,nSpecies*nNeighbors,0);

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			int n = 0;
			for(int d=nDims-1;d>=0;d--){
				n *= 3;
				n += 1 - (pos[p+d]<thresholds[d]) + (pos[p+d]>=thresholds[nDims+d]);
				// A particle at position x will use j=(int)x and j+1 for
				// interpolation. When x is integer and equal to a threshold, it
				// should migrate if on the upper threshold since it may run out
				// of ghost nodes otherwise (unless the user has specified more
				// ghost layers than necessary)
			}
			if(n!=neighborhoodCenter){
				*(migrants[n]++) = p;
				nEmigrants[n*nSpecies+s]++;
			}
		}
	}
}

// Works
// TODO: Add fault-handling in case of too small "emigrants" buffer
funPtr puExtractEmigrants3D_set(const dictionary *ini){
	return puExtractEmigrants3D;
}
void puExtractEmigrants3D(Population *pop, MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *thresholds = mpiInfo->thresholds;
	const int neighborhoodCenter = 13;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;

	// By using the dummy to hold data we won't lose track of the beginning of
	// the arrays when incrementing the pointer
	double **emigrants = mpiInfo->emigrantsDummy;
	for(int ne=0;ne<nNeighbors;ne++){
		emigrants[ne] = mpiInfo->emigrants[ne];
	}
	alSetAll(nEmigrants,nSpecies*nNeighbors,0);

	double lx = thresholds[0];
	double ly = thresholds[1];
	double lz = thresholds[2];
	double ux = thresholds[3];
	double uy = thresholds[4];
	double uz = thresholds[5];

	// adPrint(thresholds,6);

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*3;
		long int pStop = pop->iStop[s]*3;


		for(long int p=pStart;p<pStop;p+=3){
			double x = pos[p];
			double y = pos[p+1];
			double z = pos[p+2];
			int nx = - (x<lx) + (x>=ux);
			int ny = - (y<ly) + (y>=uy);
			int nz = - (z<lz) + (z>=uz);
			int ne = neighborhoodCenter + nx + 3*ny + 9*nz;

			// if(p==371*3)
			// 	msg(STATUS,"x1: %f",x);

			if(ne!=neighborhoodCenter){
				*(emigrants[ne]++) = x;
				*(emigrants[ne]++) = y;
				*(emigrants[ne]++) = z;
				*(emigrants[ne]++) = vel[p];
				*(emigrants[ne]++) = vel[p+1];
				*(emigrants[ne]++) = vel[p+2];
				nEmigrants[ne*nSpecies+s]++;

				pos[p]   = pos[pStop-3];
				pos[p+1] = pos[pStop-2];
				pos[p+2] = pos[pStop-1];
				vel[p]   = vel[pStop-3];
				vel[p+1] = vel[pStop-2];
				vel[p+2] = vel[pStop-1];

				// if(p==371*3)
				// 	msg(STATUS,"x2: %f",pos[p]);

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;


			}
		}
		// msg(STATUS,"pRange: %li-%li, iStop: %li",pStart,pStop,pop->iStop[s]);
	}
}

// Works
// TODO: Add fault-handling in case of too small "emigrants" buffer
funPtr puExtractEmigrantsND_set(const dictionary *ini){
	return puExtractEmigrantsND;
}
void puExtractEmigrantsND(Population *pop, MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *thresholds = mpiInfo->thresholds;
	int neighborhoodCenter = mpiInfo->neighborhoodCenter;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;

	// By using the dummy to hold data we won't lose track of the beginning of
	// the arrays by, say, migrants[neigh]++.
	double **emigrants = mpiInfo->emigrantsDummy;
	for(int ne=0;ne<nNeighbors;ne++){
		emigrants[ne] = mpiInfo->emigrants[ne];
	}
	alSetAll(nEmigrants,nSpecies*nNeighbors,0);

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			int ne = 0;
			for(int d=nDims-1;d>=0;d--){
				ne *= 3;
				ne += 1 - (pos[p+d]<thresholds[d]) + (pos[p+d]>=thresholds[nDims+d]);
				// A particle at position x will use j=(int)x and j+1 for
				// interpolation. When x is integer and equal to a threshold, it
				// should migrate if on the upper threshold since it may run out
				// of ghost nodes otherwise (unless the user has specified more
				// ghost layers than necessary)
			}
			if(ne!=neighborhoodCenter){
				for(int d=0;d<nDims;d++) *(emigrants[ne]++) = pos[p+d];
				for(int d=0;d<nDims;d++) *(emigrants[ne]++) = vel[p+d];
				nEmigrants[ne*nSpecies+s]++;

				for(int d=0;d<nDims;d++) pos[p+d] = pos[pStop-nDims+d];
				for(int d=0;d<nDims;d++) vel[p+d] = vel[pStop-nDims+d];
				pStop -= nDims;
				p -= nDims;
				pop->iStop[s]--;
			}
		}
	}
}

// Works
// TODO: Add fault-handling in case of too small Population struct
static inline void exchangeNMigrants(MpiInfo *mpiInfo){

	int nSpecies = mpiInfo->nSpecies;
	int nNeighbors = mpiInfo->nNeighbors;
	MPI_Request *send = mpiInfo->send;
	MPI_Request *recv = mpiInfo->recv;

	// Send number of emigrants and receive number of immigrants.
	// Order of reception is not necessarily as determined by the for-loop since
	// we're using non-blocking send/receive.
	for(int ne=0;ne<nNeighbors;ne++){
		if(ne!=mpiInfo->neighborhoodCenter){
			int rank = puNeighborToRank(mpiInfo,ne);
			int reciprocal = puNeighborToReciprocal(ne,mpiInfo->nDims);
			long int *nEmigrants  = &mpiInfo->nEmigrants[nSpecies*ne];
			long int *nImmigrants = &mpiInfo->nImmigrants[nSpecies*ne];
			MPI_Isend(nEmigrants ,nSpecies,MPI_LONG,rank,reciprocal,MPI_COMM_WORLD,&send[ne]);
			MPI_Irecv(nImmigrants,nSpecies,MPI_LONG,rank,ne        ,MPI_COMM_WORLD,&recv[ne]);
		}
	}

	MPI_Waitall(nNeighbors,send,MPI_STATUS_IGNORE);
	MPI_Waitall(nNeighbors,recv,MPI_STATUS_IGNORE);

}

// Works
static inline void shiftImmigrants(MpiInfo *mpiInfo, Grid *grid, int ne){

	double *immigrants = mpiInfo->immigrants;
	int nSpecies = mpiInfo->nSpecies;
	long int nImmigrantsTotal = alSum(&mpiInfo->nImmigrants[ne*nSpecies],nSpecies);
	int nDims = mpiInfo->nDims;

	for(int d=0;d<nDims;d++){
		int n = ne%3-1;
		ne /=3;

		double shift = n*grid->trueSize[d+1];
		for(int i=0;i<nImmigrantsTotal;i++){
			immigrants[d+2*nDims*i] += shift;

			// double pos = immigrants[d+2*nDims*i];
			// if(pos>grid->trueSize[d+1])
			// 	msg(ERROR,"particle %i skipped two domains");

		}

	}

}

// Works
static inline void importParticles(Population *pop, double *particles, long int *nParticles, int nSpecies){

	int nDims = pop->nDims;
	long int *iStop = pop->iStop;

	for(int s=0;s<nSpecies;s++){

		double *pos = &pop->pos[nDims*iStop[s]];
		double *vel = &pop->vel[nDims*iStop[s]];

		for(int i=0;i<nParticles[s];i++){
			for(int d=0;d<nDims;d++) *(pos++) = *(particles++);
			for(int d=0;d<nDims;d++) *(vel++) = *(particles++);
		}

		iStop[s] += nParticles[s];
	}

}

// Works
static inline void exchangeMigrants(Population *pop, MpiInfo *mpiInfo, Grid *grid){

	int nSpecies = mpiInfo->nSpecies;
	int nNeighbors = mpiInfo->nNeighbors;
	long int nImmigrantsAlloc = mpiInfo->nImmigrantsAlloc;
	int nDims = mpiInfo->nDims;
	double **emigrants = mpiInfo->emigrants;
	double *immigrants = mpiInfo->immigrants;
	long int *nImmigrants = mpiInfo->nImmigrants;
	MPI_Request *send = mpiInfo->send;

	for(int ne=0;ne<nNeighbors;ne++){
		if(ne!=mpiInfo->neighborhoodCenter){
			int rank = puNeighborToRank(mpiInfo,ne);
			int reciprocal = puNeighborToReciprocal(ne,nDims);
			long int *nEmigrants  = &mpiInfo->nEmigrants[nSpecies*ne];
			long int length = alSum(nEmigrants,nSpecies)*2*nDims;
			MPI_Isend(emigrants[ne],length,MPI_DOUBLE,rank,reciprocal,MPI_COMM_WORLD,&send[ne]);
		}
	}

	// Since "immigrants" is reused for every receive operation MPI_Irecv cannot
	// be used. However, in order to receive and process whichever comes first
	// MPI_ANY_SOURCE is used.
	for(int a=0;a<nNeighbors-1;a++){

		MPI_Status status;
		MPI_Recv(immigrants,nImmigrantsAlloc,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		int ne = status.MPI_TAG;	// Which neighbor it is from equals the tag

		// adPrint(mpiInfo->immigrants,6);
		shiftImmigrants(mpiInfo,grid,ne);
		// adPrint(mpiInfo->immigrants,6);
		importParticles(pop,immigrants,&nImmigrants[ne*nSpecies],nSpecies);

	}

	MPI_Waitall(nNeighbors,send,MPI_STATUS_IGNORE);

}

// Works
void puMigrate(Population *pop, MpiInfo *mpiInfo, Grid *grid){

	exchangeNMigrants(mpiInfo);
	exchangeMigrants(pop,mpiInfo,grid);

}

void puReflect(){



}

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

static void puSanity(dictionary *ini, const char* name, int dim, int order){

	int nDims = iniGetInt(ini,"grid:nDims");
	int *nGhostLayers = iniGetIntArr(ini,"grid:nGhostLayers",2*nDims);
	double *thresholds = iniGetDoubleArr(ini,"grid:thresholds",2*nDims);

	// TBD: can be improved by checking dimensions separately
	int minLayers = aiMin(nGhostLayers,2*nDims);
	double minThreshold = adMin(thresholds,2*nDims);
	double maxThreshold = adMax(thresholds,2*nDims);

	if(nDims!=dim && dim!=0)
		msg(ERROR,"%s only supports grid:nDims=%d",name,dim);

	int reqLayers = 0;
	if(order==0) reqLayers = 0;
	if(order==1) reqLayers = 1;
	if(order==2) reqLayers = 1;

	if(minLayers<1)
		msg(ERROR,"%s requires grid:nGhostLayers >=%d",name,reqLayers);

	double reqMinThreshold = 0;
	if(order==0) reqMinThreshold = -0.5;
	if(order==1) reqMinThreshold = 0;
	if(order==2) reqMinThreshold = 0.5;

	if(minThreshold<reqMinThreshold)
		msg(ERROR,"%s requires grid:thresholds >=%.1f",name,reqMinThreshold);

	if(minThreshold==reqMinThreshold)
		msg(WARNING,"%s is not very well tested for grid:thresholds of exactly %.1f",name,reqMinThreshold);

	double reqMaxThreshold = minLayers-0.5;

	if(maxThreshold>reqMaxThreshold)
		msg(ERROR,"%s requires grid:thresholds <= grid:nGhostLayers - 0.5",name);

	if(minThreshold==reqMaxThreshold)
		msg(WARNING,"%s is not very well tested for grid:thresholds of exactly equal to grid:nGhostLayers - 0.5",name);

	free(nGhostLayers);
	free(thresholds);

}

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
		result[v] =	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
							+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
					+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
							+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

static inline void puInterpND1(	double *result, const double *pos,
								const double *val, const long int *sizeProd,
								int nDims, int *integer, double *decimal,
								double *complement){

	long int p = 0;
	for(int d=0;d<nDims;d++){

		integer[d] = (int)pos[d];
		decimal[d] = pos[d]-integer[d];
		complement[d] = 1-decimal[d];

		int dd = d+1;
		p += sizeProd[dd] * integer[d];

		result[d] = 0;
	}

	puInterpND1Inner(result,val,p,&sizeProd[nDims],sizeProd[1],nDims,&decimal[nDims-1],&complement[nDims-1],1);

}

static void puInterpND1Inner(	double *result, const double *val, long int p,
								const long int *mul, long int lastMul,
								int nDims, double *decimal, double *complement,
								double factor){

	if(*mul==lastMul){
		for(int d=0;d<nDims;d++){
			result[d] += *complement*factor*val[p+d];			// stay
			result[d] += *decimal   *factor*val[p+d+*mul];		// incr.
		}
	} else {
		puInterpND1Inner(result,val,p     ,mul-1,lastMul,nDims,decimal-1,complement-1,*complement*factor);	// stay
		puInterpND1Inner(result,val,p+*mul,mul-1,lastMul,nDims,decimal-1,complement-1,*decimal*factor);		// incr.
	}

}

static inline void puInterpND0(	double *result, const double *pos,
								const double *val, const long int *sizeProd,
								int nDims){

	long int p = 0;
	for(int d=0;d<nDims;d++){
		int integer = (int)(pos[d]+0.5);
		p += sizeProd[d+1] * integer;
	}

	for(int d=0;d<nDims;d++){
		result[d] = val[p+d];
	}

}


int puNeighborToReciprocal(int neighbor, int nDims){

	int reciprocal = 0;

	for(int d=0;d<nDims;d++){
		reciprocal += (2-(neighbor % 3))*pow(3,d);
		neighbor /= 3;
	}

	return reciprocal;

}

int puNeighborToRank(MpiInfo *mpiInfo, int neighbor){

	int rank = 0;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

	for(int d=0;d<mpiInfo->nDims;d++){
		int n = (neighbor % 3) - 1;
		neighbor /= 3;

		// The additiona nSubdomains[d] in paranthesis to ensure positive modulo
		n = (subdomain[d]+n+nSubdomains[d])%nSubdomains[d];

		rank += n*nSubdomainsProd[d];
	}

	return rank;
}

int puRankToNeighbor(MpiInfo *mpiInfo, int rank){

	int neighbor = 0;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	for(int d=0;d<mpiInfo->nDims;d++){
		int n = (rank % nSubdomains[d]);

		// Additional nSubdomains[d] in paranthesis to ensure positive modulo
		n = (n-subdomain[d]+1+nSubdomains[d])%nSubdomains[d];
		rank /= nSubdomains[d];

		neighbor += n*pow(3,d);
	}

	return neighbor;

}

// Adds cross product (Cross product is only defined for 3D, assuming vector 3 long)
static inline void addCross(const double *a, const double *b, double *res){
	res[0] +=  (a[1]*b[2]-a[2]*b[1]);
	res[1] += -(a[0]*b[2]-a[2]*b[0]);
	res[2] +=  (a[0]*b[1]-a[1]*b[0]);
}
