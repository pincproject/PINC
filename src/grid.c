/**
 * @file		grid.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Grid-struct handling.
 * @date		30.10.15
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#include "core.h"
#include <mpi.h>
#include <math.h>
#include <hdf5.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



/******************************************************************************
 * DECLARAING LOCAL FUNCTIONS
 *****************************************************************************/
/**
 * @brief Returns the ND-index of this MPI node in the global reference frame
 * @param	ini		input settings
 * @return	The N-dimensional index of this MPI node
 */
static int *getSubdomain(const dictionary *ini);



/**
 * @brief Gets, sends, recieves and sets a slice, using MPI
 * @param nSlicePoints		Length of the slice array
 * @param offsetTake 		Offset of slice to get
 * @param offsetSet 		Offset of where to set slice
 * @param d					Dimension
 * @param sendTo 			mpiRank of subdomain to recieve slice
 * @param recvFrom 			mpiRank of subdomain,
							this node is to recieve from
 * @param mpiRank 			mpiRank of this node
 * @param *gridQuantity
 *
 * This exchanges the slices between two predetermined subdomains, belonging
 * to different computernodes. It uses the getSlice and setSlice functions to
 * extract and set the slices, while the communication is done with MPI_Send
 * and MPI_Recv.
 *
 * @see getSlice
 * @see setSlice
 * @see gHaloOpDim
 */

static double gPotEnergyInner(	const double **rhoVal, const double **phiVal, const int *nGhostLayersBefore,
								const int *nGhostLayersAfter, const int *trueSize, const long int *sizeProd);

static double gNeutralizeGridInner(	const double **val, const int *nGhostLayersBefore, const int *nGhostLayersAfter,
									const int *trueSize, const long int *sizeProd);


/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

static double *getSliceInner(double *nextGhost, const double **valp, const long int *mul,
											const int *points, const long int finalMul){

	if(*mul==finalMul){
		for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
		*valp += (*mul)*(*points-1);
	} else {
		for(int j=0; j<*points;j++)
			nextGhost = getSliceInner(nextGhost, valp, mul-1,points-1,finalMul);
	}
	return nextGhost;
}

void getSlice(double *slice, const Grid *grid, int d, int offset){

	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;

	const double *val = grid->val;

	val += offset*sizeProd[d];

	getSliceInner(slice, &val, &sizeProd[rank-1], &size[rank-1], sizeProd[d]);
}

static const double *setSliceInner(const double *nextGhost, double **valp, const long int *mul,
	const int *points, const long int finalMul){

		if(*mul==finalMul){
			for(int j=0;j<*mul;j++) *((*valp)++) = *(nextGhost++);
			*valp += (*mul)*(*points-1);
		} else {
			for(int j=0; j<*points;j++)
				nextGhost = setSliceInner(nextGhost, valp, mul-1,points-1,finalMul);
		}
		return nextGhost;

}

void setSlice(const double *slice, Grid *grid, int d, int offset){

	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;

	double *val = grid->val;

	val += offset*sizeProd[d];
	setSliceInner(slice, &val, &sizeProd[rank-1], &size[rank-1], sizeProd[d]);
}

static const double *addSliceInner(const double *nextGhost, double **valp, const long int *mul,
	const int *points, const long int finalMul){

		if(*mul==finalMul){
			for(int j=0;j<*mul;j++) *((*valp)++) += *(nextGhost++);
			*valp += (*mul)*(*points-1);
		} else {
			for(int j=0; j<*points;j++)
				nextGhost = addSliceInner(nextGhost, valp, mul-1,points-1,finalMul);
		}
		return nextGhost;

}

void addSlice(const double *slice, Grid *grid, int d, int offset){

	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;

	double *val = grid->val;

	val += offset*sizeProd[d];
	addSliceInner(slice, &val, &sizeProd[rank-1], &size[rank-1], sizeProd[d]);
}

static int *getSubdomain(const dictionary *ini){

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Get ini info
	int nDims;
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",&nDims);


	// Sanity check
	int totalNSubdomains = aiProd(nSubdomains,nDims);
	if(totalNSubdomains!=mpiSize)
		msg(ERROR|ONCE,"The product of grid:nSubdomains does not match the number of MPI processes");

	// Determine subdomain of this MPI node
	int *subdomain = malloc(nDims*sizeof(*subdomain));
	for(int d=0;d<nDims;d++){
		subdomain[d] = mpiRank % nSubdomains[d];
		mpiRank /= nSubdomains[d];
	}

	free(nSubdomains);
	return subdomain;

}

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

/******************************************************************************
 *	FINITE DIFFERENCE
 *****************************************************************************/

 void gFinDiff1st(const Grid *scalar, Grid *field){

 	//Performs first order centered finite difference on scalar and returns a field

 	int rank =scalar->rank;
 	long int *sizeProd = scalar->sizeProd;
 	long int *fieldSizeProd = field->sizeProd;

 	double *scalarVal = scalar->val;
 	double *fieldVal = field->val;

 	//Scalar indexes
 	long int sNext, sPrev;
 	long int f;
 	int fNext = fieldSizeProd[1];

 	//Centered Finite difference
 	for(int d = 1; d < rank; d++){
 		sNext = sizeProd[d];
 		sPrev = -sizeProd[d];
 		f = d-1;
 		for(int g = 0; g < sizeProd[rank]; g++){
 			fieldVal[f] = 0.5*(scalarVal[sNext] - scalarVal[sPrev]);
 			sNext++;
 			sPrev++;
 			f += fNext;
 		}
 	}

 	return;
 }

 void gFinDiff2nd2D(Grid *result, const Grid *object){

 	//Load
 	int rank = object->rank;
 	long int *sizeProd = object->sizeProd;

 	double *resultVal = result->val;
 	double *objectVal = object->val;

 	// Index of neighboring nodes
	long int g = sizeProd[1] + sizeProd[2];
 	long int gj = g + sizeProd[1];
 	long int gjj= g - sizeProd[1];
 	long int gk = g + sizeProd[2];
 	long int gkk= g - sizeProd[2];

	long int end = sizeProd[rank] - 2*g;

 	//Laplacian
 	for(int q = 0; q < end; q++){
 		resultVal[g] = -4.*objectVal[g];
 		resultVal[g] += objectVal[gj] + objectVal[gjj]
 						+objectVal[gk] + objectVal[gkk];

 		//Increment indexes
		g++;
 		gj++;
 		gjj++;
 		gk++;
 		gkk++;
 	}

 	return;
 }

void gFinDiff2nd3D(Grid *result, const  Grid *object){

 	//Load
 	int rank = object->rank;
 	long int *sizeProd = object->sizeProd;

 	double *resultVal = result->val;
 	double *objectVal = object->val;

 	// Index of neighboring nodes
 	long int g = sizeProd[1] + sizeProd[2] + sizeProd[3];
 	long int gj = g + sizeProd[1];
 	long int gjj= g - sizeProd[1];
 	long int gk = g + sizeProd[2];
 	long int gkk= g - sizeProd[2];
 	long int gl = g + sizeProd[3];
 	long int gll= g - sizeProd[3];

	long int end = sizeProd[rank] - 2*g;

 	//Laplacian
 	for(int q = 0; q < end; q++){
 		resultVal[g] = -6.*objectVal[g];
 		resultVal[g] += objectVal[gj] + objectVal[gjj]
 						+objectVal[gk] + objectVal[gkk]
 						+objectVal[gl] + objectVal[gll];

 		//Increment indexes
 		g++;
 		gj++;
 		gjj++;
 		gk++;
 		gkk++;
 		gl++;
 		gll++;
 	}

 	return;
 }

/******************************************************************************
 *	HALO FUNCTIONS
 *****************************************************************************/

void gHaloOp(SliceOpPointer sliceOp, Grid *grid, const MpiInfo *mpiInfo, int reverse){

	int rank = grid->rank;
	for(int d = 1; d < rank; d++){
		gHaloOpDim(sliceOp, grid, mpiInfo, d, reverse);
	}

}

void gHaloOpDim(SliceOpPointer sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, int inverse){

 	//Load MpiInfo
 	int mpiRank = mpiInfo->mpiRank;
 	int *subdomain = mpiInfo->subdomain;
 	int *nSubdomains = mpiInfo->nSubdomains;
 	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

	//Load
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	double *sendSlice = grid->sendSlice;
	double *recvSlice = grid->recvSlice;

	// Normal operation: take 2nd outermost layer and place it outermost
	// Inverse operation: take outermost layer and place it 2nd outermost
	int offsetUpperTake  = size[d]-2+inverse;
	int offsetUpperPlace = size[d]-1-inverse;
	int offsetLowerTake  =         1-inverse;
	int offsetLowerPlace =           inverse;

	//Dimension used for subdomains, 1 less entry than grid dimensions
	int dd = d - 1;
	int nSlicePoints = sizeProd[rank]/size[d];

	int firstElem = mpiRank - subdomain[dd]*nSubdomainsProd[dd];

 	// msg(STATUS,"lowerSubdomain: %i",lowerSubdomain);

	int upperSubdomain = firstElem
		+ ((subdomain[dd] + 1)%nSubdomains[dd])*nSubdomainsProd[dd];
	int lowerSubdomain = firstElem
		+ ((subdomain[dd] - 1 + nSubdomains[dd])%nSubdomains[dd])*nSubdomainsProd[dd];

	MPI_Status 	status;

	// TBD: Ommitting this seems to yield race condition between consecutive
	// calls to gHaloOpDim(). I'm not quite sure why so this should be
	// investigated further.
	MPI_Barrier(MPI_COMM_WORLD);

	// Send and recieve upper (tag 1)
	getSlice(sendSlice, grid, d, offsetUpperTake);
	MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 1,
                 recvSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 1,
                 MPI_COMM_WORLD, &status);
	sliceOp(recvSlice, grid, d, offsetLowerPlace);

	// Send and recieve lower (tag 0)
	getSlice(sendSlice, grid, d, offsetLowerTake);
	MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 0,
                 recvSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 0,
                 MPI_COMM_WORLD, &status);
	sliceOp(recvSlice, grid, d, offsetUpperPlace);

}


/*****************************************************************************
 *		ALLOC/DESTRUCTORS
 ****************************************************************************/

Grid *gAlloc(const dictionary *ini, int nValues){

	//Sanity check
	iniAssertEqualNElements(ini, 2,"grid:trueSize", "grid:stepSize");

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims, nBnd, nBoundaries;
	int *trueSizeTemp = iniGetIntArr(ini, "grid:trueSize", &nDims);
	int *nGhostLayersTemp = iniGetIntArr(ini, "grid:nGhostLayers", &nBoundaries);
	double *stepSizeTemp = iniGetDoubleArr(ini, "grid:stepSize", &nDims);
	char **boundaries = iniGetStrArr(ini, "grid:boundaries" , &nBnd);

	//Sanity check
	if(!(nBnd==(nDims)*2) && !(nBnd == 1)){
		msg(ERROR, "%d boundary edges specified, need %d, or 1",nBnd, nDims*2 );
	}

	//More sanity check
	if(nBoundaries != 2*nDims){
		msg(ERROR|ONCE, "Need ghost cells depth for all the boundaries: 2*nDims");
	}

	// Calculate the number of grid points (True points + ghost points)
	int rank = nDims+1;
	int *size 			= malloc(rank*sizeof(*size));
	int *trueSize 		= malloc(rank*sizeof(*trueSize));
	double *stepSize	= malloc(rank*sizeof(*stepSize));
	int *nGhostLayers 	= malloc(2*rank*sizeof(*nGhostLayers));

	size[0] = nValues;
	trueSize[0] = nValues;
	stepSize[0] = 1;
	nGhostLayers[0] = 0;
	nGhostLayers[rank] = 0;

	for(int d = 1 ; d < rank; d++){
		trueSize[d] = trueSizeTemp[d-1];
		stepSize[d] = stepSizeTemp[d-1];
		nGhostLayers[d] = nGhostLayersTemp[d-1];
		nGhostLayers[d+rank] = nGhostLayersTemp[d+nDims-1];

		size[d] = trueSize[d] + nGhostLayers[d] + nGhostLayers[d+rank];
	}
	free(trueSizeTemp);
	free(stepSizeTemp);
	free(nGhostLayersTemp);

	//Cumulative products
	long int *sizeProd = malloc((rank+1)*sizeof(*sizeProd));
	ailCumProd(size,sizeProd,rank);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=0;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//Memory for values and a slice
	double *val = malloc(sizeProd[rank]*sizeof(*val));
	double *sendSlice = malloc(nSliceMax*sizeof(*sendSlice));
	double *recvSlice = malloc(nSliceMax*sizeof(*recvSlice));
	double *bndSlice = malloc(2*rank*nSliceMax*sizeof(*bndSlice)); //Maybe seek
	//a different solution where it is only stored where needed

	//Set boundary conditions, should be cleaned up
	int inc = 1;
	bndType *bnd = malloc(2*rank*sizeof(*bnd));

	if(nBnd==1){
		for(int b = 0; b<2*rank; b++){
			if(!strcmp(boundaries[0], "PERIODIC")){
				bnd[b] = PERIODIC;
			} else if(!strcmp(boundaries[0], "DIRICHLET")){
				bnd[b] = DIRICHLET;
			} else if(!strcmp(boundaries[0], "NEUMANN")){
				bnd[b] = NEUMANN;
			}
		}
	} else {
		for(int b = 0; b < nBnd; b++){
			if(!strcmp(boundaries[b], "PERIODIC")){
				bnd[b+inc] = PERIODIC;
			} else if(!strcmp(boundaries[b], "DIRICHLET")){
				bnd[b+inc] = DIRICHLET;
			} else if(!strcmp(boundaries[b], "NEUMANN")){
				bnd[b+inc] = NEUMANN;
			}
			if(b == rank-2) inc = 2;
		}
	}
	bnd[0] = NONE;
	bnd[rank] = NONE;

	/* Store in Grid */
	Grid *grid = malloc(sizeof(*grid));

	grid->rank = rank;
	grid->size = size;
	grid->trueSize = trueSize;
	grid->sizeProd = sizeProd;
	grid->nGhostLayers = nGhostLayers;
	grid->stepSize = stepSize;
	grid->val = val;
	grid->h5 = 0;	// Must be activated separately
	grid->sendSlice = sendSlice;
	grid->recvSlice = recvSlice;
	grid->bndSlice = bndSlice;
	grid->bnd = bnd;

	return grid;
}

MpiInfo *gAllocMpi(const dictionary *ini){
	//Sanity check
	iniAssertEqualNElements(ini, 2,"grid:nSubdomains","grid:trueSize");

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims, temp;
	int nSpecies = iniGetNElements(ini, "population:nParticles");
	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", &nDims);
	int *nGhostLayers = iniGetIntArr(ini, "grid:nGhostLayers", &temp);
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", &nDims);
	int *nSubdomainsProd = malloc((nDims+1)*sizeof(*nSubdomainsProd));
	aiCumProd(nSubdomains,nSubdomainsProd,nDims);

	//Position of the subdomain in the total domain
	int *subdomain = getSubdomain(ini);
	int *offset = malloc(nDims*sizeof(*offset));
	double *posToSubdomain = malloc(nDims*sizeof(*posToSubdomain));

	for(int d = 0; d < nDims; d++){
		// offset[d] = subdomain[d]*trueSize[d];
		offset[d] = subdomain[d]*trueSize[d]-nGhostLayers[d];
		posToSubdomain[d] = (double)1/trueSize[d];
	}

    MpiInfo *mpiInfo = malloc(sizeof(*mpiInfo));
	mpiInfo->subdomain = subdomain;
	mpiInfo->nSubdomains = nSubdomains;
	mpiInfo->nSubdomainsProd = nSubdomainsProd;
	mpiInfo->offset = offset;
	mpiInfo->nDims = nDims;
	mpiInfo->posToSubdomain = posToSubdomain;
	mpiInfo->mpiSize = mpiSize;
	mpiInfo->mpiRank = mpiRank;

	mpiInfo->nSpecies = nSpecies;
	mpiInfo->nNeighbors = 0;	// Neighbourhood not created

	free(trueSize);

    return mpiInfo;
}

void gFreeMpi(MpiInfo *mpiInfo){

	free(mpiInfo->subdomain);
	free(mpiInfo->nSubdomains);
	free(mpiInfo->nSubdomainsProd);
	free(mpiInfo->offset);
	free(mpiInfo->posToSubdomain);
	free(mpiInfo);

}

void gFree(Grid *grid){

	free(grid->size);
	free(grid->trueSize);
	free(grid->sizeProd);
	free(grid->nGhostLayers);
	free(grid->stepSize);
	free(grid->val);
	free(grid->sendSlice);
	free(grid->recvSlice);
	free(grid->bnd);
	free(grid);

}

int *gGetGlobalSize(const dictionary *ini){

	int nDims;
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",&nDims);
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",&nDims);
	char *bnd = iniGetStr(ini,"grid:boundaries");

	int *L = malloc(nDims*sizeof(*L));

	if(!strcmp(bnd,"PERIODIC")){
		for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d];
	} else {
		msg(ERROR,"Only PERIODIC grid:boundaries supported yet");
	}

	free(trueSize);
	free(nSubdomains);
	free(bnd);

	return L;
}

void gSetBndSlices(Grid *grid,MpiInfo *mpiInfo){

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=0;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	double constant1 = 1.;
	double constant2 = 2.;

	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET)
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant1;
				}
			if(bnd[d] == NEUMANN)
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET)
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant1;

				}
			if(bnd[d] == NEUMANN)
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
		}
	}

	// adPrint(bndSlice, nSliceMax*rank);

	return;
}

/****************************************************************************
 *	CONVENIENCE GRID OPERATIONS
 ***************************************************************************/

void gMul(Grid *grid, double num){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	for(long int p=0;p<nElements;p++) grid->val[p] *= num;
}

void gAdd(Grid *grid, double num){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	for(long int p=0;p<nElements;p++) grid->val[p] += num;
}

void gSub(Grid *grid, double num){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	for(long int p=0;p<nElements;p++) grid->val[p] -= num;
}


void gZero(Grid *grid){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	for(long int p=0;p<nElements;p++) grid->val[p] = 0;
}

void gSet(Grid *grid, const double *value){

	long int *sizeProd = grid->sizeProd;
	int *size = grid->size;
	int rank = grid->rank;

	for(long int p=0;p<sizeProd[rank];p+=size[0])
		for(int pp=0;pp<size[0];pp++)
			grid->val[p+pp] = value[pp];
}

void gNormalizeE(const dictionary *ini, Grid *E){

	int nSpecies, nDims;
	double *q = iniGetDoubleArr(ini,"population:charge",&nSpecies);
	double *m = iniGetDoubleArr(ini,"population:mass",&nSpecies);
	double timeStep = iniGetDouble(ini,"time:timeStep");
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",&nDims);
	gMul(E,pow(timeStep,2)*(q[0]/m[0]));
	for(int p=0;p<E->sizeProd[E->rank];p++){
		E->val[p] /= stepSize[p%E->size[0]];
	}

}

//Not that well tested
void gNeutralizeGrid(Grid *grid, const MpiInfo *mpiInfo){

	const double *val = grid->val;
	long int *sizeProd = grid->sizeProd;
	int *trueSize = grid->trueSize;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	int mpiSize = mpiInfo->mpiSize;


	double myCharge = gNeutralizeGridInner(&val,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);
	double totCharge = 0;

	// MPI_Barrier(MPI_COMM_WORLD);

	// msg(STATUS, "Hello before");
	MPI_Allreduce(&myCharge, &totCharge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// msg(STATUS, "HELLO AFTER");

	// MPI_Barrier(MPI_COMM_WORLD);


	double avgCharge = totCharge/((double)aiProd(&trueSize[1] , rank-1)*mpiSize);

	gSub(grid, avgCharge);

	// avgCharge = gNeutralizeGridInner(&val,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);
	return;
}

static double gNeutralizeGridInner(	const double **val, const int *nGhostLayersBefore, const int *nGhostLayersAfter,
									const int *trueSize, const long int *sizeProd){

	double charge = 0.;

	if(*sizeProd==1){

		*val += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			charge += (*(*val)++);

		*val += *sizeProd**nGhostLayersAfter;

	} else {

		*val += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			charge += gNeutralizeGridInner(val,nGhostLayersBefore-1,nGhostLayersAfter-1,trueSize-1,sizeProd-1);

		*val += *sizeProd**nGhostLayersAfter;
	}

	return charge;
}



void gAddTo(Grid *result, Grid *addition){

	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	double *resultVal = result->val;
	double *addVal = addition->val;

	for(long int g = 0; g < sizeProd[rank]; g++)	resultVal[g] += addVal[g];

}

void gSubFrom(Grid *result, Grid *subtraction){

	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	double *resultVal = result->val;
	double *subVal = subtraction->val;

	for(long int g = 0; g < sizeProd[rank]; g++)	resultVal[g] -= subVal[g];

}


/****************************************************************
 *		Boundary conditions
 ***************************************************************/

static void gPeriodic(Grid *phi, const  MpiInfo *mpiInfo){
	// msg(STATUS, "Hello");
	// gNeutralizeGrid(phi, mpiInfo);

	return;
}

void gDirichlet(Grid *grid, const int boundary,  const  MpiInfo *mpiInfo){

	// msg(STATUS, "Hello from Dirichlet");

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice;

	//Compute dimensions and size of slice
	int d = boundary%rank;
	int offset = 1 + (boundary>rank)*(size[d]-2);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=0;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset);

	// adPrint(&bndSlice[(boundary-1)*nSliceMax], 10);

	return;

}

void gNeumann(Grid *grid, const int boundary, const MpiInfo *mpiInfo){

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice;
	double *slice = grid->sendSlice;

	//Compute dimensions and slicesize
	int d = boundary%rank;
	int offset = (boundary>rank)*(size[d]-1);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=0;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}
	//Compute d/dx u(x) = u(x_2) - 2A
	// constant *=-2;
	getSlice(slice, grid, d, offset + 2 - 4*(boundary>rank));

	for(int s = 0; s < nSliceMax; s++) slice[s] -=2*bndSlice[s+boundary*nSliceMax];

	setSlice(slice, grid, d, offset);

	return;
}

void gBnd(Grid *grid, const MpiInfo *mpiInfo){

	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == PERIODIC)	gPeriodic(grid, mpiInfo);
			else if(bnd[d] == DIRICHLET) gDirichlet(grid, d, mpiInfo);
			else if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == PERIODIC)	gPeriodic(grid, mpiInfo);
			else if(bnd[d] == DIRICHLET) gDirichlet(grid, d, mpiInfo);
			else if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}

	return;
}

/*****************************************************************************
 *		NEIGHBORHOOD
 ****************************************************************************/

void gCreateNeighborhood(const dictionary *ini, MpiInfo *mpiInfo, Grid *grid){

	// RETRIEVE NECESSARY VARIABLES

	int nDims = mpiInfo->nDims;
	int nSpecies = mpiInfo->nSpecies;
	int *size = grid->size;

	// COMPUTE SIMPLE VARIABLES

	int nNeighbors = pow(3,nDims);
	int neighborhoodCenter = 0;
	for(int i=0;i<nDims;i++) neighborhoodCenter += pow(3,i);

	// ALLOCATE FOR MIGRANTS AND FIND NUMBER TO ALLOCATE FOR

	int nTest;
	long int *nEmigrantsAllocTemp = iniGetLongIntArr(ini,"grid:nEmigrantsAlloc",&nTest);
	if(nTest!=nNeighbors && nTest!=1 && nTest!=nDims){
		msg(ERROR|ONCE,"grid:nEmigrantsAlloc must consist of 1, nDims=%i or 3^nDims=%i elements",nDims,nNeighbors);
	}
	long int *nEmigrantsAlloc = malloc(nNeighbors*sizeof(*nEmigrantsAlloc));

	// Set all migrant-buffers to the same size
	if(nTest==1){
		alSetAll(nEmigrantsAlloc,nNeighbors,nEmigrantsAllocTemp[0]);
		nEmigrantsAlloc[neighborhoodCenter] = 0;
	}

	// User has manually specified each buffer in lexicographical order
	if(nTest==nNeighbors){
		memcpy(nEmigrantsAlloc,nEmigrantsAllocTemp,nNeighbors*sizeof(*nEmigrantsAlloc));
		nEmigrantsAlloc[neighborhoodCenter] = 0;
	}

	// User has specified the buffers according to how many dimensions the
	// interface to the neighbour is (e.g. 0 for corners, 1 for edges, 2 for
	// faces) in increasing order
	if(nTest==nDims){
		for(int neigh=0;neigh<nNeighbors;neigh++){
			if(neigh==neighborhoodCenter) nEmigrantsAlloc[neigh] = 0;
			else {
				int temp = neigh;
				int interfaceDims = nDims;
				for(int d=nDims-1;d>=0;d--){
					int power = pow(3,d);
					if(temp/power!=1) interfaceDims--;
					temp %= power;
				}
				nEmigrantsAlloc[neigh] = nEmigrantsAllocTemp[interfaceDims];
			}
		}
	}

	long int **migrants = malloc(nNeighbors*sizeof(**migrants));
	long int **migrantsDummy = malloc(nNeighbors*sizeof(**migrantsDummy));
	double **emigrants = malloc(nNeighbors*sizeof(**emigrants));
	double **emigrantsDummy = malloc(nNeighbors*sizeof(**emigrantsDummy));
	for(int i=0;i<nNeighbors;i++)
		if(i!=neighborhoodCenter){
			migrants[i] = malloc(nEmigrantsAlloc[i]*sizeof(*migrants));
			emigrants[i] = malloc(2*nDims*nEmigrantsAlloc[i]*sizeof(*emigrants));
		}

	double *thresholds = iniGetDoubleArr(ini,"grid:thresholds",&nTest);
	if(nTest!=2*nDims){
		msg(ERROR|ONCE,"grid:threshold must be 2*nDims=%i elements", 2*nDims);
	}
	for(int i=0;i<2*nDims;i++){
		// if(thresholds[i]<0) thresholds[i] = size[i%nDims+1] + thresholds[i];
		if(thresholds[i]<0) thresholds[i] += (size[i%nDims+1]-1);
	}

	// ALLOCATE SIMPLE ARRAYS AND STORE IN STRUCT

	//long int *nMigrants = malloc(nNeighbors*nSpecies*sizeof(*nMigrants));
	long int *nEmigrants = malloc(nNeighbors*nSpecies*sizeof(*nEmigrants));
	long int *nImmigrants = malloc(nNeighbors*nSpecies*sizeof(*nImmigrants));

	long int nImmigrantsAlloc = 2*nDims*alMax(nEmigrantsAlloc,nNeighbors);
	double *immigrants = malloc(nImmigrantsAlloc*sizeof(*immigrants));

	MPI_Request *send = malloc(nNeighbors*sizeof(*send));
	MPI_Request *recv = malloc(nNeighbors*sizeof(*recv));
	for(int ne=0;ne<nNeighbors;ne++){
		send[ne] = MPI_REQUEST_NULL;
		recv[ne] = MPI_REQUEST_NULL;
	}


	mpiInfo->send = send;
	mpiInfo->recv = recv;
	mpiInfo->nNeighbors = nNeighbors;
	mpiInfo->migrants = migrants;
	mpiInfo->migrantsDummy = migrantsDummy;
	mpiInfo->emigrants = emigrants;
	mpiInfo->emigrantsDummy = emigrantsDummy;
	mpiInfo->nEmigrants = nEmigrants;
	mpiInfo->nImmigrants = nImmigrants;
	mpiInfo->nEmigrantsAlloc = nEmigrantsAlloc;
	mpiInfo->nImmigrantsAlloc = nImmigrantsAlloc;
	mpiInfo->thresholds = thresholds;
	mpiInfo->immigrants = immigrants;
	mpiInfo->neighborhoodCenter = neighborhoodCenter;

}

void gDestroyNeighborhood(MpiInfo *mpiInfo){

	long int **migrants = mpiInfo->migrants;
	double **emigrants = mpiInfo->emigrants;
	for(int neigh=0;neigh<mpiInfo->nNeighbors;neigh++){
		if(neigh!=mpiInfo->neighborhoodCenter){
			free(migrants[neigh]);
			free(emigrants[neigh]);
		}
	}
	free(migrants);
	free(emigrants);
	free(mpiInfo->migrantsDummy);
	free(mpiInfo->emigrantsDummy);
	mpiInfo->nNeighbors = 0;
	free(mpiInfo->nEmigrantsAlloc);
	free(mpiInfo->thresholds);
	free(mpiInfo->immigrants);
	free(mpiInfo->nImmigrants);
	free(mpiInfo->send);
	free(mpiInfo->recv);
}

/******************************************************************************
 * H5 FUNCTIONS
 *****************************************************************************/

void gWriteH5(const Grid *grid, const MpiInfo *mpiInfo, double n){

	hid_t fileSpace = grid->h5FileSpace;
	hid_t memSpace = grid->h5MemSpace;
	hid_t file = grid->h5;
	double *val = grid->val;

	// Enable collective datawriting
	hid_t pList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);

	char name[64];
	sprintf(name,"/n=%.1f",n);

	hid_t dataset = H5Dcreate(file,name,H5T_IEEE_F64LE,fileSpace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);

	H5Dclose(dataset);
	H5Pclose(pList);
}

void gReadH5(Grid *grid, const MpiInfo *mpiInfo, double n){

	hid_t fileSpace = grid->h5FileSpace;
	hid_t memSpace = grid->h5MemSpace;
	hid_t file = grid->h5;
	double *val = grid->val;

	// Enable collective datawriting
	hid_t pList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);

	char name[64];
	sprintf(name,"/n=%.1f",n);

	hid_t dataset = H5Dopen(file,name,H5P_DEFAULT);
	H5Dread(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);

	H5Dclose(dataset);
	H5Pclose(pList);

}

void gCloseH5(Grid *grid){
	H5Sclose(grid->h5MemSpace);
	H5Sclose(grid->h5FileSpace);
	H5Fclose(grid->h5);
}

void gOpenH5(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo,
						  const double *denorm, const double *dimen, const char *fName){

	int rank = grid->rank;
	int nDims = rank-1;
	int *size = grid->size;
	int *trueSize = grid->trueSize;
	int	*nGhostLayers = grid->nGhostLayers;
	int *nSubdomains = mpiInfo->nSubdomains;
	int *subdomain = mpiInfo->subdomain;

	/*
	 * CREATE FILE
	 */

	hid_t file = openH5File(ini,fName,"grid");

	/*
	 * CREATE ATTRIBUTES
	 */

	double *debye = malloc(nDims*sizeof(*debye));
	debye[0] = iniGetDouble(ini,"grid:debye");
	for(int d=1;d<nDims;d++) debye[d]=debye[0];

	setH5Attr(file,"Axis denormalization factor",&grid->stepSize[1],nDims);
	setH5Attr(file,"Axis dimensionalizing factor",debye,nDims);
	setH5Attr(file,"Quantity denormalization factor",denorm,size[0]);
	setH5Attr(file,"Quantity dimensionalizing factor",denorm,size[0]);

	/*
	 * HDF5 HYPERSLAB DEFINITION
	 */

	hsize_t *fileDims 	= malloc(rank*sizeof(*fileDims));
	hsize_t *memDims 	= malloc(rank*sizeof(*memDims));
	hsize_t *memOffset 	= malloc(rank*sizeof(*memOffset));
	hsize_t *fileOffset = malloc(rank*sizeof(*fileOffset));

	for(int d=0;d<rank;d++){
		// HDF5 indices needs to be reversed compared to ours due to non-C ordering.
		memDims[d]		= (hsize_t)size[rank-d-1];
		memOffset[d]	= (hsize_t)nGhostLayers[rank-d-1];
		fileDims[d]		= (hsize_t)trueSize[rank-d-1]*nSubdomains[rank-d-2];
		fileOffset[d]	= (hsize_t)trueSize[rank-d-1]*subdomain[rank-d-2];
	}

	fileDims[rank-1] = (hsize_t)trueSize[0];
	fileOffset[rank-1] = (hsize_t)0.;

	hid_t memSpace = H5Screate_simple(rank,memDims,NULL);
	for(int d=0;d<nDims;d++) memDims[d] = trueSize[rank-d-1];
	H5Sselect_hyperslab(memSpace,H5S_SELECT_SET,memOffset,NULL,memDims,NULL);

	hid_t fileSpace = H5Screate_simple(rank,fileDims,NULL);
	H5Sselect_hyperslab(fileSpace,H5S_SELECT_SET,fileOffset,NULL,memDims,NULL);

	free(fileDims);
	free(memDims);
	free(memOffset);
	free(fileOffset);

	grid->h5 = file;
	grid->h5MemSpace = memSpace;
	grid->h5FileSpace = fileSpace;

}

/******************************************************************************
 * ENERGY FUNCTIONS
 *****************************************************************************/

void gPotEnergy(const Grid *rho, const Grid *phi, Population *pop){

	const double *rhoVal = rho->val;
	const double *phiVal = phi->val;
	long int *sizeProd = rho->sizeProd;
	int *trueSize = rho->trueSize;
	int *nGhostLayers = rho->nGhostLayers;
	int rank = rho->rank;

	double energy = gPotEnergyInner(&rhoVal,&phiVal,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);

	int nSpecies = pop->nSpecies;
	pop->potEnergy[nSpecies] = energy;

}

static double gPotEnergyInner(	const double **rhoVal, const double **phiVal, const int *nGhostLayersBefore,
								const int *nGhostLayersAfter, const int *trueSize, const long int *sizeProd){

	double energy = 0;

	if(*sizeProd==1){

		*rhoVal += *sizeProd**nGhostLayersBefore;
		*phiVal += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			energy += (*(*rhoVal)++)*(*(*phiVal)++);

		*rhoVal += *sizeProd**nGhostLayersAfter;
		*phiVal += *sizeProd**nGhostLayersAfter;

	} else {

		*rhoVal += *sizeProd**nGhostLayersBefore;
		*phiVal += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			energy += gPotEnergyInner(rhoVal,phiVal,nGhostLayersBefore-1,nGhostLayersAfter-1,trueSize-1,sizeProd-1);

		*rhoVal += *sizeProd**nGhostLayersAfter;
		*phiVal += *sizeProd**nGhostLayersAfter;
	}

	return energy;
}

/*****************************************************************************
 *		MISC
 ****************************************************************************/

void gValDebug(Grid *grid, const MpiInfo *mpiInfo){

	int mpiRank = mpiInfo->mpiRank;
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	double *v = grid->val;

	for(long int p=0;p<sizeProd[rank];p++){
		v[p] = mpiRank*1000;
		long int temp = p;

		for(int d=0;d<rank;d++){
			v[p] += (temp%size[d])*pow(10,d-1);
			temp/=size[d];
		}
	}
}

/**************************************************************
*			TEMP, to reading of h5 files are ready
*************************************************************/

void fillHeaviside(Grid *grid, const MpiInfo *mpiInfo){

   //Load
   int *size = grid->size;
   int *trueSize = grid->trueSize;
   long int *sizeProd = grid->sizeProd;
   int *subdomain = mpiInfo->subdomain;
   int *nSubdomains = mpiInfo->nSubdomains;

   double *val = grid->val;

   gZero(grid);

   //Hardcoding try
   long int ind = 0;
   if(nSubdomains[1]==1){
	   //One core
	   for(int j = 1; j < size[1]-1; j++){
		   for (int k = 1; k<size[2]-1; k++) {
			   for(int l = 1; l < size[3]-1; l++){
				   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				   if(k < (trueSize[2]+1)/2.) val[ind] = 1.;
				   // else if (k == trueSize[2]/2 || k == trueSize[2]) val[ind] = 0.;
				   else val[ind] = -1.;
			   }
		   }
	   }
   } else {
	   msg(STATUS|ONCE, "Hello");
	   //Multi core
	   for(int j = 1; j < size[1]-1; j++){
		   for (int k = 1; k<size[2]-1; k++) {
			   for(int l = 1; l < size[3]-1; l++){
				   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				//    val[ind] = 1.;
				//    msg(STATUS|ONCE, "%d", ind);
				   if(subdomain[1]<nSubdomains[1]/2) val[ind] = 1.0;
				   else val[ind] = -1.;
			   }
		   }
	   }
	//    Set in 0 at between domains (sendSlice is safe to use, since it is reset every time it is used)
	   double *slice = grid->sendSlice;
	   for(int j = 0; j < size[1]; j++)	slice[j] = 0;
	   setSlice(slice, grid, 1, size[1]-2);
   }



   // msg(STATUS, "%f", adSum(val, sizeProd[4]));

   return;
}

void fillHeaviSol(Grid *grid, const MpiInfo *mpiInfo){

   //Load
   int *size = grid->size;
   int *trueSize = grid->trueSize;
   long int *sizeProd = grid->sizeProd;
   double *val =grid->val;

   //Hardcoding try
   long int ind = 0;
   //One core in z dir
   for(int j = 1; j < size[1]-1; j++){
	   for (int k = 1; k<size[2]-1; k++) {
		   for(int l = 1; l < size[3]-1; l++){
			   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
			   if(l < trueSize[3]/2)  val[ind] = -128. + (double)(l - 16.)*(l - 16.)/2.;
			   else if (l == trueSize[3]/2. || l == trueSize[3]) val[ind] = 0.;
			   else val[ind] = -128.;
		   }
	   }
   }

   return;
}

void fillPolynomial(Grid *grid , const MpiInfo *mpiInfo){

   int *size = grid->size;
   long int *sizeProd = grid->sizeProd;

   //Load GridQuantity
   double *val = grid->val;

   long int ind = 0;
   for(int j = 1; j < size[1]-1; j++){
	   for (int k = 1; k<size[2]-1; k++) {
		   for(int l = 1; l < size[3]-1; l++){
		   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
		   val[ind] = (double) j*j;
		   // val[ind] = 1.;

		   }
		   // msg(STATUS, "j = %d, k =%d", j,k);
	   }
   }

   return;

}



void fillPointCharge(Grid *grid , const MpiInfo *mpiInfo){

   //Grid info
   int *size = grid->size;
   // int *nGhostLayers = grid->nGhostLayers;
   long int *sizeProd = grid->sizeProd;
   // int rank = grid->rank;
   double *val = grid->val;

   double value = (double ) -1e2;

   gZero(grid);

   //Mpi info
   int mpiRank = mpiInfo->mpiRank;

   if(mpiRank == 0){	//Always placed at center of node 0;
	   long int ind =(size[1]/2)*sizeProd[1]
				   + (size[2]/2)*sizeProd[2]
				   + (size[3]/2)*sizeProd[3];
	   val[ind] = value;
   }

   return;

}

void fillPointSol(Grid *grid, const MpiInfo *mpiInfo){

   //Grid info
   int *size = grid->size;
   // int *nGhostLayers = grid->nGhostLayers;
   long int *sizeProd = grid->sizeProd;
   // int rank = grid->rank;
   double *val = grid->val;

   //Point of charge
   double x = (double)(size[1]/2);
   double y = (double)(size[2]/2);
   double z = (double)(size[3]/2);

   // msg(STATUS, "x,y,z = [%f, %f, %f]", x,y,z);

   //Mpi info
   int mpiRank = mpiInfo->mpiRank;

   if(mpiRank == 0){	//Always placed at center of node 0;
	   //Find node at center
	   for(int j = 1; j < size[1]-1; j ++){
		   for(int k = 1; k < size[2]-1; k ++){
			   for(int l = 1; l < size[3]-1; l++){
				   double distance = (j-x)*(j-x) + (k-y)*(k-y) + (l-z)*(l-z);
				   distance = sqrt(distance);
				   long int ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];

				   if(distance > 0.00001)
					   val[ind] = 1./distance;
				   else val[ind] = 0.;
			   }
		   }
	   }
   }

   return;
}

void fillSin(Grid *grid, const MpiInfo *mpiInfo){

   //Load MPI stuff
   int *subdomain = mpiInfo->subdomain;
   int *nSubdomains = mpiInfo->nSubdomains;

   //Load grid info
   int *trueSize = grid->trueSize;
   int *size = grid->size;
   int *nGhostLayers = grid->nGhostLayers;
   long int *sizeProd = grid->sizeProd;
   double *val = grid->val;

   //Temp quick functions and constants
   double sin(double);

   long int ind;
   double coeffX = 2*PI/((trueSize[1])*nSubdomains[0]);
   double coeffY = 2*PI/(trueSize[2]*nSubdomains[1]);
   // double coeffZ = PI/(size[3]*nSubdomains[2]);


   for(int j = 1; j < size[1]-nGhostLayers[5]; j++){
	   for (int k = 0; k<size[2]-nGhostLayers[4]; k++) {
		   for(int l = 0; l < size[3]-nGhostLayers[3]; l++){
			   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3]; //Position in subgrid
			   // ind += (size[1]*subdomain[0]) + (size[2]*subdomain[1]) + (size[3]*subdomain[2]); //Adding position in whole grid
			   val[ind] = 	sin( ( (j - nGhostLayers[1]) +trueSize[1]*subdomain[0] )*coeffX)
						   *sin( ( (k - nGhostLayers[2]) +trueSize[2]*subdomain[1] )*coeffY);
						   // *sin((l+size[3]*subdomain[2])*coeffZ);

		   }
	   }
   }


   return;
}

void fillSinSol(Grid *grid, const MpiInfo *mpiInfo){

   //Load MPI stuff
   int *subdomain = mpiInfo->subdomain;
   int *nSubdomains = mpiInfo->nSubdomains;

   //Load grid info
   int *trueSize = grid->trueSize;
   int *size = grid->size;
   int *nGhostLayers = grid->nGhostLayers;
   long int *sizeProd = grid->sizeProd;
   double *val = grid->val;

   //Temp quick functions and constants
   double sin(double);

   long int ind;
   double coeffX = 2*PI/(trueSize[1]*nSubdomains[0]);
   double coeffY = 2*PI/(trueSize[2]*nSubdomains[1]);
   // double coeffZ = PI/(size[3]*nSubdomains[2]);


   for(int j = 1; j < size[1]-nGhostLayers[5]; j++){
	   for (int k = 0; k<size[2]-nGhostLayers[4]; k++) {
		   for(int l = 0; l < size[3]-nGhostLayers[3]; l++){
			   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3]; //Position in subgrid
			   // ind += (size[1]*subdomain[0]) + (size[2]*subdomain[1]) + (size[3]*subdomain[2]); //Adding position in whole grid
			   val[ind] = 	-1/(coeffX*coeffX)*sin( ( (j-nGhostLayers[1]) +trueSize[1]*subdomain[0] )*coeffX)
						   -1/(coeffY*coeffY)*sin( ( (k-nGhostLayers[2]) +trueSize[2]*subdomain[1] )*coeffX);
						   // *sin((l+size[3]*subdomain[2])*coeffZ);

		   }
	   }
   }


   return;
}

void fillExp(Grid *grid, const MpiInfo *mpiInfo){

   //Load Mpi
   // int *subdomain = mpiInfo->subdomain;
   // int *nSubdomains = mpiInfo->nSubdomains;

   //Load
   int *trueSize = grid->size;
   long int *sizeProd = grid->sizeProd;

   //Load GridQuantity
   double *val = grid->val;

   //Function
   double exp(double x);
   double halfPoint = trueSize[1]/2;
   double normalize = 1./(trueSize[1]*trueSize[1]);

   long int ind = 0;
   for(int j = 0; j < trueSize[1]; j++){
	   for (int k = 0; k<trueSize[2]; k++) {
		   for(int l = 0; l < trueSize[3]; l++){
			   ind = (j+1)*sizeProd[1] + (k+1)*sizeProd[2] + (l+1)*sizeProd[3];
			   val[ind] = exp(-10*((j-halfPoint)*(j-halfPoint)*normalize +
								   (k-halfPoint)*(k-halfPoint)*normalize +
									   (l-halfPoint)*(l-halfPoint)*normalize));
		   }
	   }
   }

   return;
}

void fillRng(Grid *grid, const MpiInfo *mpiInfo, const gsl_rng *rng){

   //Load
   double *val = grid->val;
   long int *sizeProd = grid->sizeProd;
   int rank = grid->rank;
   for(int g = 0; g < sizeProd[rank]; g++) val[g] = gsl_ran_gaussian_ziggurat (rng,1.);

   return;
}

void fillCst(Grid *grid, const MpiInfo *mpiInfo){

   int rank = grid->rank;
   long int *sizeProd = grid->sizeProd;
   double *val = grid->val;

   for(int g = 0; g < sizeProd[rank]; g++)  val[g] = 1.;


   return;
}
