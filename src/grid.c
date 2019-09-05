/**
 * @file		grid.c
 * @brief		Grid-struct handling.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
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
 * LOCAL FUNCTION DECLARATIONS
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

static double gPotEnergyInner(	const double **rhoVal, const double **phiVal,
								const int *nGhostLayersBefore, const int *nGhostLayersAfter,
								const int *trueSize, const long int *sizeProd);

static double gNeutralizeGridInner(	const double **val,
									const int *nGhostLayersBefore, const int *nGhostLayersAfter,
									const int *trueSize, const long int *sizeProd);

static void gContractInner(	const double **in, double **out,
							const int *layersBefore, const int *layersAfter,
	 						const int *trueSize, const long int *sizeProd);

static void gExpandInner(	const double **in, double **out,
							const int *layersBefore, const int *layersAfter,
	 						const int *trueSize, const long int *sizeProd);

/******************************************************************************
 * LOCAL FUNCTION DEFINITIONS
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
	int nDims = iniGetInt(ini,"grid:nDims");
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);


	// Sanity check
	int totalNSubdomains = aiProd(nSubdomains,nDims);
	if(totalNSubdomains!=mpiSize)
		msg(ERROR,"The product of grid:nSubdomains does not match the number of MPI processes");

	// Determine subdomain of this MPI node
	int *subdomain = malloc(nDims*sizeof(*subdomain));
	for(int d=0;d<nDims;d++){
		subdomain[d] = mpiRank % nSubdomains[d];
		mpiRank /= nSubdomains[d];
	}

	free(nSubdomains);
	return subdomain;

}

static void gContractInner(	const double **in, double **out,
							const int *layersBefore, const int *layersAfter,
	 						const int *trueSize, const long int *sizeProd){

	*in += *sizeProd**layersBefore;
	if(*sizeProd==1){
		for(int j=0;j<*trueSize;j++){
			*(*out)++ = *(*in)++;
		}
	} else {
		for(int j=0;j<*trueSize;j++){
			gContractInner(	in, out,
							layersBefore-1, layersAfter-1,
							trueSize-1, sizeProd-1);
		}
	}
	*in += *sizeProd**layersAfter;

}

static void gExpandInner(	const double **in, double **out,
							const int *layersBefore, const int *layersAfter,
	 						const int *trueSize, const long int *sizeProd){

	*out -= *sizeProd**layersAfter;
	if(*sizeProd==1){
		for(int j=0;j<*trueSize;j++){
			*(*out)-- = *(*in)--;
		}
	} else {
		for(int j=0;j<*trueSize;j++){
			gExpandInner(	in, out,
							layersBefore-1, layersAfter-1,
							trueSize-1, sizeProd-1);
		}
	}
	*out -= *sizeProd**layersBefore;

}

/******************************************************************************
 * GLOCAL FUNCTION DEFINITIONS
 *****************************************************************************/

/******************************************************************************
 *	FINITE DIFFERENCE
 *****************************************************************************/

 void gFinDiff1st(const Grid *scalar, Grid *field){

	// Performs first order centered finite difference on scalar and returns a field

	int rank = scalar->rank;
	// int *size = scalar->size;
	long int *sizeProd = scalar->sizeProd;
	long int *fieldSizeProd = field->sizeProd;

	double *scalarVal = scalar->val;
	double *fieldVal = field->val;

 	// Scalar indices
	long int sNext, sPrev;
	long int f;
	int fNext = fieldSizeProd[1];

	long int start = alSum(&sizeProd[1], rank-1 );
	long int end = sizeProd[rank]-start;


	// Centered Finite difference
	for(int d = 1; d < rank; d++){
		sNext = start + sizeProd[d];
		sPrev = start - sizeProd[d];
		f = start*fieldSizeProd[1] + (d-1);


		for(int g = start; g < end; g++){
			fieldVal[f] = 0.5*(scalarVal[sNext] - scalarVal[sPrev]);
			sNext++;
			sPrev++;
			f += fNext;
		}
	}
}


void gFinDiff2ndND(Grid *result, const Grid *object){

	// Load
	int rank = object->rank;
	long int *sizeProd = object->sizeProd;

	double *resultVal = result->val;
	double *objectVal = object->val;

	// Index of neighboring nodes
	long int g = alSum(&sizeProd[1], rank-1 );
	long int end = sizeProd[rank] - 2*g;
	int gStep;

	double coeff = 2.*(rank-1);

	// Laplacian
	for(int q = 0; q < end +1; q++){
		resultVal[g] = -coeff*objectVal[g];
		for(int r = 1; r < rank; r++){
			gStep = sizeProd[r];
			resultVal[g] += objectVal[g + gStep] + objectVal[g - gStep];
		}


		// Increment indices
		g++;
	}

	return;
}

void gFinDiff2nd3D(Grid *result, const  Grid *object){

 	// Load
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

 	// Laplacian
 	for(int q = 0; q < end; q++){
 		resultVal[g] = -6.*objectVal[g];
 		resultVal[g] += objectVal[gj] + objectVal[gjj]
 						+objectVal[gk] + objectVal[gkk]
 						+objectVal[gl] + objectVal[gll];

 		// Increment indices
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

void gHaloOp(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, opDirection dir){

	int rank = grid->rank;
	for(int d = 1; d < rank; d++){
		gHaloOpDim(sliceOp, grid, mpiInfo, d, dir);
	}

}

void gHaloOpDim(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, opDirection dir){

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

	// dir=TOHALO=0: take 2nd outermost layer and place it outermost
	// dir=FROMHALO=1: take outermost layer and place it 2nd outermost
	int offsetUpperTake  = size[d]-2+dir;
	int offsetUpperPlace = size[d]-1-dir;
	int offsetLowerTake  =         1-dir;
	int offsetLowerPlace =           dir;

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

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims = iniGetInt(ini, "grid:nDims");
	int *trueSizeTemp = iniGetIntArr(ini, "grid:trueSize", nDims);
	int *nGhostLayersTemp = iniGetIntArr(ini, "grid:nGhostLayers", 2*nDims);
	char **boundaries = iniGetStrArr(ini, "grid:boundaries" , 2*nDims);

	// Calculate the number of grid points (True points + ghost points)
	int rank = nDims+1;
	int *size 			= malloc(rank*sizeof(*size));
	int *trueSize 		= malloc(rank*sizeof(*trueSize));
	int *nGhostLayers 	= malloc(2*rank*sizeof(*nGhostLayers));

	if(nValues==VECTOR) nValues = nDims; // VECTOR equals -1

	size[0] = nValues;
	trueSize[0] = nValues;
	nGhostLayers[0] = 0;
	nGhostLayers[rank] = 0;

	for(int d = 1 ; d < rank; d++){
		trueSize[d] = trueSizeTemp[d-1];
		nGhostLayers[d] = nGhostLayersTemp[d-1];
		nGhostLayers[d+rank] = nGhostLayersTemp[d+nDims-1];

		size[d] = trueSize[d] + nGhostLayers[d] + nGhostLayers[d+rank];
	}
	free(trueSizeTemp);
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

	// Memory for values and a slice
	double *val = malloc(sizeProd[rank]*sizeof(*val));
	double *sendSlice = malloc(nSliceMax*sizeof(*sendSlice));
	double *recvSlice = malloc(nSliceMax*sizeof(*recvSlice));
	double *bndSlice = malloc(2*rank*nSliceMax*sizeof(*bndSlice));
	// Maybe seek a different solution where it is only stored where needed

	bndType *bnd = malloc(2*rank*sizeof(*bnd));
	int b = 0;
	for(int r=0; r<2*rank; r++){
		if(r%rank==0){
			bnd[r] = NONE;
		} else {
			if(		!strcmp(boundaries[b], "PERIODIC"))		bnd[r] = PERIODIC;
			else if(!strcmp(boundaries[b], "DIRICHLET"))	bnd[r] = DIRICHLET;
			else if(!strcmp(boundaries[b], "NEUMANN"))		bnd[r] = NEUMANN;
			else msg(ERROR,"%s invalid value for grid:boundaries",boundaries[b]);
			b++;
		}
	}

	/* Store in Grid */
	Grid *grid = malloc(sizeof(*grid));

	grid->rank = rank;
	grid->size = size;
	grid->trueSize = trueSize;
	grid->sizeProd = sizeProd;
	grid->nGhostLayers = nGhostLayers;
	grid->val = val;
	grid->h5 = 0;	// Must be activated separately
	grid->sendSlice = sendSlice;
	grid->recvSlice = recvSlice;
	grid->bndSlice = bndSlice;
	grid->bnd = bnd;

	return grid;
}

MpiInfo *gAllocMpi(const dictionary *ini){

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims = iniGetInt(ini, "grid:nDims");
	int nSpecies = iniGetInt(ini, "population:nSpecies");
	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", nDims);
	int *nGhostLayers = iniGetIntArr(ini, "grid:nGhostLayers", 2*nDims);
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims);
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
	free(grid->val);
	free(grid->sendSlice);
	free(grid->recvSlice);
	free(grid->bnd);
	free(grid);

}

int *gGetGlobalSize(const dictionary *ini){

	int nDims = iniGetInt(ini,"grid:nDims");
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
	char *bnd = iniGetStr(ini,"grid:boundaries");

	int *L = malloc(nDims*sizeof(*L));

	if(!strcmp(bnd,"PERIODIC")){
		for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d];
	} else {
		msg(ERROR,"Only all PERIODIC grid:boundaries supported by gGetGlobalSize() yet");
	}

	free(trueSize);
	free(nSubdomains);
	free(bnd);

	return L;
}

long int gGetGlobalVolume(const dictionary *ini){

	int nDims = iniGetInt(ini,"grid:nDims");
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	free(trueSize);

	int *L = gGetGlobalSize(ini);
	long int V = aiProd(L,nDims);
	free(L);

	return V;

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

void gSquare(Grid *grid){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	double *val = grid->val;
	for(long int g=0;g<nElements;g++) val[g] = val[g]*val[g];

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

void gCopy(const Grid *original, Grid *copy){
	//Load
	int rank = original->rank;
	long int *sizeProd = original->sizeProd;

	double *origVal =	original->val;
	double *copyVal=	copy->val;

	for(int g = 0; g < sizeProd[rank]; g++) copyVal[g] = origVal[g];

}

//Not that well tested
void gNeutralizeGrid(Grid *grid, const MpiInfo *mpiInfo){
	const double *val = grid->val;
	long int *sizeProd = grid->sizeProd;
	int *trueSize = grid->trueSize;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	int mpiSize = mpiInfo->mpiSize;



	double myCharge = gNeutralizeGridInner(&val,&nGhostLayers[rank-1],
					&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);
	double totCharge = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&myCharge, &totCharge, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double avgCharge = totCharge/((double)aiProd(&trueSize[1] , rank-1)*mpiSize);

	gSub(grid, avgCharge);

	return;
}

static inline double gNeutralizeGridInner(	const double **val, const int *nGhostLayersBefore, const int *nGhostLayersAfter,
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

void gSubFrom(Grid *result, const Grid *subtraction){

	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	double *resultVal = result->val;
	double *subVal = subtraction->val;

	for(long int g = 0; g < sizeProd[rank]; g++)	resultVal[g] -= subVal[g];

}

static double gSumTruegridInner(const double **val, const int *nGhostLayersBefore,
								const int *nGhostLayersAfter, const int *trueSize,
								const long int *sizeProd){

	double sum = 0;

	if(*sizeProd==1){

		*val += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			sum += (*(*val)++);

		*val += *sizeProd**nGhostLayersAfter;

	} else {

		*val += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			sum += gSumTruegridInner(val,nGhostLayersBefore-1,
									nGhostLayersAfter-1,trueSize-1,sizeProd-1);

		*val += *sizeProd**nGhostLayersAfter;
	}

	return sum;
}

double gSumTruegrid(const Grid *grid){

	//Load
	const double *val = grid->val;
	long int *sizeProd = grid->sizeProd;
	int *trueSize = grid->trueSize;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;

	double sum = gSumTruegridInner(&val,&nGhostLayers[rank-1],
					&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);

	return sum;

}

long int gTotTruesize(const Grid *grid, const MpiInfo *mpiInfo){

	int *nSubdomains = mpiInfo->nSubdomains;
	int *trueSize = grid->trueSize;
	int rank = grid->rank;

	long int totTruesize = 1;

	for(int r = 1; r < rank; r++) 	totTruesize *= nSubdomains[r-1]*trueSize[r];

	return totTruesize;
}

void gAssertNeutralGrid(const Grid *rho, const MpiInfo *mpiInfo){

	double sum = gSumTruegrid(rho);
	double totSum = 1.;
	MPI_Allreduce(&sum, &totSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	if( totSum < -0.001 || totSum > 0.001) msg(ERROR, "Total charge is %f", totSum);
}


void gRemoveHalo(Grid *grid){

	double *oldVal = grid->val;
	double *newVal = grid->val;
	int rank = grid->rank;
	int *nGhostLayers = grid->nGhostLayers;
	int *trueSize = grid->trueSize;

	long int *sizeProd = grid->sizeProd;

	gContractInner(	(const double **)&oldVal, &newVal,
					&nGhostLayers[rank-1], &nGhostLayers[2*rank-1],
					&trueSize[rank-1], &sizeProd[rank-1]);

	memcpy(grid->size,trueSize,rank*sizeof(*trueSize));
	ailCumProd(trueSize,sizeProd,rank);
	aiSetAll(nGhostLayers,2*rank,0);

}

void gInsertHalo(Grid *grid, const int *nGhostLayers){

	double *val = grid->val;
	int rank = grid->rank;
	int *size = grid->size;
	int *trueSize = grid->trueSize;

	memcpy(grid->nGhostLayers, nGhostLayers, 2*rank*sizeof(*nGhostLayers));

	for(int r=0;r<rank;r++){
		size[r] += nGhostLayers[r] + nGhostLayers[r+rank];
	}

	double *lastOldVal = &val[grid->sizeProd[rank]-1];
	ailCumProd(size,grid->sizeProd,rank);
	double *lastNewVal = &val[grid->sizeProd[rank]-1];

	long int *sizeProd = grid->sizeProd;

	gExpandInner(	(const double **)&lastOldVal, &lastNewVal,
					&nGhostLayers[rank-1], &nGhostLayers[2*rank-1],
					&trueSize[rank-1], &sizeProd[rank-1]);

}


/****************************************************************
 *		Boundary conditions
 ***************************************************************/

static void gPeriodic(Grid *phi, const  MpiInfo *mpiInfo){
	// msg(STATUS, "Hello");
	gNeutralizeGrid(phi, mpiInfo);

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

	//If periodic neutralize phi
	int periodic = 0;
	for(int d = 1; d < rank; d++){
		if(bnd[d] == PERIODIC)	periodic = 1;
	}
	if(periodic)	gPeriodic(grid, mpiInfo);

	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET) gDirichlet(grid, d, mpiInfo);
			else if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET) gDirichlet(grid, d, mpiInfo);
			if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
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

	int nTest = iniGetNElements(ini,"grid:nEmigrantsAlloc");
	if(nTest!=nNeighbors && nTest!=1 && nTest!=nDims){
		msg(ERROR,"grid:nEmigrantsAlloc must consist of 1, nDims=%i or 3^nDims=%i elements",nDims,nNeighbors);
	}

	long int *nEmigrantsAllocTemp = iniGetLongIntArr(ini,"grid:nEmigrantsAlloc",nTest);
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

	double *thresholds = iniGetDoubleArr(ini,"grid:thresholds",2*nDims);

	// upper thresholds should be counted from upper edge
	for(int i=nDims;i<2*nDims;i++){
		thresholds[i] = (size[i%nDims+1]-1) - thresholds[i];
	}

	// ALLOCATE SIMPLE ARRAYS AND STORE IN STRUCT

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
			 const Units *units, double denorm, const char *fName){

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

	setH5Attr(file,"Axis denormalization factor",&units->length,1);
	setH5Attr(file,"Quantity denormalization factor",&denorm,1);

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

	double energy = gPotEnergyInner(&rhoVal, &phiVal,
									&nGhostLayers[rank-1],
									&nGhostLayers[2*rank-1],
									&trueSize[rank-1], &sizeProd[rank-1]);
	energy *= 0.5;

	int nSpecies = pop->nSpecies;
	pop->potEnergy[nSpecies] = energy;

}

static double gPotEnergyInner(	const double **rhoVal, const double **phiVal,
								const int *nGhostLayersBefore, const int *nGhostLayersAfter,
								const int *trueSize, const long int *sizeProd){

	double energy = 0;

	*rhoVal += *sizeProd**nGhostLayersBefore;
	*phiVal += *sizeProd**nGhostLayersBefore;

	if(*sizeProd==1){
		for(int j=0;j<*trueSize;j++){
			energy += (*(*rhoVal)++)*(*(*phiVal)++);
		}
	} else {
		for(int j=0;j<*trueSize;j++){
			energy += gPotEnergyInner(	rhoVal, phiVal,
										nGhostLayersBefore-1, nGhostLayersAfter-1,
										trueSize-1, sizeProd-1);
		}
	}

	*rhoVal += *sizeProd**nGhostLayersAfter;
	*phiVal += *sizeProd**nGhostLayersAfter;

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

void gFillHeavi(Grid *grid, int d, const MpiInfo *mpiInfo){

   //Load
   int *size = grid->size;
   int rank = grid->rank;
   long int *sizeProd = grid->sizeProd;
   int *subdomain = mpiInfo->subdomain;
   int *nSubdomains = mpiInfo->nSubdomains;

   double *val = grid->val;

   gZero(grid);

   //Hardcoding try
   long int ind = 0;
   if(nSubdomains[d-1]==1){
	   //SetSlices
	   double *slice = grid->sendSlice;

	   int nSlicePoints = sizeProd[rank]/size[d];

	   for(int j = 0; j < nSlicePoints; j++)
	   		slice[j] = 1.;

	   for(int j = 2; j < size[d]/2; j++)
	   		setSlice(slice, grid, d, j);

	   for(int j = 0; j < nSlicePoints; j++)
	   		slice[j] = -1.;

	   for(int j = size[d]/2 + 1; j < size[d]; j++)
	   		setSlice(slice, grid, d, j);

   } else {
	   //Multi core
	   for(int j = 1; j < size[1]-1; j++){
		   for (int k = 1; k<size[2]-1; k++) {
			   for(int l = 1; l < size[3]-1; l++){
				   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				   if(subdomain[d-1]<nSubdomains[d-1]/2) val[ind] = 1.0;
				   else val[ind] = -1.;
			   }
		   }
	   }
	// Set in 0 at between domains (sendSlice is safe to use, since it is reset every time it is used)
	   double *slice = grid->sendSlice;
	   for(int j = 0; j < sizeProd[4]; j++)	slice[j] = 0.;
	   if(subdomain[d-1] == 0)	setSlice(slice, grid, d, 1);

	   for(int j = 0; j < sizeProd[4]; j++)	slice[j] = -0.;
	   if(subdomain[d-1] == nSubdomains[d-1]/2)	setSlice(slice, grid, d, 1);
   }

	gHaloOp(setSlice, grid, mpiInfo, 0);

   return;
}

void gFillHeaviSol(Grid *grid, int d ,const MpiInfo *mpiInfo){

	//Load
    int *size = grid->size;
    int *trueSize = grid->trueSize;
	int rank = grid->rank;
    long int *sizeProd = grid->sizeProd;
    int *subdomain = mpiInfo->subdomain;
    int *nSubdomains = mpiInfo->nSubdomains;

	double *val = grid->val;


   //Hardcoding try
   long int ind = 0;
   if(nSubdomains[d-1]==1){
	   //Smart setSlice-use
	   double half = 0.5*trueSize[d];
	   double *sol = malloc(trueSize[d]*sizeof(*sol));
	   double *slice = grid->sendSlice;
	   long int nSlicePoints = sizeProd[rank]/size[d];

	   //First half  f = -(a-x)*x/2
	   for(int j = 0; j < trueSize[d]/2; j++)
	   		sol[j] = 0.5*(half - j)*j;

		//Second half	f = (a - x)*x/2
 	   for(int j = trueSize[d]/2; j < trueSize[d]; j++)
			sol[j] = -0.5*(half - (j-half))*(j-half);

	   for(int j = 1; j < trueSize[d] + 1; j++){
		   for(int k = 0; k < nSlicePoints; k ++){
			   slice[k] = sol[j-1];
		   }
		   setSlice(slice, grid, d, j);
	   }

   } else {
	   //Multi core
	   long int advance;	//Total domain position
	   int half = nSubdomains[d-1]/2 * trueSize[d];
	   for(int j = 1; j < size[1]-1; j++){
		   for (int k = 1; k<size[2]-1; k++) {
			   for(int l = 1; l < size[3]-1; l++){
				   ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				   if(subdomain[d-1]<nSubdomains[d-1]/2){
					   //Rewrite to use d
					   advance = (j*(d==1) + k*(d==2) + l*(d==3))-1
					   			+ subdomain[d-1]*trueSize[d];
					   val[ind] = -0.5*(half - advance)*advance;
				   }
				   else{
					   advance = (j*(d==1) + k*(d==2) + l*(d==3))-1
					    		+ trueSize[1]* (subdomain[0]%(nSubdomains[0]/2));
					   val[ind] = -0.5*(advance - half)*advance;
				   }
			   }
		   }
	   }
   }

   gHaloOp(setSlice, grid, mpiInfo, 0);
   return;
}

void gFillPolynomial(Grid *grid , const MpiInfo *mpiInfo){

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
	   }
   }

   return;

}

void gFillPoint(Grid *grid , const MpiInfo *mpiInfo){

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

void gFillPointSol(Grid *grid, const MpiInfo *mpiInfo){

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

void gFillSin(Grid *grid,int d, const MpiInfo *mpiInfo, int norm){

	//MPI
	int *nSubdomains 	= mpiInfo->nSubdomains;
	int *subdomain 		= mpiInfo->subdomain;

	//Load
	int *size = grid->size;
	int *trueSize = grid->trueSize;
	int rank = grid->rank;
	long int *sizeProd = grid->sizeProd;

	double *sol = malloc(trueSize[d]*sizeof(*sol));
	double *slice = grid->sendSlice;
	long int nSlicePoints = sizeProd[rank]/size[d];

	//f = sin(x*2pi/L)
	double coeff = 2*PI/((trueSize[d])*nSubdomains[d-1]);
	double coeffNorm;

	//If it should be normalized to give E or phi correct
	if(norm){
		coeffNorm = coeff;
	} else {
		coeffNorm = coeff*coeff;
	}


	int J = subdomain[d-1]*trueSize[d];

	for(int j = 0; j < trueSize[d]; j++){
		sol[j] = coeffNorm*sin(J*coeff);
		J++;
	}

	for(int j = 1; j < trueSize[d]+1; j++){
		for(int k = 0; k < nSlicePoints; k ++){
			slice[k] = sol[j-1];
		}
		setSlice(slice, grid, d, j);
	}

   gHaloOp(setSlice, grid,mpiInfo, TOHALO);

   return;
}

void gFillSinSol(Grid *grid, int d ,const MpiInfo *mpiInfo){

	//MPI
	int *nSubdomains 	= mpiInfo->nSubdomains;
	int *subdomain 		= mpiInfo->subdomain;

	//Load
	int *size = grid->size;
	int *trueSize = grid->trueSize;
	int rank = grid->rank;
	long int *sizeProd = grid->sizeProd;

	double *sol = malloc(trueSize[d]*sizeof(*sol));
	double *slice = grid->sendSlice;
	long int nSlicePoints = sizeProd[rank]/size[d];

	//f = sin(x/2piL)
	double coeff = 2*PI/((trueSize[d])*nSubdomains[d-1]);

	int J = subdomain[d-1]*trueSize[d];

	for(int j = 0; j < trueSize[d]; j++){
		sol[j] = sin(J*coeff);
		J++;
	}

	for(int j = 1; j < trueSize[d] + 1; j++){
		for(int k = 0; k < nSlicePoints; k ++){
			slice[k] = sol[j-1];
		}
		setSlice(slice, grid, d, j);
	}

	gHaloOp(setSlice, grid,mpiInfo, TOHALO);

   return;
}

void gFillSinESol(Grid *grid, int d ,const MpiInfo *mpiInfo){

	//MPI
	int *nSubdomains 	= mpiInfo->nSubdomains;
	int *subdomain 		= mpiInfo->subdomain;

	//Load
	int *size = grid->size;
	int *trueSize = grid->trueSize;
	int rank = grid->rank;
	long int *sizeProd = grid->sizeProd;

	double *sol = calloc(sizeProd[1]*trueSize[d],sizeof(*sol));
	double *slice = grid->sendSlice;
	long int nSlicePoints = sizeProd[rank]/size[d];

	//f = cos(x/2piL)
	double coeff = 2*PI/((trueSize[d])*nSubdomains[d-1]);

	int J = subdomain[d-1]*trueSize[d];

	msg(STATUS, "sizeProd[1] = %d", sizeProd[1]);

	for(int j = 0; j < trueSize[d]; j++){
		sol[j] = cos(J*coeff);
		J++;
	}

	for(int j = 1; j < trueSize[d] + 1; j++){
		// int j = 1;
		for(int k = 0; k <  nSlicePoints; k ++){
			if((k+d-1)%3==0)	slice[k] = sol[j-1];
			else slice[k] = 0;
		}
		setSlice(slice, grid, d, j);
	}

	gHaloOp(setSlice, grid,mpiInfo, TOHALO);

	free(sol);
}

void gFillExp(Grid *grid, const MpiInfo *mpiInfo){

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

void gFillRng(Grid *grid, const MpiInfo *mpiInfo, const gsl_rng *rng){

   //Load
   double *val = grid->val;
   long int *sizeProd = grid->sizeProd;
   int rank = grid->rank;
   for(int g = 0; g < sizeProd[rank]; g++) val[g] = gsl_ran_gaussian_ziggurat (rng,1.);

   return;
}

void gFillCst(Grid *grid, const MpiInfo *mpiInfo){

   int rank = grid->rank;
   long int *sizeProd = grid->sizeProd;
   double *val = grid->val;

   for(int g = 0; g < sizeProd[rank]; g++)  val[g] = 1.;


   return;
}
