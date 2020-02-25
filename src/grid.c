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

static double *getSliceInner(double *nextGhost, double **valp, const long int *mul,
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

	double *val = grid->val;

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



// test


static const double *addSliceInnerAvg(const double *nextGhost, double **valp, const long int *mul,
	const int *points, const long int finalMul){

		if(*mul==finalMul){
			for(int j=0;j<*mul;j++) *((*valp)) = (*((*valp)++) + *(nextGhost++))/2.;
			*valp += (*mul)*(*points-1);
		} else {
			for(int j=0; j<*points;j++)
				nextGhost = addSliceInner(nextGhost, valp, mul-1,points-1,finalMul);
		}
		return nextGhost;

}

void addSliceAvg(const double *slice, Grid *grid, int d, int offset){

	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;

	double *val = grid->val;

	val += offset*sizeProd[d];
	addSliceInner(slice, &val, &sizeProd[rank-1], &size[rank-1], sizeProd[d]);
}

//test end


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
			//printf("node: %li, using %li, and %li \n",f,sNext,sPrev);
		}
	}

	return;
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
 		gHaloOpDim(sliceOp, grid, mpiInfo, d, dir); //lower
 		//printf("EXCHANGING \n");
 	}

 }

 void gHaloOpDim(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, opDirection dir){

  	//Load MpiInfo
  	int mpiRank = mpiInfo->mpiRank;
  	int *subdomain = mpiInfo->subdomain;
  	int *nSubdomains = mpiInfo->nSubdomains;
  	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;
 	bndType *bnd = grid->bnd;

 	//printf("in gHaloOpDimOpen rank = %i \n",mpiRank);

 	//MPI_Barrier(MPI_COMM_WORLD);
 	//Load
 	int rank = grid->rank;
 	int *size = grid->size;
 	long int *sizeProd = grid->sizeProd;
 	double *sendSlice = grid->sendSlice;
 	double *recvSlice = grid->recvSlice;

 	// dir=TOHALO=0: take 2nd outermost layer and place it outermost
 	// dir=FROMHALO=1: take outermost layer and place it 2nd outermost
 	//if (boundary == 1){ //upper
 		int offsetUpperTake  = size[d]-2+dir;
 		int offsetUpperPlace = size[d]-1-dir;
 	//}
 	//if (boundary == 1){ //lower
 		int offsetLowerTake  =         1-dir;
 		int offsetLowerPlace =           dir;
 	//}
 	//Dimension used for subdomains, 1 less entry than grid dimensions
 	int dd = d - 1;
 	int nSlicePoints = sizeProd[rank]/size[d];

 	int firstElem = mpiRank - subdomain[dd]*nSubdomainsProd[dd];

  	// msg(STATUS,"lowerSubdomain: %i",lowerSubdomain);

 	int upperSubdomain = firstElem
 		+ ((subdomain[dd] + 1)%nSubdomains[dd])*nSubdomainsProd[dd];
 	int lowerSubdomain = firstElem
 		+ ((subdomain[dd] - 1 + nSubdomains[dd])%nSubdomains[dd])*nSubdomainsProd[dd];

   //msg(STATUS,"lowerSubdomain: %i, upperSubdomain: %i",lowerSubdomain,upperSubdomain);
 	MPI_Status 	status;

 	// TBD: Ommitting this seems to yield race condition between consecutive
 	// calls to gHaloOpDim(). I'm not quite sure why so this should be
 	// investigated further.


 	MPI_Barrier(MPI_COMM_WORLD);

 	// Send and recieve upper (tag 1)
 	//printf("exchanging slices");
 	//printf("d = %i, bnd[rank+d] = %i, bnd[d] = %u \n",d,bnd[rank+d],bnd[d]);
 	//printf("d = %i, bnd[d] = %u \n",d,bnd[d]);

 	if(bnd[rank+d] == PERIODIC){
 		// Send and recieve lower (tag 0)

 			//printf("rank = %i, sending on edge %i \n",mpiInfo->mpiRank,d);

 		getSlice(sendSlice, grid, d, offsetLowerTake);
 		MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 0,
 	                 recvSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 0,
 	                 MPI_COMM_WORLD, &status);
 		sliceOp(recvSlice, grid, d, offsetUpperPlace);
 	}else{
 		getSlice(sendSlice, grid, d, offsetLowerTake);
 		MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 0,
 	                 recvSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 0,
 	                 MPI_COMM_WORLD, &status);//handshake
 	}
 	if(bnd[d] == PERIODIC){

 			//printf("rank = %i, sending on edge %i \n",mpiInfo->mpiRank,rank+d);

 		getSlice(sendSlice, grid, d, offsetUpperTake);
 		MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 1,
                  recvSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 1,
                  MPI_COMM_WORLD, &status);
 	  sliceOp(recvSlice, grid, d, offsetLowerPlace);
 	}else{
 		getSlice(sendSlice, grid, d, offsetUpperTake);
 		MPI_Sendrecv(sendSlice, nSlicePoints, MPI_DOUBLE, upperSubdomain, 1,
                  recvSlice, nSlicePoints, MPI_DOUBLE, lowerSubdomain, 1,
                  MPI_COMM_WORLD, &status);//handshake
 	}
 	//printf("in gHaloOpDimOpen rank = %i \n",mpiRank);
 }




/*****************************************************************************
 *		ALLOC/DESTRUCTORS
 ****************************************************************************/

 Grid *gAlloc(const dictionary *ini, int nValues, const MpiInfo *mpiInfo){

 	// Get MPI info
 	//int mpiRank = mpiInfo->mpiRank;
 	//int mpiSize = mpiInfo->mpiSize;
 	int *subdomain = mpiInfo->subdomain;
 	int *nSubdomains = mpiInfo->nSubdomains;
 	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

 	int mpiSize, mpiRank;
 	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
 	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

 	// Load data from ini
 	int nDims = iniGetInt(ini, "grid:nDims");
 	int *trueSizeTemp = iniGetIntArr(ini, "grid:trueSize", nDims);
 	int *nGhostLayersTemp = iniGetIntArr(ini, "grid:nGhostLayers", 2*nDims);
 	char **boundaries = iniGetStrArr(ini, "grid:boundaries" , 2*nDims);

 	//printf("boundaries =%s,%s,%s,%s,%s,%s, \n",boundaries[0],boundaries[1],boundaries[2],boundaries[3],boundaries[4],boundaries[5]);
	//printf("truesize =%i,%i,%i, \n",trueSizeTemp[0],trueSizeTemp[1],trueSizeTemp[2]);
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
 		//printf("size[d] =%i trueSize[d] =%i nGhostLayers[d] =%i nGhostLayers[d+rank] = %i \n",size[d],trueSize[d],nGhostLayers[d],nGhostLayers[d+rank]);
 	}
 	free(trueSizeTemp);
 	free(nGhostLayersTemp);

 	//Cumulative products
 	long int *sizeProd = malloc((rank+1)*sizeof(*sizeProd));
 	ailCumProd(size,sizeProd,rank);



 	//Number of elements in slice
 	long int nSliceMax = 0;
 	for(int d=1;d<rank;d++){
 		long int nSlice = 1;
 		for(int dd=0;dd<rank;dd++){
 			if(dd!=d) nSlice *= size[dd]*sizeProd[1];
 		}
 		if(nSlice>nSliceMax) nSliceMax = nSlice;
 	}
 	//msg(STATUS,"nSliceMax = %li",nSliceMax);
	//printf("sizeProd[1] = %li\n",sizeProd[1] );

 	// Memory for values and a slice
 	double *val = malloc(sizeProd[rank]*sizeof(*val));
 	double *sendSlice = malloc(nSliceMax*sizeof(*sendSlice));
 	double *recvSlice = malloc(nSliceMax*sizeof(*recvSlice));
 	double *bndSlice = malloc(2*rank*nSliceMax*sizeProd[1]*sizeof(*bndSlice));
	adSetAll(bndSlice,2*rank*nSliceMax*sizeProd[1],0);
	//double *bndSolution = malloc(2*rank*nSliceMax*sizeof(*bndSolution));
 	//printf("alloc sizeProd[rank] = %li\n",sizeProd[rank]);
 	// Maybe seek a different solution where it is only stored where needed


 	//Load
 	//int rank = grid->rank;
 	//long int *sizeProd = grid->sizeProd;

 	//Dimension used for subdomains, 1 less entry than grid dimensions
 	//int dd = d - 1;

 	bndType *bnd = malloc(2*rank*sizeof(*bnd));

 	bnd[0] = NONE; //initialize
 	bnd[rank] = NONE; //initialize


 	int b = 0;
 	for (int d = 1; d<rank;d++){
 		//msg(STATUS,"b=%i, d = %i, rank = %i",b,d,rank);
 		int dd = d - 1;
 		//msg(STATUS,"mpiRank = %i, subdomain[dd] = %i nSubdomainsProd[dd] = %i",mpiRank, subdomain[dd], nSubdomainsProd[dd]);
 		int firstElem = mpiRank - subdomain[dd]*nSubdomainsProd[dd];

 		// msg(STATUS,"lowerSubdomain: %i",lowerSubdomain);

 		int upperSubdomain = firstElem
 			+ ((subdomain[dd] + 1)%nSubdomains[dd])*nSubdomainsProd[dd];
 		int lowerSubdomain = firstElem
 			+ ((subdomain[dd] - 1 + nSubdomains[dd])%nSubdomains[dd])*nSubdomainsProd[dd];

 		//printf("rank: %i, lowerSubdomain: %i, upperSubdomain = %i \n",mpiRank,lowerSubdomain,upperSubdomain);

 		//lower
 		//for(int r=dd; r<dd+1; r++){
 		//printf("b = %i \n",b);
 		int r=dd+1;
 		if(lowerSubdomain>=mpiRank){
 				if(		!strcmp(boundaries[b], "PERIODIC"))		bnd[r] = PERIODIC;
 				else if(!strcmp(boundaries[b], "DIRICHLET"))	bnd[r] = DIRICHLET;
 				else if(!strcmp(boundaries[b], "NEUMANN"))		bnd[r] = NEUMANN;
 				else msg(ERROR,"%s invalid value for grid:boundaries",boundaries[b]);

 			}else if (lowerSubdomain<mpiRank){
 				//printf("YEPS!");
 				bnd[r] = PERIODIC;
 			} else {
 				bnd[r] = NONE; //initialize
 				//printf("bnd[%i] = NONE",r);

 		}
 		//upper
 		//for(int r=rank+dd; r<rank+dd+1; r++){
 		r = rank+dd+1;
 		//printf("r = %i \n",r);
 		if(upperSubdomain<=mpiRank){
 				if(		!strcmp(boundaries[b+rank-1], "PERIODIC"))		bnd[r] = PERIODIC;
 				else if(!strcmp(boundaries[b+rank-1], "DIRICHLET"))	bnd[r] = DIRICHLET;
 				else if(!strcmp(boundaries[b+rank-1], "NEUMANN"))		bnd[r] = NEUMANN;
 				else msg(ERROR,"%s invalid value for grid:boundaries",boundaries[b]);

 			}else if(upperSubdomain>mpiRank){
 				bnd[r] = PERIODIC;
 			}else {
 				bnd[r] = NONE; //initialize
 				//printf("bnd[%i] = NONE",r);
 			}
 			b++;
 		//}
 	}
 	//msg(ERROR,"%d, %d, %d, %d, %d, %d,%d, %d, ",bnd[0], bnd[1],bnd[2],bnd[3],bnd[4],bnd[5],bnd[6],bnd[7]);
 	//printf("in Grid rank = %i, %d, %d, %d, %d, %d, %d,%d, %d \n",mpiRank,bnd[0], bnd[1],bnd[2],bnd[3],bnd[4],bnd[5],bnd[6],bnd[7]);
 	//printf("rank = %i,EXITING \n",mpiRank);
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
 	//grid->bndSolution = bndSolution;
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
	//int *trueSizeAlloc = malloc(nDims*sizeof(*trueSizeAlloc));
	bool periodic = false;

	for(int d = 0; d < nDims; d++){
		// offset[d] = subdomain[d]*trueSize[d];
		offset[d] = subdomain[d]*trueSize[d]-nGhostLayers[d];
		posToSubdomain[d] = (double)1/trueSize[d];
		//trueSizeAlloc[d] = trueSize[d];
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
	mpiInfo->trueSize = trueSize;

	mpiInfo->nSpecies = nSpecies;
	mpiInfo->nNeighbors = 0;	// Neighbourhood not created
	mpiInfo->periodic = periodic;

	//free(trueSize);

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
	//free(grid->bndSolution);
	free(grid);

}

int *gGetGlobalSize(const dictionary *ini){

	int nDims = iniGetInt(ini,"grid:nDims");
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
	char **bnd = iniGetStrArr(ini, "grid:boundaries" , 2*nDims);

	int *L = malloc(nDims*sizeof(*L));

	//TODO: Handling of mixed boundaries
	if(!strcmp(bnd[0],"PERIODIC")){
		for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d];

	}else if(!strcmp(bnd[0],"NEUMANN")){
		for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d];
	}else if(!strcmp(bnd[0],"DIRICHLET")){
		for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d];
	} else {
		msg(ERROR,"Only all PERIODIC, DIRICHLET or NEUMANN grid:boundaries supported by gGetGlobalSize() yet");
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

void gSetBndSlices(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo){

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *nGhostLayers = grid->nGhostLayers;
	//double *bndSolution = grid->bndSolution;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.
	int nSpecies = iniGetInt(ini,"population:nSpecies");
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
	double *B = iniGetDoubleArr(ini,"fields:BExt",nDims);

	double veld[nDims];

	for (int s=0;s<nSpecies;s++){
		for (int d = 0;d<nDims;d++){
			//msg(STATUS,"d = %i, d+s*nDims = %i",d,d+s*nDims);
			veld[d] += (1./nSpecies)*velDrift[d+s*nDims];
			}
	}

	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veld[0],veld[1],veld[2]);

    //double B[3] = {1., 0., 0.};
	//double veld[3] = {0., 1., 1.};
	double veldCrossB[3] = {0., 0., 0.};
	adCrossProd(veld, B, veldCrossB);

	//printf("B[0] = %f, B[1] = %f, B[2] = %f \n",B[0],B[1],B[2]);
	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veldCrossB[0],veldCrossB[1],veldCrossB[2]);
	//printf("veldCrossB = %f,%f,%f",veldCrossB[0],veldCrossB[1],veldCrossB[2]);


	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//double constant1 = 0.; //Dirichlet
	double constant2 = 0.; // solution to bnd cond = 0.

	// The dirichlet part in this function ended up realy ugly. I am sorry for
	// that, but it works.

	// indices are indices of slice perp to d
	long int indices[rank];
	long int edge[rank]; //tells us what dimensions are perp (slice)

	//Lower edge
	for(int d = 1; d < rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-1));
			//printf("q = %li, indices[q] = %li \n",q,indices[q]);
		}
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = nGhostLayers[q+1];
					} else{
						// start indices at minimum local indice
						indices[q] = 0;//
					}
				}
				for(int s = 0; s < nSliceMax; s++){

					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = nGhostLayers[q+1];
						}
					}

					bndSlice[s + (nSliceMax * d)] = 0;
					for(int dd = 0;dd<rank-1;dd++){ //dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + (nSliceMax * d)] += veldCrossB[dd]*(indices[dd]+(subdomain[dd]*size[dd+1])-(2*subdomain[dd]-(nSubdomains[dd]-1))) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
						//																veldCrossB[dd]*(indices[dd]+(subdomain[dd]*size[dd+1])+(!edge[dd])*(subdomain[dd]+(nSubdomains[dd]-1))+0.5-3*edge[dd]*subdomain[dd]) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
						//printf("subdomain[%i] = %i, nSubdomains[%i] = %i \n",dd,subdomain[dd],dd,nSubdomains[dd]);
						if(veldCrossB[dd]*indices[dd]*edge[dd]!=0){
						}
					}
					// counter to increment only in the slice dims, and not the
					// dim perp to slice
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 domain
						if(indices[dd]<(size[dd+1])-1 && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//nGhostLayers[dd+1];
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}
			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					//Solution to equation. constant for now

				}
			}
		}
	}

	//Upper edge
	for(int d = rank+1; d < 2*rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-rank-1));
			//printf("d = %i, edge[q] = %li \n",d,edge[q]);
		}
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
						//printf("nSubdomains = %i\n",nSubdomains[q]);
					} else{
						// start indices at minimum
						indices[q] = 0;
					}
				}
				for(int s = 0; s < nSliceMax; s++){
					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
							//printf("nSubdomains = %i\n",nSubdomains[q]);
						}
					}
					bndSlice[s + (nSliceMax * (d))] = 0;
					for(int dd = 0;dd<rank-1;dd++){
						bndSlice[s + (nSliceMax * (d))] +=  veldCrossB[dd]*((indices[dd])+(subdomain[dd]*size[dd+1])-(subdomain[dd])-(subdomain[dd]-(nSubdomains[dd]-1))) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
																								//veldCrossB[dd]*(indices[dd]+(subdomain[dd]*size[dd+1])+0.5-(subdomain[dd])) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
					}
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 local domain
						if(indices[dd]<(size[dd+1]-1) && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//(size[dd+1]-nGhostLayers[dd]);
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}

			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	//msg(STATUS,"nSliceMax = %li",nSliceMax);
	//adPrint(&bndSlice[nSliceMax], nSliceMax*(rank));
	free(velDrift);
	//free(B);
	return;
}











void gSetBndSlicesE(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo){

	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *nGhostLayers = grid->nGhostLayers;
	//double *bndSolution = grid->bndSolution;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.
	int nSpecies = iniGetInt(ini,"population:nSpecies");
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
	double *B = iniGetDoubleArr(ini,"fields:BExt",nDims);

	double veld[nDims];

	for (int s=0;s<nSpecies;s++){
		for (int d = 0;d<nDims;d++){
			//msg(STATUS,"d = %i, d+s*nDims = %i",d,d+s*nDims);
			veld[d] += (1./nSpecies)*velDrift[d+s*nDims];
			}
	}

	double veldCrossB[3] = {0., 0., 0.};
	adCrossProd(veld, B, veldCrossB);

	//exit(0);
	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd]*sizeProd[1];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//double constant1 = 0.; //Dirichlet
	double constant2 = 0.; // solution to bnd cond = 0.

	//Lower edge
	for(int d = 1; d < rank; d++){

		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){

				for(int s = 0; s < nSliceMax; s+=sizeProd[1]){

					for(int dd = 0;dd<rank-1;dd++){ //dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + (nSliceMax * d)+dd] = -veldCrossB[dd];
						//																veldCrossB[dd]*(indices[dd]+(subdomain[dd]*size[dd+1])+(!edge[dd])*(subdomain[dd]+(nSubdomains[dd]-1))+0.5-3*edge[dd]*subdomain[dd]) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
						//printf("subdomain[%i] = %i, nSubdomains[%i] = %i \n",dd,subdomain[dd],dd,nSubdomains[dd]);
					}
				}
			}
			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					//Solution to equation. constant for now

				}
			}
		}
	}

	//Upper edge
	for(int d = rank+1; d < 2*rank; d++){

		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET){

				for(int s = 0; s < nSliceMax; s+=sizeProd[1]){

					for(int dd = 0;dd<rank-1;dd++){ //dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + (nSliceMax * d)+dd] = -veldCrossB[dd];
						//																veldCrossB[dd]*(indices[dd]+(subdomain[dd]*size[dd+1])+(!edge[dd])*(subdomain[dd]+(nSubdomains[dd]-1))+0.5-3*edge[dd]*subdomain[dd]) -0.5*veldCrossB[dd]*size[dd+1]*nSubdomains[dd];
						//printf("subdomain[%i] = %i, nSubdomains[%i] = %i \n",dd,subdomain[dd],dd,nSubdomains[dd]);
					}
				}
			}

			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	//msg(STATUS,"nSliceMax = %li",nSliceMax);
	//adPrint(&bndSlice[nSliceMax], nSliceMax*(rank));
	free(velDrift);
	//free(B);
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
	//msg(STATUS,"neutralize grid");
	const double *val = grid->val;
	long int *sizeProd = grid->sizeProd;
	int *trueSize = grid->trueSize;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	int mpiSize = mpiInfo->mpiSize;



	double myCharge = gNeutralizeGridInner(&val,&nGhostLayers[rank-1],
					&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1]);
	double totCharge = 0;

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

void gMulTo(Grid *result, Grid *addition){

	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	double *resultVal = result->val;
	double *addVal = addition->val;
	for(long int g = 0; g < sizeProd[rank]; g++)	resultVal[g] *= addVal[g];

}

void gDivTo(Grid *result, Grid *addition){

	int rank = result->rank;
	long int *sizeProd = result->sizeProd;
	double *resultVal = result->val;
	double *addVal = addition->val;
	for(long int g = 0; g < sizeProd[rank]; g++)	resultVal[g] /= addVal[g];

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

		for(int j=0;j<*trueSize;j++){
			//printf("*val = %f \n",**val);
			sum += (*(*val)++);
		}


		*val += *sizeProd**nGhostLayersAfter;

	} else {

		*val += *sizeProd**nGhostLayersBefore;

		for(int j=0;j<*trueSize;j++)
			sum += gSumTruegridInner(val,nGhostLayersBefore-1,
									nGhostLayersAfter-1,trueSize-1,sizeProd-1);

		*val += *sizeProd**nGhostLayersAfter;
	}
	//exit(0);

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

	//msg(STATUS, "Hello from Dirichlet");

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	double *bndSlice = grid->bndSlice;

	//Compute dimensions and size of slice
	int d = boundary%rank;
	int offset = 1 + (boundary>rank)*(size[d]-3);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd]*sizeProd[1];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}
	//msg(STATUS,"offset. eg. index to set slice in perp. direction %i",offset);
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset); //edge before halo
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset - 1 + (boundary>rank)*2); //halo
	//setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset + 1 - (boundary>rank)*2); //edge before edge

	//adPrint(&bndSlice[(boundary)*nSliceMax], nSliceMax);

	return;

}


void gDirichletE(Grid *grid, const int boundary,  const  MpiInfo *mpiInfo){

	//msg(STATUS, "Hello from Dirichlet");

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	double *bndSlice = grid->bndSlice;

	//Compute dimensions and size of slice
	int d = boundary%rank;
	int offset = 1 + (boundary>rank)*(size[d]-3);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd]*sizeProd[1];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}
	//msg(STATUS,"offset. eg. index to set slice in perp. direction %i",offset);
	//setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset); //edge before halo
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset - 1 + (boundary>rank)*2); //halo
	//setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset + 1 - (boundary>rank)*2); //edge before edge

	//adPrint(&bndSlice[(boundary)*nSliceMax], nSliceMax);

	return;

}



void gNeumann(Grid *grid, const int boundary, const MpiInfo *mpiInfo){

	//msg(STATUS, "Hello from NEUMANN");

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice; // two slices in each dim
	double *slice = grid->sendSlice;
	//double *bndSolution = grid->bndSolution; // two slices in each dim


	//Compute dimensions and slicesize
	int d = boundary%rank;
	int offset = 1 + (boundary>rank)*(size[d]-3);

	//int offset = 0;
	//if (boundary>rank) offset = 1 + (boundary>rank)*(size[d]-1);
	//if (boundary<rank) offset = (boundary>rank)*(size[d]-1);

	//Number of elements in slice
	long int nSliceMax = 0; //TODO: this can/should be stored.
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}
	//  /// OLD comment: Compute d/dx u(x) = u(x_2) - 2A
	// constant *=-2;

	// if (boundary>rank){ //Upper
	// 	//getSlice(slice, grid, d, offset - 1);
	// 	//setSlice(slice, grid, d, offset + 1);
	// 	getSlice(&bndSlice[boundary*nSliceMax], grid, d, offset);
	// }
	// if (boundary<rank) { //Lower
	// 	//getSlice(slice, grid, d, offset + 1);
	// 	//setSlice(slice, grid, d, offset - 1);
	// 	getSlice(&bndSlice[boundary*nSliceMax], grid, d, offset);
	// }
	// //adPrint(bndSolution,2*4*nSliceMax);
	getSlice(&bndSlice[boundary*nSliceMax], grid, d, offset); //edge before halo
	//getSlice(&bndSlice[boundary*nSliceMax], grid, d, offset + 1 - (boundary>rank)*2); //edge before edge
	// for(long int k = boundary*nSliceMax;k<nSliceMax;k++){
	// 	bndSlice[k] *= 2*bndSlice[k];
	// }

	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset - 1 + (boundary>rank)*2); //halo

	// if (boundary>rank){
	// 	for(int s = 0; s < nSliceMax; s++){
	// 		slice[s] -= 2.*0;//bndSolution[s+2*rank*nSliceMax-nSliceMax]; //bndSlice[s];
	// 		//msg(STATUS,"2*rank*nSliceMax = %li",2*rank*nSliceMax);
	// 		bndSlice[s + d*nSliceMax] = slice[s];
	// 		//msg(STATUS,"slice[s] = %f",slice[s]);
	// 	}
	// 	setSlice(slice, grid, d, offset + 1); //halo
	// 	//for(int s = 0; s < nSliceMax; s++){
	//  	//bndSlice[s + 2*d*nSliceMax] = 0;//slice[s];
	//  	//}
	// 	//adPrint(&bndSlice[2*d*nSliceMax],nSliceMax);
	// }
	// if (boundary<rank){
	// 	for(int s = 0; s < nSliceMax; s++){
	// 		slice[s] = -2.*0;//bndSolution[s+d*nSliceMax-nSliceMax]; //bndSlice[s];
	// 		bndSlice[s + d*nSliceMax] = slice[s];
	// 		//msg(STATUS,"slice[s] = %f",slice[s]);
	// 	}
	//	setSlice(slice, grid, d, offset - 1); //halo
		//for(int s = 0; s < nSliceMax; s++){
		//bndSlice[s + d*nSliceMax] = 0;//slice[s];

	// }
	 //adPrint(&bndSlice[d*nSliceMax],nSliceMax);
	//}
	//gNeutralizeGrid(grid, mpiInfo);
	//adPrint(bndSlice,2*4*nSliceMax);
	//adPrint(grid->val,grid->sizeProd[4]);
	//msg(STATUS,"done");
	//exit(0);
	return;
}




void gBnd(Grid *grid, const MpiInfo *mpiInfo){


	//printf("before boundary cond, rank %i\n",mpiInfo->mpiRank);
	//msg(STATUS,"phi size = %i",grid->sizeProd[4]);
	//adPrint(grid->val,grid->sizeProd[4]);
	//adPrint(grid->val,grid->sizeProd[4] );
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	bool periodic = mpiInfo->periodic;

	//If periodic neutralize phi

	// bool periodic =  true;
	// for(int d = 1; d < rank; d++){
	// 	//msg(STATUS,"d = %i",d);
	// 	if(bnd[d] != PERIODIC){
	// 		//msg(STATUS,"bnd[d] != PERIODIC, d = %i",d);
	// 		periodic = false;
	// 		}
	// }
	// for(int d = rank+1; d < 2*rank; d++){
	// 	//msg(STATUS,"d = %i",d);
	// 	if(bnd[d] != PERIODIC){
	// 		//msg(STATUS,"bnd[d] != PERIODIC, d = %i",d);
	// 		periodic = false;
	// 		}
	// }
	//printf("periodic = %d \n",mpiInfo->periodic);
	if(periodic == true){
		//printf("PERIODIC cond, rank %i\n",mpiInfo->mpiRank);
		gPeriodic(grid, mpiInfo);

	}
	//gHaloOp(setSlice, grid,mpiInfo,TOHALO);
	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				//msg(STATUS,"bnd[d] = DIRICHLET, giving d = %i, rank = %i",d,rank);
				gDirichlet(grid, d, mpiInfo);
				//msg(STATUS,"bnd[d] = DIRICHLET, d = %i",d);
			}
			else if(bnd[d] == NEUMANN){
				//msg(STATUS,"bnd[d] = NEUMANN, d = %i",d);
				gNeumann(grid, d, mpiInfo);
			}
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET) gDirichlet(grid, d, mpiInfo);
			if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}
	//printf("after boundary cond, rank %i\n",mpiInfo->mpiRank);
	//msg(STATUS,"phi size = %i",grid->sizeProd[4]);
	//if (mpiInfo->mpiRank == 7){
		//adPrint(grid->val,grid->sizeProd[4] );
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//exit(1);
	return;
}


void gBndE(Grid *grid, const MpiInfo *mpiInfo){


	//printf("before boundary cond, rank %i\n",mpiInfo->mpiRank);
	//msg(STATUS,"phi size = %i",grid->sizeProd[4]);
	//adPrint(grid->val,grid->sizeProd[4]);
	//adPrint(grid->val,grid->sizeProd[4] );
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	bool periodic = mpiInfo->periodic;


	//gHaloOp(setSlice, grid,mpiInfo,TOHALO);
	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				//msg(STATUS,"bnd[d] = DIRICHLET, giving d = %i, rank = %i",d,rank);
				//gDirichletE(grid, d, mpiInfo);
				gNeumann(grid, d, mpiInfo);
				//msg(STATUS,"bnd[d] = DIRICHLET, d = %i",d);
			}
			else if(bnd[d] == NEUMANN){
				//msg(STATUS,"bnd[d] = NEUMANN, d = %i",d);
				gNeumann(grid, d, mpiInfo);
			}
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET) gNeumann(grid, d, mpiInfo);//gDirichletE(grid, d, mpiInfo);
			if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}
	//printf("after boundary cond, rank %i\n",mpiInfo->mpiRank);
	//msg(STATUS,"phi size = %i",grid->sizeProd[4]);
	//if (mpiInfo->mpiRank == 7){
		//adPrint(grid->val,grid->sizeProd[4] );
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//exit(1);
	return;
}


/*****************************************************************************
 *		NEIGHBORHOOD
 ****************************************************************************/

void gCreateNeighborhood(const dictionary *ini, MpiInfo *mpiInfo, Grid *grid){

	// RETRIEVE NECESSARY VARIABLES

	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	int nDims = mpiInfo->nDims;
	int nSpecies = mpiInfo->nSpecies;
	int *size = grid->size;

	int *periodicAll = malloc(mpiSize*sizeof(*periodicAll));
	bndType *bnd = grid->bnd;

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
			emigrants[i] = malloc(2*nSpecies*nDims*nEmigrantsAlloc[i]*sizeof(*emigrants));
		}

	double *thresholds = iniGetDoubleArr(ini,"grid:thresholds",2*nDims);

	// upper thresholds should be counted from upper edge
	for(int i=nDims;i<2*nDims;i++){
		thresholds[i] = (size[i%nDims+1]-1) - thresholds[i];
	}

	// ALLOCATE SIMPLE ARRAYS AND STORE IN STRUCT

	long int *nEmigrants = malloc(nNeighbors*(nSpecies+1)*sizeof(*nEmigrants));
	long int *nImmigrants = malloc(nNeighbors*(nSpecies+1)*sizeof(*nImmigrants));

	long int nImmigrantsAlloc = 2*alMax(nEmigrantsAlloc,nNeighbors);
	// TODO: Investigate: It seems to be more stable with the 2*, not shure why
	double *immigrants = malloc(2*nSpecies*nDims*nImmigrantsAlloc*sizeof(*immigrants));

	MPI_Request *send = malloc(nNeighbors*sizeof(*send));
	MPI_Request *recv = malloc(nNeighbors*sizeof(*recv));
	for(int ne=0;ne<nNeighbors;ne++){
		send[ne] = MPI_REQUEST_NULL;
		recv[ne] = MPI_REQUEST_NULL;
	}

	// determine if the global domain is periodic
	bool periodic = true;
	mpiInfo->periodic = true;
	periodicAll[mpiInfo->mpiRank] = 1; //true
	for(int i=1;i<nDims+1;i++){
		if(bnd[i] != PERIODIC || bnd[i+nDims+1] != PERIODIC) periodicAll[mpiInfo->mpiRank] = 0;
		//printf("i+nDims = %i, bnd[i] = %d, bnd[2*i] = %d \n",i+nDims+2,bnd[i],bnd[i+nDims+1]);

	}



	//if(periodic == true){
		//MPI_Bcast(&periodicAll[mpiInfo->mpiRank], 1, MPI_INT, mpiInfo->mpiRank, MPI_COMM_WORLD);
		MPI_Allgather(&periodicAll[mpiInfo->mpiRank], 1, MPI_INT,
                  periodicAll, 1, MPI_INT,  MPI_COMM_WORLD);
	//}
	for(int i = 0;i<mpiSize;i++){
		if(periodicAll[i] != 1) mpiInfo->periodic = false;
		// if(mpiInfo->mpiRank == 0){
		// 	printf("periodicAll[i] = %i \n",periodicAll[i]);
		// }
	}

	//printf("periodic = %d \n",periodicAll[mpiInfo->mpiRank]);
	//printf("periodic = %d \n",mpiInfo->periodic);

	free(periodicAll);
	//exit(0);
	//mpiInfo->periodic = periodic;
	mpiInfo->send = send;
	mpiInfo->recv = recv;
	mpiInfo->nNeighbors = nNeighbors;
	mpiInfo->migrants = migrants; // depricated ?
	mpiInfo->migrantsDummy = migrantsDummy; // depricated ?
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
	//free(mpiInfo->periodic);
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
