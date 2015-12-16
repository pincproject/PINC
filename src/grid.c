/**
 * @file		grid.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Grid-struct handling.
 * @date		30.10.15
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#include "pinc.h"
#include <mpi.h>
#include <math.h>
#include <hdf5.h>


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
 * @brief Extracts a (dim-1) dimensional slice of grid values.
 * @param	slice 		Return array
 * @param	grid		Grid
 * @param	d			Perpendicular direction to slice
 * @param	offset 		Offset of slice
 * @return				Void
 *
 * This function gets extracts a slice from a N dimensional grid. The integer d
 * decides in which direction the slice is perpendicular to, and the offset decides
 * which which slice it picks out. It needs a preallocated slice array where
 * the extracted slice will be stored.
 *
 * 2D example: Here we have a 5x4 grid and we want to extract a slice corresponding
 * to the second row, where x (d=0) is a constant 1.
 *
 * @code
 * 15   16   17   18   19
 *
 * 10   11   12   13   14
 *
 *  5    6    7    8    9

 *  0    1    2    3    4
 * @endcode
 *
 * @code
	 getSlice(slice, grid, 0, 1);
 * @endcode
 * After running this the slice array consists of
 * slice = \f( [1, 6, 11, 16] \f)
 *
 * @see setSlice
 * @see swapHalo
 **/

void getSlice(double *slice, const Grid *grid, int d, int offset);

/**
 * @brief places a (dim-1) dimensional slice onto a selected slice on the grid.
 * @param	slice		Slice containing a layer of values
 * @param	grid		Grid
 * @param	d 			Perpendicular direction to slice
 * @param	offset 		Offset of slice
 * @return				Void
 *
 * This function places a a slice on a grid. If we have a slice and want to
 * insert it onto a grid this function is used.
 *
 * Example: We have a 1D slice consisting of 6 2s and want to insert it onto
 * the third row, from the bottom.
 * @code
 *	111111
 *	111111
 *	111111
 *	111111
 *
 *	setSlice(slice, grid, 0, 2);
 *
 *	111111
 *	222222
 *	111111
 *	111111
 * @endcode
 *
 * @see setSlice
 * @see swapHalo
 */
void setSlice(const double *slice, Grid *grid, int d, int offset);

/**
 * @brief Adds a slice to a slice in a Grid
 * @param	slice		Slice of values to add into grid
 * @param	grid		Grid
 * @param	d			Perpendicular direction to slice grid
 * @param	offset		Offset of slice in grid
 * @return				void
 *
 * Similar to setSlice() but adds slice to existing values rather than replacing
 * them.
 */
void addSlice(const double *slice, Grid *grid, int d, int offset);

/**
 * @brief Gets, sends, recieves and sets a slice, using MPI
 * @param nSlicePoints		Length of the slice array
 * @param offsetTake 		Offset of slice to get
 * @param offsetSet 		Offset of where to set slice
 * @param d					Dimension
 * @param reciever 			mpiRank of subdomain to recieve slice
 * @param sender 			mpiRank of subdomain,
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
 * @see swapHalo
 */
void getSendRecvSetSlice(const int nSlicePoints, const int offsetTake,
					const int offsetPlace, const int d, const int reciever,
					const int sender, const int mpiRank, Grid *grid);

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
				nextGhost = setSliceInner(nextGhost, valp, mul-1,points-1,finalMul);
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

void getSendRecvSetSlice(const int nSlicePoints, const int offsetTake,
					const int offsetPlace, const int d, const int reciever,
					const int sender, const int mpiRank, Grid *grid){

	double *slice = grid->slice;

	getSlice(slice, grid, d, offsetTake);
	//Send and recieve (Need to check out if using one of the more sophisticated send
	//functinos to MPI could be used)
	MPI_Send(slice, nSlicePoints, MPI_DOUBLE, reciever, mpiRank, MPI_COMM_WORLD);
	MPI_Recv(slice, nSlicePoints, MPI_DOUBLE, sender, sender, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	setSlice(slice, grid, d, offsetPlace);
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
	int totalNSubdomains = intArrProd(nSubdomains,nDims);
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

// void swapHalo(dictionary *ini, GridQuantity *gridQuantity, MpiInfo *mpiInfo, int d){
//
// 	//Load MpiInfo
// 	int mpiRank = mpiInfo->mpiRank;
// 	int *subdomain = mpiInfo->subdomain;
// 	int *nSubdomains = mpiInfo->nSubdomains;
// 	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;
//
// 	//Load
// 	Grid *grid = gridQuantity->grid;
// 	int nDims = grid->nDims;
// 	int *nGPoints = grid->nGPoints;
//
// 	//Local temporary variables
// 	int reciever, sender;
// 	int nSlicePoints = 1;
// 	int offsetTake, offsetPlace;
// 	for(int d = 0; d < nDims ; d++) nSlicePoints *=nGPoints[d];
// 	nSlicePoints *= 1./nGPoints[d];
//
// 	/*****************************************************
// 	 *			Sending and recieving upper
// 	******************************************************/
// 	offsetTake = nGPoints[d]-1;
// 	offsetPlace = 0;
// 	reciever = (mpiRank + nSubdomainsProd[d]);
// 	sender = (mpiRank - nSubdomainsProd[d]);
//
// 	//Here we need an implementation of the boundary conditions, I will probably put in a function called boundaryCond(...)
// 	//here, or alternatively deal with the boundary some other place
// 	if(subdomain[d] == nSubdomains[d] - 1) reciever -= 2*nSubdomainsProd[d];
// 	if(subdomain[d] == 0) sender += 2*nSubdomainsProd[d];
//
// 	getSendRecvSetSlice(nSlicePoints, offsetTake, offsetPlace, d, reciever,
// 					sender, mpiRank, gridQuantity);
//
// 	/*****************************************************
// 	 *			Sending and recieving lower
// 	******************************************************/
// 	offsetTake = 1;
// 	offsetPlace = nGPoints[d]-1;
// 	reciever = (mpiRank - nSubdomainsProd[d]);
// 	sender = (mpiRank + nSubdomainsProd[d]);
//
// 	//Boundary
// 	if(subdomain[d] == nSubdomains[d] - 1) sender -= 2*nSubdomainsProd[d];
// 	if(subdomain[d] == 0) reciever += 2*nSubdomainsProd[d];
//
// 	getSendRecvSetSlice(nSlicePoints, offsetTake, offsetPlace, d, reciever,
// 					sender, mpiRank, gridQuantity);
//
//
// 	return;
// }

Grid *gAlloc(const dictionary *ini, int nValues){

	//Sanity check
	iniAssertEqualNElements(ini, 2,"grid:trueSize", "grid:stepSize");

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims, nBoundaries;
	int *trueSizeTemp = iniGetIntArr(ini, "grid:trueSize", &nDims);
	int *nGhostLayersTemp = iniGetIntArr(ini, "grid:nGhostLayers", &nBoundaries);
	double *stepSizeTemp = iniGetDoubleArr(ini, "grid:stepSize", &nDims);

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
	long int *sizeProd = longIntArrCumProd(size,rank);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=0;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<nDims;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//Memory for values and a slice
	double *val = malloc(sizeProd[rank]*sizeof(*val));
	double *slice = malloc(nSliceMax*sizeof(*slice));

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
	grid->slice = slice;

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
	int nDims;
	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", &nDims);
	int *trueSize = iniGetIntArr(ini, "grid:trueSize", &nDims);
	int *nSubdomainsProd = intArrCumProd(nSubdomains,nDims);

	//Position of the subdomain in the total domain
	int *subdomain = getSubdomain(ini);
	int *offset = malloc(nDims*sizeof(*offset));
	double *posToSubdomain = malloc(nDims*sizeof(*posToSubdomain));

	for(int d = 0; d < nDims; d++){
		offset[d] = subdomain[d]*trueSize[d];
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
	free(grid->slice);
	free(grid);

}

void gValDebug(Grid *grid, const MpiInfo *mpiInfo){

	int mpiRank = mpiInfo->mpiRank;
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	double *v = grid->val;

	for(long int p=0;p<sizeProd[rank];p++){

		v[p] = 0;
		long int temp = p;

		for(int d=0;d<rank;d++){
			v[p] += (temp%size[d])*pow(10,d-1) + mpiRank*1000;
			temp/=size[d];
		}

	}
}

void gWriteH5(const Grid *grid, const MpiInfo *mpiInfo, double n){

	hid_t fileSpace = grid->h5FileSpace;
	hid_t memSpace = grid->h5MemSpace;
	hid_t file = grid->h5;
	double *val = grid->val;

	/*
	 * STORE DATA COLLECTIVELY
	 */
	hid_t pList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);

	char name[64];
	sprintf(name,"/n=%.1f",n);

	hid_t dataset = H5Dcreate(file,name,H5T_IEEE_F64LE,fileSpace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);
	H5Dclose(dataset);

	H5Pclose(pList);
}

void gCloseH5(Grid *grid){
	H5Sclose(grid->h5MemSpace);
	H5Sclose(grid->h5FileSpace);
	H5Fclose(grid->h5);
}

void gCreateH5(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo,
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

	hid_t file = createH5File(ini,fName,"grid");

	/*
	 * CREATE ATTRIBUTES
	 */

	double *debye = malloc(nDims*sizeof(*debye));
	debye[0] = iniparser_getdouble((dictionary *)ini,"grid:debye",0);
	for(int d=1;d<nDims;d++) debye[d]=debye[0];

	hsize_t attrSize;
    hid_t attrSpace;
    hid_t attribute;

	attrSize = (hsize_t)nDims;
	attrSpace = H5Screate_simple(1,&attrSize,NULL);

	attribute = H5Acreate(file, "Axis denormalization factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, grid->stepSize);
    H5Aclose(attribute);

	attribute = H5Acreate(file, "Axis dimensionalizing factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, debye);
    H5Aclose(attribute);

    H5Sclose(attrSpace);
	free(debye);

	attrSize = (hsize_t)size[0];
	attrSpace = H5Screate_simple(1,&attrSize,NULL);

	attribute = H5Acreate(file, "Quantity denormalization factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, denorm);
	H5Aclose(attribute);

	attribute = H5Acreate(file, "Quantity dimensionalizing factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, dimen);
	H5Aclose(attribute);

	H5Sclose(attrSpace);

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
		fileDims[d]		= (hsize_t)trueSize[rank-d-1]*nSubdomains[rank-d-1];
		fileOffset[d]	= (hsize_t)trueSize[rank-d-1]*subdomain[rank-d-1];
	}

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

void gMulDouble(Grid *grid, double num){

	int rank = grid->rank;
	long int nElements = grid->sizeProd[rank];
	for(long int p=0;p<nElements;p++) grid->val[p] *= num;
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
	double *q = iniGetDoubleArr(ini,"population:q",&nSpecies);
	double *m = iniGetDoubleArr(ini,"population:m",&nSpecies);
	double timeStep = iniparser_getdouble((dictionary *)ini,"time:timeStep",0.0);
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",&nDims);
	gMulDouble(E,pow(timeStep,2)*(q[0]/m[0]));
	for(int p=0;p<E->sizeProd[E->rank];p++){
		E->val[p] /= stepSize[p%E->size[0]];
	}

}
