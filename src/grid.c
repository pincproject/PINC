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
/*******TEMPORARY*************/
#include "test.h"
/*******TEMPORARY*************/



/******************************************************************************
* DECLARATIONS,
*****************************************************************************/
/**
 * @brief Returns the ND-index of this MPI node in the global reference frame
 * @param	ini		input settings
 * @return	The N-dimensional index of this MPI node
 */
static int *getSubdomain(const dictionary *ini);

/**
 * @brief Extracts a (dim-1) dimensional slice of grid values.
 * @param double *slice 	empty array of size slice
 * @param GridQuantity *gridQuantity
 * @param int d 					perpendicular direction to slice
 * @param int offset 			offset of slice
 * @return double *slice
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
	 getSlice(slice, gridQuantity, 0, 1);
 * @code
 * After running this the slice array consists of
 * slice = \f( [1, 6, 11, 16] \f)
 *
 * @see setSlice
 * @see swapHalo
 **/

void getSlice(double *slice, const GridQuantity *gridQuantity, int d, int offset);

/**
 * @brief places a (dim-1) dimensional slice onto a selected slice on the  grid.
 * @param const double *slice
 * @param GridQuantity *gridQuantity
 * @param const int d 				(perpendicular direction to slice)
 * @param const int offset 		(offset of slice)
 * @return Gridquantity *gridQuantity
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
 *	setSlice(slice, gridQuantity, 0, 2);
 *
 *	111111
 *	222222
 *	111111
 *	111111
 * @endcode
 *
 * @see setSlice
 * @see swapHalo
 **/
void setSlice(const double *slice, GridQuantity *gridQuantity, int d, int offset);

/**********************************************************
 *	Local functions
 *********************************************************/

double *getSliceInner(double *nextGhost, const double **valp, const long int *mul,
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

void getSlice(double *slice, const GridQuantity *gridQuantity, int d, int offset){

	Grid *grid = gridQuantity->grid;
	int nDims = grid->nDims;
	long int *nGPointsProd = grid->nGPointsProd;
	int *nGPoints = grid->nGPoints;

	const double *val = gridQuantity->val;

	val += offset*nGPointsProd[d];

	getSliceInner(slice, &val, &nGPointsProd[nDims-1], &nGPoints[nDims-1],nGPointsProd[d]);
}

const double *setSliceInner(const double *nextGhost, double **valp, const long int *mul,
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

void setSlice(const double *slice, GridQuantity *gridQuantity, int d, int offset){

	Grid *grid = gridQuantity->grid;
	int nDims = grid->nDims;
	long int *nGPointsProd = grid->nGPointsProd;
	int *nGPoints = grid->nGPoints;

	double *val = gridQuantity->val;

	val += offset*nGPointsProd[d];
	setSliceInner(slice, &val, &nGPointsProd[nDims-1], &nGPoints[nDims-1],nGPointsProd[d]);
}



/******************************************************************************
* DEFINITIONS
*****************************************************************************/

static int *getSubdomain(const dictionary *ini){

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Get ini info
	int nDims;
	int *nSubdomains = (int*)iniGetIntArr(ini,"grid:nSubdomains",&nDims);

	// Sanity check
	int totalNSubdomains = intArrProd(nSubdomains,nDims);
	if(totalNSubdomains!=mpiSize)
		msg(ERROR|ONCE,"The product of grid:nSubdomains does not match the number of MPI processes");

	// Determine subdomain of this MPI node
	int *subdomain = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++){

		subdomain[d] = mpiRank % nSubdomains[d];
		mpiRank /= nSubdomains[d];

	}

	free(nSubdomains);
	return subdomain;

}


void swapHalo(dictionary *ini, GridQuantity *gridQuantity){

	//Load
	Grid *grid = gridQuantity->grid;
	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;
	// long int *nGPointsProd = grid->nGPointsProd;

	int sliceDim = 1;
	int offset = 0;
	int offset2= nGPoints[sliceDim] -1;
	int nSlicePoints = 1;
	for(int d = 0; d < nDims ; d++) nSlicePoints *=nGPoints[d];
	nSlicePoints /= nGPoints[sliceDim];

	double *slice = malloc(nSlicePoints*sizeof(double));
	double *slice2 = malloc(nSlicePoints*sizeof(double));

	getSlice(slice, gridQuantity, sliceDim, offset);
	getSlice(slice2, gridQuantity, sliceDim, offset2);
	setSlice(slice2, gridQuantity, sliceDim, offset);
	setSlice(slice, gridQuantity, sliceDim, offset2);

	// msg(STATUS, "**Slice obtained**");
	// for(int p = 0; p < nSlicePoints; p++) msg(STATUS, "%f", slice[p]);
	

	return;
}


Grid *allocGrid(const dictionary *ini){

	//Sanity check
	iniAssertEqualNElements(ini, 2,"grid:nTGPoints", "grid:dr");

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims, nBoundaries;
	int *nTGPoints = iniGetIntArr(ini, "grid:nTGPoints", &nDims);
	int *nGhosts = iniGetIntArr(ini, "grid:nGhosts", &nBoundaries);
	double *dr = iniGetDoubleArr(ini, "grid:dr", &nDims);

	//More sanity check
	if(nBoundaries != 2*nDims){
		msg(ERROR|ONCE, "Need ghost cells depth for all the boundaries: 2*nDims");
	}

	// Calculate the number of grid points (True points + ghost points)
	int *nGPoints = malloc(nDims *sizeof(int));

	for(int d = 0 ; d < nDims; d ++){
		nGPoints[d] = nTGPoints[d];
		nGPoints[d] += nGhosts[d];
		nGPoints[d] += nGhosts[nDims + d];
	}

	//Cumulative products
	long int *nGPointsProd = longIntArrCumProd(nGPoints,nDims);

	//Free temporary variables
	free(nTGPoints);

	/* Store in Grid */
	Grid *grid = malloc(sizeof(Grid));

	grid->nDims = nDims;
	grid->nGPoints = nGPoints;
	grid->nGPointsProd = nGPointsProd;
	grid->nGhosts = nGhosts;
	grid->dr = dr;

	return grid;
}

MpiInfo *allocMpiInfo(const dictionary *ini){
	//Sanity check
	iniAssertEqualNElements(ini, 2,"grid:nSubdomains","grid:nTGPoints");

	// Get MPI info
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Load data from ini
	int nDims;
	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", &nDims);
	int *nTGPoints = iniGetIntArr(ini, "grid:nTGPoints", &nDims);
	int *nSubdomainsProd = intArrCumProd(nSubdomains,nDims);

	//Position of the subdomain in the total domain
	int *subdomain = getSubdomain(ini);
	int *offset = malloc(nDims*sizeof(int));
	double *posToSubdomain = malloc(nDims*sizeof(double));

	for(int d = 0; d < nDims; d++){
		offset[d] = subdomain[d]*nTGPoints[d];
		posToSubdomain[d] = (double)1/nTGPoints[d];
	}

    /* Store in Grid */
    MpiInfo *mpiInfo = malloc(sizeof(MpiInfo));

	mpiInfo->subdomain = subdomain;
	mpiInfo->nSubdomains = nSubdomains;
	mpiInfo->nSubdomainsProd = nSubdomainsProd;
	mpiInfo->offset = offset;
	mpiInfo->posToSubdomain = posToSubdomain;
	mpiInfo->mpiSize = mpiSize;
	mpiInfo->mpiRank = mpiRank;

	free(nTGPoints);

    return mpiInfo;
}

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues){

	//Load data
	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;
	int *nGhosts = grid->nGhosts;
	long int *nGPointsProd = grid->nGPointsProd;
	int nTotPoints = 1;			//#Grid points in all dimensions
	int nGhostPoints = 0;		//#Total ghost points


	//Total grid points N^d
	for(int d = 0; d < nDims; d++){
		nTotPoints *= nGPoints[d];
	}
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g+nDims])*nGPoints[g];
	}

	//Memory for values
	double *val = malloc(nTotPoints*nValues*sizeof(*val));
	double *slice = malloc(nGPointsProd[nDims-1]*nValues*sizeof(*slice));

	/* Store in gridQuantity */
	GridQuantity *gridQuantity = malloc(sizeof(*gridQuantity));

	gridQuantity->grid = grid;
	gridQuantity->nValues = nValues;
	gridQuantity->val = val;
	gridQuantity->h5 = 0;	// Must be activated separately
	gridQuantity->slice = slice;

	return gridQuantity;
}

void freeMpiInfo(MpiInfo *mpiInfo){

	free(mpiInfo->subdomain);
	free(mpiInfo->nSubdomains);
	free(mpiInfo->nSubdomainsProd);
	free(mpiInfo->offset);
	free(mpiInfo->posToSubdomain);
	free(mpiInfo);

}

void freeGrid(Grid *grid){

	free(grid->nGPoints);
	free(grid->nGPointsProd);
	free(grid->nGhosts);
	free(grid->dr);
	free(grid);

}

void freeGridQuantity(GridQuantity *gridQuantity){

	free(gridQuantity->val);
	free(gridQuantity->slice);
	free(gridQuantity);

}
