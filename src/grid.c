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

/******************************************************************************
 * DECLARATIONS
 *****************************************************************************/
/**
 * @brief Returns the ND-index of this MPI node in the global reference frame
 * @param	ini		input settings
 * @return	The N-dimensional index of this MPI node
 */
static int *getSubdomain(const dictionary *ini);

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
	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", &nDims);
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
	int *nGPointsProd = intArrCumProd(nGPoints,nDims);
	int *nSubdomainsProd = intArrCumProd(nSubdomains,nDims);

	//Position of the subdomain in the total domain
	int *subdomain = getSubdomain(ini);
	int *offset = malloc(nDims*sizeof(int));
	double *posToSubdomain = malloc(nDims*sizeof(double));

	for(int d = 0; d < nDims; d++){
		offset[d] = subdomain[d]*nTGPoints[d];
		posToSubdomain[d] = (double)1/nTGPoints[d];
	}

	//Free temporary variables
	free(nTGPoints);

    /* Store in Grid */
    Grid *grid = malloc(sizeof(Grid));

	grid->nDims = nDims;
	grid->nGPoints = nGPoints;
	grid->nGPointsProd = nGPointsProd;
	grid->nGhosts = nGhosts;
	grid->subdomain = subdomain;
	grid->nSubdomains = nSubdomains;
	grid->nSubdomainsProd = nSubdomainsProd;
	grid->offset = offset;
	grid->posToSubdomain = posToSubdomain;
	grid->dr = dr;
	grid->mpiSize = mpiSize;
	grid->mpiRank = mpiRank;

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
	mpiInfo->dr = dr;
	mpiInfo->mpiSize = mpiSize;
	mpiInfo->mpiRank = mpiRank;

    return mpiInfo;
}

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues){

	//Load data
	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;
	int nTotPoints = 1;		//#Grid points in all dimensions

	//Total grid points N^d
	for(int d = 0; d < nDims; d++){
		nTotPoints *= nGPoints[d];
	}

	//Memory for values
	double *val = malloc(nTotPoints*nValues*sizeof(double));

	/* Store in gridQuantity */
	GridQuantity *gridQuantity = malloc(sizeof(GridQuantity));

	gridQuantity->grid = grid;
	gridQuantity->nValues = nValues;
	gridQuantity->val = val;
	gridQuantity->h5 = 0;	// Must be activated separately

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
	free(grid->subdomain);
	free(grid->nSubdomains);
	free(grid->nSubdomainsProd);
	free(grid->offset);
	free(grid->posToSubdomain);
	free(grid->dr);
	free(grid);

	return;
}

void freeGridQuantity(GridQuantity *gridQuantity){

	free(gridQuantity->val);

	return;
}

double *getGhostEdge(dictionary *ini, GridQuantity *gridQuantity){

	//Picking up data
	Grid *grid = gridQuantity->grid;
	int *nGhosts = grid->nGhosts;
	int *nGPoints = grid->nGPoints;
	int *nGPointsProd = grid->nGPointsProd;
	int nDims = grid->nDims;

	//Allocate space for ghost vector
	int nGhostPoints = 0;
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g +nDims])*nGPoints[g];
	}

	double *ghostEdge = malloc(nGhostPoints*sizeof(double));

	for(int w = 0; w < nGhostPoints; w++){
		ghostEdge[w] = 0.;
	}

	//Gather lower edge ghosts
	int l;
	int w = 0;
	for(int d = 0; d < nDims; d++){
		l = 0;
		for(int g = 0; g < nGPoints[d]; g++){
			ghostEdge[w] = gridQuantity->val[l];

			l += nGPointsProd[d];
			w++;
		}
	}

	/*
	 *	NOTE! Look for a clearer way to do higher edge,
	 *  and not sure if it works for all dimensions
	 */

	//Gather higher edge ghosts
	int h;
	int temp = 1;
	for(int d = 0; d < nDims; d++){

		temp *= nGPoints[d];
		h = (nGPointsProd[nDims] - 1) - (temp - nGPointsProd[d]);

		for(int g = 0; g < nGPoints[d]; g++){
			ghostEdge[w] = gridQuantity->val[h];

			h += nGPointsProd[d];
			w++;
		}
	}

	return ghostEdge;

}


void gridParseDump(dictionary *ini, Grid *grid, GridQuantity *gridQuantity){
	/******************************************
	*	Writing information to the parsedump
	*******************************************/
	fMsg(ini,"parsedump", "Grids: \n");

	fMsg(ini,"parsedump", "#Computational Nodes: ");
	for(int d = 0; d < grid->nDims; d++){
		fMsg(ini,"parsedump", "%d " , grid->nSubdomains[d]);
	}

	fMsg(ini,"parsedump", "\nTotal true grid points: ");
	for(int d = 0; d < grid->nDims; d++){
		fMsg(ini, "parsedump", "%d ", (grid->nGPoints[d]- \
			(grid->nGhosts[d] + grid->nGhosts[grid->nDims +1]))*grid->nSubdomains[d]);
	}


	fMsg(ini,"parsedump", "\n \n");


	/*
	 *         	TEST AREA
	 */
	fMsg(ini, "parsedump", "TEST AREA \n \n");
	fMsg(ini, "parsedump", "Values in the grid in first x values, not sorted: \t");
	for(int i = 0; i < 5; i++){
		fMsg(ini, "parsedump", "%f ", gridQuantity->val[0]);
	}

	fMsg(ini,"parsedump", "\n \n");
	return;
}
