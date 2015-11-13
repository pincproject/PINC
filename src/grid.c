/**
* @file		grid.c
* @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
*				Gullik Vetvik Killie <gullikvk@gys.uio.no>
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
* DECLARATIONS
*****************************************************************************/
/**
* @brief Returns the ND-index of this MPI node in the global reference frame
* @param	ini		input settings
* @return	The N-dimensional index of this MPI node
*/
static int *getNode(const dictionary *ini);

/******************************************************************************
* DEFINITIONS
*****************************************************************************/

static int *getNode(const dictionary *ini){

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Get ini info
	int nDims;
	int *nNodes = (int*)iniGetIntArr(ini,"grid:nNodes",&nDims);

	// Sanity check
	int totalNNodes = intArrProd(nNodes,nDims);
	if(totalNNodes!=size)
	msg(ERROR,"The product of grid:nNodes does not match the number of MPI processes");

	// Determine node
	int *node = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++){

		node[d] = rank % nNodes[d];
		rank /= nNodes[d];

	}

	free(nNodes);
	return node;

}

Grid *allocGrid(const dictionary *ini){

	//Sanity check
	iniAssertEqualNElements(ini, 3,"grid:nNodes","grid:nTGPoints", "grid:dr");

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data from ini
	int nDims, nBoundaries;
	int *nTGPoints = iniGetIntArr(ini, "grid:nTGPoints", &nDims);
	int *nGhosts = iniGetIntArr(ini, "grid:nGhosts", &nBoundaries);
	int *nNodes = iniGetIntArr(ini, "grid:nNodes", &nDims);
	// double *dr = iniGetDoubleArr(ini, "grid:dr", &nDims);

	//More sanity check
	if(nBoundaries != 2*nDims){
		msg(ERROR, "Need ghost cells depth for all the boundaries: 2*nDims");
	}

	// Calculate the number of grid points (True points + ghost points)
	int *nGPoints = malloc(nDims *sizeof(int));

	for(int d = 0 ; d < nDims; d ++){
		nGPoints[d] = nTGPoints[d];
		nGPoints[d] += nGhosts[d];
		nGPoints[d] += nGhosts[nDims + d];
	}

	//Cumulative products
	int *nGPointsProd = malloc ((nDims+1)*sizeof(int));
	nGPointsProd[0] = 1;
	for(int d = 1; d < nDims+1; d++){
		nGPointsProd[d] = nGPointsProd[d-1]*nGPoints[d-1];
	}

	//Position of the subdomain in the total domain
	int *node = getNode(ini);
	int *offset = malloc(nDims*sizeof(int));
	double *posToNode = malloc(nDims*sizeof(double));

	for(int d = 0; d < nDims; d++){
		offset[d] = node[d]*nTGPoints[d];
		posToNode[d] = (double)1/nTGPoints[d];
	}

	//Free temporary variables
	free(nTGPoints);

	/* Store in Grid */
	Grid *grid = malloc(sizeof(Grid));

	grid->nDims = nDims;
	grid->nGPoints = nGPoints;
	grid->nGPointsProd = nGPointsProd;
	grid->nGhosts = nGhosts;
	grid->node = node;
	grid->nNodes = nNodes;
	grid->offset = offset;
	grid->posToNode = posToNode;

	return grid;
}

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues){

	//Load data
	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;
	int *nGhosts = grid->nGhosts;
	int nTotPoints = 1;			//#Grid points in all dimensions
	int nGhostPoints = 0;		//#Total ghost points


	//Total grid points N^d
	for(int d = 0; d < nDims; d++){
		nTotPoints *= nGPoints[d];
	}
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g +nDims])*nGPoints[g];
	}

	//Memory for values
	double *val = malloc(nTotPoints*nValues*sizeof(double));
	double *halo = malloc(nGhostPoints*nValues*sizeof(double));

	/* Store in gridQuantity */
	GridQuantity *gridQuantity = malloc(sizeof(GridQuantity));

	gridQuantity->grid = grid;
	gridQuantity->nValues = nValues;
	gridQuantity->val = val;
	gridQuantity->halo = halo;

	return gridQuantity;
}


void freeGrid(Grid *grid){

	free(grid->nGPoints);
	free(grid->nGPointsProd);
	free(grid->nGhosts);
	free(grid->node);
	free(grid->nNodes);
	free(grid->offset);
	free(grid->posToNode);

	return;
}

void freeGridQuantity(GridQuantity *gridQuantity){

	free(gridQuantity->val);

	return;
}

double *getHalo(dictionary *ini, GridQuantity *gridQuantity){

	//Picking up data
	Grid *grid = gridQuantity->grid;
	int *nGhosts = grid->nGhosts;
	int *nGPoints = grid->nGPoints;
	int *nGPointsProd = grid->nGPointsProd;
	double *halo = gridQuantity->halo;
	int nDims = grid->nDims;

	//Allocate space for ghost vector
	int nGhostPoints = 0;
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g +nDims])*nGPoints[g];
	}


	//Gather lower edge halo
	int l;
	int w = 0;
	for(int d = 0; d < nDims; d++){
		l = 0;
		for(int g = 0; g < nGPoints[d]; g++){
			halo[w] = gridQuantity->val[l];

			l += nGPointsProd[d];
			w++;
		}
	}

	/*
	*	Note! Look for a clearer way to do higher edge,
	*  and not sure if it works for all dimensions
	*/

	//Gather higher edge halo
	int h;
	int temp = 1;
	for(int d = 0; d < nDims; d++){

		temp *= nGPoints[d];
		h = (nGPointsProd[nDims] - 1) - (temp - nGPointsProd[d]);

		for(int g = 0; g < nGPoints[d]; g++){
			halo[w] = gridQuantity->val[h];

			h += nGPointsProd[d];
			w++;
		}
	}

	return halo;

}

void distributeHalo(dictionary *ini, GridQuantity *gridQuantity){
	msg(STATUS|ONCE, "Distributing halo");

	//Picking up data
	double *halo = gridQuantity->halo;
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

	//Gather lower edge halo
	int l;
	int w = 0;
	for(int d = 0; d < nDims; d++){
		l = 0;
		for(int g = 0; g < nGPoints[d]; g++){
			gridQuantity->val[l] = halo[w];

			l += nGPointsProd[d];
			w++;
		}
	}

	/*
	*	NOTE! Look for a clearer way to do higher edge,
	*  and not sure if it works for all dimensions
	*/

	//Gather higher edge halo
	int h;
	int temp = 1;
	for(int d = 0; d < nDims; d++){

		temp *= nGPoints[d];
		h = (nGPointsProd[nDims] - 1) - (temp - nGPointsProd[d]);

		for(int g = 0; g < nGPoints[d]; g++){
			gridQuantity->val[h] = halo[w];
			if((grid->node[0] == 0) & (grid->node[1] == 0)){
				msg(STATUS|ONCE, "%f", gridQuantity->val[h]);
			}
			h += nGPointsProd[d];
			w++;
		}
	}

	return;
}

void swapHalo(dictionary *ini, GridQuantity *gridQuantity){

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//Load
	Grid *grid = gridQuantity->grid;
	int *node = grid->node;
	int *nNodes = grid->nNodes;
	int nDims = grid->nDims;
	int *nGhosts = grid->nGhosts;
	int *nGPoints = grid->nGPoints;

	int nGhostPoints = 0;
	for(int d = 0; d < nDims; d++){
		nGhostPoints += (nGhosts[d]+nGhosts[d +nDims])*nGPoints[d];
	}

	double *halo = getHalo(ini, gridQuantity);

	//Test stuffs
	if(rank==0){
		dumpHalo(ini, gridQuantity);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==1){
		dumpHalo(ini, gridQuantity);
	}




	if(rank == 0){
		MPI_Send(halo, nGhostPoints, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		MPI_Recv(halo, nGhostPoints, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		}

		if(rank == 1){
			MPI_Send(halo, nGhostPoints, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
			MPI_Recv(halo, nGhostPoints, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);

			}
			MPI_Barrier(MPI_COMM_WORLD);
			//Test stuffs
			if(rank==1){
				fMsg(ini,"parsedump", "\n*****\nSwap\n*****\n");
			}
			if(rank==0){
				dumpHalo(ini, gridQuantity);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank==1){
				dumpHalo(ini, gridQuantity);
			}

			distributeHalo(ini, gridQuantity);


			return;
		}



		void gridParseDump(dictionary *ini, Grid *grid, GridQuantity *gridQuantity){
			/******************************************
			*	Writing information to the parsedump
			*******************************************/
			fMsg(ini,"parsedump", "Grids: \n");

			fMsg(ini,"parsedump", "#Computational Nodes: ");
			for(int d = 0; d < grid->nDims; d++){
				fMsg(ini,"parsedump", "%d " , grid->nNodes[d]);
			}

			fMsg(ini,"parsedump", "\nTotal true grid points: ");
			for(int d = 0; d < grid->nDims; d++){
				fMsg(ini, "parsedump", "%d ", (grid->nGPoints[d]- \
					(grid->nGhosts[d] + grid->nGhosts[grid->nDims +1]))*grid->nNodes[d]);
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
