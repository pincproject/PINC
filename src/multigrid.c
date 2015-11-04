/**
 * @file		multigrid.c
 * @author		Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Poisson Solver.
 * @date		26.10.15
 *
 *
 * Functions dealing with the initialisation and destruction of multigrid structures and
 * a multigrid solver containing restriction, prolongation operatorors and smoothers 
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "multigrid.h"
#include "grid.c"
#include "pinc.h"


Grid *allocGrid(const dictionary *ini){
	/**
	* This function initalises the grid struct, and sets all the parameters from
	* the input file
	* Input: 	ini: 	initilization structure
	* Output:	Grid *grid
	*/

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
	double *dr = iniGetDoubleArr(ini, "grid:dr", &nDims);

	//More sanity check
	if(nBoundaries != 2*nDims){
		msg(ERROR, "Need ghost cells depth for all the boundaries");
	}
	
	// Calculate the number of grid points (True points + ghost points)
	int *nGPoints = malloc(nDims *sizeof(int));

	for(int i = 0 ; i < nDims; i ++){
		nGPoints[i] = nTGPoints[i];
		nGPoints[i] += nGhosts[i];
		nGPoints[i] += nGhosts[nDims + i];
	}

	//Cumulative products
	int *nGPointsProd = malloc (nDims*sizeof(int));
	int tempPoints = 1;
	for(int i = 0; i < nDims; i++){
		nGPointsProd[i] = tempPoints*nGPoints[i];
		tempPoints = nGPointsProd[i];
	}

	//Position of the subdomain in the total domain
	int *node = getNode(ini);
	int *offset = malloc(nDims*sizeof(int));
	double *posToNode = malloc(nDims*sizeof(double));

	for(int i = 0; i < nDims; i++){
		offset[i] = node[i]*nTGPoints[i];
		posToNode[i] = (double) offset[i] * dr[i];
	}

/*	msg(WARNING, "node[%d, %d, %d] has offset[] = [%d, %d, %d] and position[] = [%f, %f, %f]", \
		node[0], node[1],node[2],offset[0], offset[1], offset[2], posToNode[0], posToNode[1], posToNode[2]);
*/
/*	printf("True Grid Points = %d \n", nTGPoints[0]);
	printf("nGhosts = %d \n", nGhosts[0]);
	printf("nBoundaries = %d \n", nBoundaries);
	printf("Grid points: = %d \n", nGPoints[0]);
	printf("nGPoints = %d \n", nGPoints[0]);

	printf("I process %d, am node [%d, %d, %d] \n", rank, node[0], node[1],node[2]);

	printf("nGPointsProd[nDims] = %d \n", nGPointsProd[nDims-1]);

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


	//Free all the variables used
/*	free(nDims);
*/


    return grid;
}



void freeGrid(Grid *grid){
	/**
	 * This functions frees the memory stored in the grid
	*/

	/*free(grid->nValues);*/

	return;
}

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid){

	GridQuantity *gridQuantity = malloc(sizeof(GridQuantity));

	return gridQuantity;
}


Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity){
	/**
	* This function initalises the MultiGrid struct,
	* Input: 	grid: 	pointer to finegrid struct
	*			nNodes:	...
	*
	* This function takes a the fine grid, and creates several grids,
	* each with half the number of points and stores them as an array 
	* of grids
	*/

	/*Multigrid *multigrid = malloc (sizeof(Multigrid));
*/
	/*int i, dimension;
	int nNodes[2] = {4, 4};
	int nDim = fineGrid->nDim;
	int nValues = fineGrid->nValues;

	printf("Hello from multigrid \n");

	multigrid->grids = (Grid *)malloc(nLevels * sizeof(Grid));

	multigrid->grids[0] = *fineGrid;



	printf("%d \n", fineGrid->nNodes[0]);

	
	for(i = 1; i < nLevels; i++){
		for(dimension = 0; dimension < nDim; dimension++){
			nNodes[dimension] *= 0.5;
		}
    	struct Grid grid;
    	gridInit(&grid, nNodes, nDim, nValues);
    	multigrid->grids[i] = grid;
	}*/

    Multigrid *multigrid = malloc(sizeof(Multigrid));

    /*
     * Test area pointer function
    */
    /*int strElements;
    char **preSmoothName = iniGetStrArr(ini,"algorithms:preSmooth", &strElements);
  	
    printf("%s  \n",preSmoothName[0]);*/
    preSmooth = &gaussSeidel;

  	preSmooth();

	return multigrid;
}

void jacobian(void){
	printf("Hello from Jacobian \n");
	return;
}

void gaussSeidel(void){

	printf("Hello from Gauss Seidel\n");
	return;
}
