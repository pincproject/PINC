#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "multigrid.h"
#include "pinc.h"


Grid *allocGrid(const dictionary *ini){
	/**
	* This function initalises the grid struct, and sets all the parameters from
	* the input file
	* Input: 	ini: 	initilization structure
	* Output:	Grid *grid
	*/

	//Sanity check
	//TBD

	// Get MPI info
	int nNodes, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nNodes);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data from ini
	int nDims = iniGetNElements(ini, "grid:nCells");


/*    grid->nDim = nDim;
    grid->nValues = nValues;

    grid->nNodes = (int *)malloc(nDim * sizeof(int));
    grid->nNodes = nNodes;

    grid->nNodesProd = (int *)malloc(nDim * sizeof(int));

    //Next part needs to be valid for different dimensions
    grid->values = (int *)malloc(nNodes[0]*nNodes[1] * sizeof(int));*/

    /* Store in Grid */
    Grid *grid = malloc(sizeof(Grid));

	grid->nDims = nDims;

    return grid;
}



void freeGrid(Grid *grid){
	/**
	 * This functions frees the memory stored in the grid
	*/

	/*free(grid->nValues);*/

	return;
}

Multigrid *allocMultigrid(const dictionary *ini){
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
    int strElements;
    char **preSmoothName = iniGetStrArr(ini,"algorithms:preSmooth", &strElements);
  	
    printf("%s  \n",preSmoothName[0]);
    preSmooth = &jacobian;

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
