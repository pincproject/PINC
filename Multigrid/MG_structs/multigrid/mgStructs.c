#include "mgStructs.h"
#include <stdio.h>
#include <stdlib.h>


void gridInit(Grid *grid, int *nNodes, int nDim,int nValues){
	/**
	* This function initalises the grid struct,
	* Input: 	grid: 	pointer to grid struct
	*			nNodes:	...
	*/

    grid->nDim = nDim;
    grid->nValues = nValues;

    grid->nNodes = (int *)malloc(nDim * sizeof(int));
    grid->nNodes = nNodes;

    grid->nNodesProd = (int *)malloc(nDim * sizeof(int));

    //Next part needs to be valid for different dimensions
    grid->values = (int *)malloc(nNodes[0]*nNodes[1] * sizeof(int));

    return;
}

void multigridInit(Multigrid *multigrid, Grid *fineGrid, int nLevels){
	/**
	* This function initalises the MultiGrid struct,
	* Input: 	grid: 	pointer to finegrid struct
	*			nNodes:	...
	*
	* This function takes a the fine grid, and creates several grids,
	* each with half the number of points and stores them as an array 
	* of grids
	*/

	int i, dimension;
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
	}




	return;
}
