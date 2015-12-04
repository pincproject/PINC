/**
 * @file		multigrid.c
 * @author		Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Poisson Solver, multigrid.
 * @date		26.10.15
 *
 *
 * Functions dealing with the initialisation and destruction of multigrid structures and
 * a multigrid solver containing restriction, prolongation operatorors and smoothers
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "pinc.h"
#include "multigrid.h"

/******************************************************************************
 * 				Local functions
 *****************************************************************************/
void setSolvers(const dictionary *ini, Multigrid *multigrid){

	char *preSmoothName = iniparser_getstring((dictionary*)ini, "algorithms:preSmooth", "\0");
    char *postSmoothName = iniparser_getstring((dictionary*)ini, "algorithms:postSmooth", "\0");
    char *coarseSolverName = iniparser_getstring((dictionary*)ini, "algorithms:coarseSolv", "\0");

    if(!strcmp(preSmoothName,"gaussSeidel")){
    	multigrid->preSmooth = &gaussSeidel;
    } else {
    	msg(ERROR, "No Presmoothing algorithm specified");
    }

    if(!strcmp(postSmoothName,"gaussSeidel")){
    	multigrid->postSmooth = &gaussSeidel;
    } else {
    	msg(ERROR, "No Postsmoothing algorithm specified");
    }

    if(!strcmp(coarseSolverName,"gaussSeidel")){
    	multigrid->coarseSolv = &gaussSeidel;
    } else {
    	msg(ERROR, "No Coarse Grid Solver algorithm specified");
    }

    //Free!

//    printf("If this is the last I say, I'm bad and segfaults at freeing stuff \n");
    // free(coarseSolverName);
    // free(postSmoothName);
    // free(preSmoothName);

	return;
}

GridQuantity **allocSubGrids(const dictionary *ini, GridQuantity *gridQuantity,
							const int nLevels){

	GridQuantity **gridQuantities = malloc(nLevels * sizeof(GridQuantity));
	int nDims;
	int nBoundaries;

	//Set first grid to point to fine grid
	gridQuantities[0] = gridQuantity;

	int *nTGPoints = iniGetIntArr(ini, "grid:nTGPoints", &nDims);
	int *nGhosts = iniGetIntArr(ini, "grid:nGhosts", &nBoundaries);
	int nValues = gridQuantity->nValues;
	for(int q = 1; q < nLevels; q++){

		//The subgrid needs half the grid points
		for(int d = 0; d < nDims; d++)	nTGPoints[d] /= 2;

		// Calculate the number of grid points (True points + ghost points)
		int *nGPoints = malloc(nDims *sizeof(int));

		for(int d = 0 ; d < nDims; d ++){
			nGPoints[d] = nTGPoints[d];
			nGPoints[d] += nGhosts[d];
			nGPoints[d] += nGhosts[nDims + d];
		}

		//Cumulative products
		long int *nGPointsProd = malloc (nDims*sizeof(int));

		int tempPoints = 1;
		for(int d = 0; d < nDims; d++){
			nGPointsProd[d] = tempPoints*nGPoints[d];
			tempPoints = nGPointsProd[d];
		}

		//Assign to grid
		Grid *grid = malloc(sizeof(Grid));

		grid->nDims = nDims;
		grid->nGPoints = nGPoints;
		grid->nGPointsProd = nGPointsProd;
		grid->nGhosts = nGhosts;
		//The rest picks the information the fine grid
		GridQuantity *gridQuantity = allocGridQuantity(ini, grid, nValues);

		gridQuantities[q] = gridQuantity;

	}


	return gridQuantities;
}


/*************************************************
 *		DEFINITIONS
 ************************************************/

Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity){

	//Multigrid
	int nLevels = iniparser_getint((dictionary *) ini, "multigrid:mgLevels", 0);
	int nCycles = iniparser_getint((dictionary *) ini, "multigrid:mgCycles", 0);

	//Sanity checks
	if(!nLevels){
		msg(ERROR, "Multi Grid levels is 0, direct solver not implemented yet \n");
	}

	if(!nCycles){
		msg(ERROR, "MG cycles is 0 \n");
	}


	// // Sanity check (true grid points need to be a multiple of 2^(multigrid levels)
	// for(int d = 0; d < nDims; d++){
	// 	if(nTGPoints[d] % (int) pow(2, nLevels)){ //Sloppy and wrong
	// 		msg(ERROR, "The number of True Grid Points needs to be a multiple of 2^nLevels");
	// 	}
	// }

	GridQuantity **gridQuantities = allocSubGrids(ini, gridQuantity, nLevels);

	//Store in multigrid struct
    Multigrid *multigrid = malloc(sizeof(Multigrid));

    multigrid->nLevels = nLevels;
    multigrid->nCycles = nCycles;
    multigrid->gridQuantities = gridQuantities;

    //Setting the algorithms to be used, pointer functions
	setSolvers(ini, multigrid);

  	return multigrid;

}

void freeMultigrid(Multigrid *multigrid){

	GridQuantity **gridQuantities = multigrid->gridQuantities;
	int nLevels = multigrid->nLevels;

	for(int n = 1; n < nLevels; n++)
	{
		freeGrid(gridQuantities[n]->grid);
		freeGridQuantity(gridQuantities[n]);
	}
	free(multigrid);

	return;
}


void jacobian(void){
	printf("Hello from Jacobian \n");
	return;
}

void gaussSeidel(void){

	printf("Hello from Gauss Seidel\n");
	return;
}
