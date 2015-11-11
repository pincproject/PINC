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
#include "multigrid.h"
#include "pinc.h"



Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity){
	

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//Load data
	int nDims, nBoundaries;
	
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

	//Grid
	int *nTGPoints = iniGetIntArr(ini, "grid:nTGPoints", &nDims);
	int *nGhosts = iniGetIntArr(ini, "grid:nGhosts", &nBoundaries);
	int nValues = gridQuantity->nValues;

	// Sanity check (true grid points need to be a multiple of 2^(multigrid levels)
	for(int d = 0; d < nDims; d++){
		if(nTGPoints[d] % (int) pow(2, nLevels)){ //Sloppy and wrong
			msg(ERROR, "The number of True Grid Points needs to be a multiple of 2^nLevels");
		}
	}

	//Declare gridQuantity array
	GridQuantity **gridQuantities = malloc(nLevels * sizeof(GridQuantity));

	//Set first grid to point to fine grid
	gridQuantities[0] = gridQuantity;

	/*
	 * Here the true grid points is divided by 2, and then a new grid + gridquantity
	 * allocated
	 */

	//Making the smaller grids (uses a lot of the logic in allocGrid and allocGridQuantity)
	//Should think of a better solution here, but grid needs to have different values for each subgrid
	
	for(int i = 1; i < nLevels; i++){ 
	//Did we agree on a standard index when we aren't dealing with dim, spatial or particle ind?
	//Maybe q would be okay, if that isn't used in population struct

			//The subgrid needs half the grid points
			for(int d = 0; d < nDims; d++){
				nTGPoints[d] /= 2;
			}

			// Calculate the number of grid points (True points + ghost points)
			int *nGPoints = malloc(nDims *sizeof(int));

			for(int d = 0 ; d < nDims; d ++){
				nGPoints[d] = nTGPoints[d];
				nGPoints[d] += nGhosts[d];
				nGPoints[d] += nGhosts[nDims + d];
			}

			//Cumulative products
			int *nGPointsProd = malloc (nDims*sizeof(int));
			
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
			grid->node = gridQuantities[0]->grid->node;		
			grid->nNodes = gridQuantities[0]->grid->nNodes;
			grid->offset = gridQuantities[0]->grid->offset;
			grid->posToNode = gridQuantities[0]->grid->posToNode; 	//Not 100% sure this is correct

			GridQuantity *gridQuantity = allocGridQuantity(ini, grid, nValues);

			gridQuantities[i] = gridQuantity;

	}

	//Store in multigrid struct
    Multigrid *multigrid = malloc(sizeof(Multigrid));

    multigrid->nLevels = nLevels;
    multigrid->nCycles = nCycles;
    multigrid->gridQuantities = gridQuantities;


    //Setting the algorithms to be used, pointer functions
    /*
     *	Will be slightly messy when more options are added later, consider moving to a 
     *	new function/file 
     */
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
			
/*    printf("If this is the last I say, I'm bad and segfaults at freeing stuff \n");
*//*    free(coarseSolverName);
    free(postSmoothName);
    free(preSmoothName);
*/  	
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

void multigridParseDump(dictionary *ini, Multigrid *multigrid){

	fMsg(ini,"parsedump", "Multigrids: \n");

	fMsg(ini, "parsedump", "nCycles \t %d \n", multigrid->nCycles);
	fMsg(ini, "parsedump", "nLevels \t %d \n", multigrid->nLevels);

	/*
	 *         	TEST AREA
	 */
	fMsg(ini, "parsedump", "TEST AREA \n \n");


	fMsg(ini, "parsedump", "Values in the grid in first x values, not sorted: \t");
	for(int i = 0; i < 5; i++){
		fMsg(ini, "parsedump", "%f ", multigrid->gridQuantities[0]->val[0]);
	}

	return;
}


