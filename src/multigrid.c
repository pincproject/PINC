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
#include "pinc.h"



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
    multigrid->preSmooth = &gaussSeidel;

    /*
     * Test area pointer function
    */
    /*int strElements;
    char **preSmoothName = iniparser_getstring(ini, "multigrid:");*/
  	
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

void multigridParseDump(dictionary *ini, Multigrid *multigrid){

	fMsg(ini,"parsedump", "Multigrids: \n");


	return;
}


