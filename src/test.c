/**
 * @file	    test.c
 * @author	    Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC temp tests
 * @date        09.10.15
 *
 * Test Area for PINC model, a place where some small test are made to see that small
 * parts of the code works as suspected. Here in case one of the develpment test is needed
 * later, then it can be useful to have it stored here.
 *
 * 		NOT TO BE INCLUDED IN FINAL PRODUCT
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include "multigrid.h"
#include "test.h"


void testGridAndMGStructs(dictionary *ini, GridQuantity *gridQuantity, Multigrid *multigrid){

	Grid *grid = gridQuantity->grid;

	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}

	msg(STATUS|ONCE, "Performing a manual test of grid structs");

	//Writing information to parsedump (debugging)
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);

	//Changing them from the multigrid struct
	for(int i = 0; i < 5; i++){
		multigrid->gridQuantities[0]->val[i] = 5.;
	}

	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES FROM MULTIGRID \n");

	//Check that both are changed
	//Writing information to parsedump (debugging)
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);


	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES BACK FROM GRIDQUANTITY \n");

	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);

	return;
}

void testBoundarySendRecieve(dictionary *ini, GridQuantity *gridQuantity, Multigrid *multigrid){

	msg(STATUS|ONCE, "Performing a manual test of boundary communication. Check parsedump");

	//Gathering data from grid
	Grid *grid = gridQuantity->grid;
	int nDims = grid->nDims;
	// int *nGPoints = grid->nGPoints;
	long int *nGPointsProd = grid->nGPointsProd;
	// int *node = grid->node;

	msg(STATUS|ONCE, "Total grid points: \t %d, nGPointsProd = [%d , %d]",\
	 nGPointsProd[nDims], nGPointsProd[0], nGPointsProd[1]);


	//Get rank to put as grid values
	//rank 1 grid has all values 1 and so on
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	//Populate grid as rank.
	for(int g = 0; g < nGPointsProd[nDims]; g++){
		gridQuantity->val[g]= (double) rank;
	}

	/*int l = 0;	//Lower edge
	int h = 0;	//Higher edge
	int temp = 1;

	for(int d = 0; d < nDims; d++){

		l = 0;
		temp *= nGPoints[d];
		h = (nGPointsProd[nDims] - 1) - (temp - nGPointsProd[d]);

		for(int g = 0; g < nGPoints[d]; g++){
			gridQuantity->val[l] = (double) b;
			gridQuantity->val[h] = (double) b + 1;

			l += nGPointsProd[d];
			h += nGPointsProd[d];
		}

		b += 2;
	}*/
	for(int r = 0; r < size; r++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == r){
			dumpGrid(ini, gridQuantity);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	swapHalo(ini, gridQuantity);
	for(int r = 0; r < size; r++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == r){
			dumpGrid(ini, gridQuantity);
		}
	}

	return;
}

void testGetSlice(dictionary *ini, GridQuantity *gridQuantity){
	msg(STATUS|ONCE, "Performing manual test of getSlice, see parsefile");

	Grid *grid = gridQuantity->grid;
	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;
	msg(STATUS, "nGPoints = [%d,%d,%d]", nGPoints[0],nGPoints[1],nGPoints[2]);

	int nLayerPoints = nGPoints[0]*nGPoints[1];
	int ind = 0;

	for(int l = 0; l < nDims; l++){
		for(int p = 0; p < nLayerPoints; p++){
			gridQuantity->val[ind] = ind;
			ind++;
		}
	}
	dumpGrid(ini, gridQuantity);

	swapHalo(ini, gridQuantity);

	dumpGrid(ini, gridQuantity);

}

void dumpGrid(dictionary *ini, GridQuantity *gridQuantity){

	Grid *grid = gridQuantity->grid;
	int *nGPoints = grid->nGPoints;
	long int *nGPointsProd = grid->nGPointsProd;
	int nDims = grid->nDims;

	msg(STATUS|ONCE, "Dumps 2D grid to parsefile");

	fMsg(ini,"parsedump", "\nDump of 3D grid: (%dx%dx%d) \n \n",
	 			nGPoints[0], nGPoints[1], nGPoints[2]);

	//Cycles trough and prints the grid (not optimized)
	int p;
	for(int l = 0; l < nGPoints[2]; l++){
		fMsg(ini, "parsedump", "\t\t\t\t\t l = %d \n", l);
		for(int k = nGPoints[1] - 1; k > -1; k--){ //y-rows
			for(int j = 0; j < nGPoints[0]; j++){ //x-rows
				p = j*nGPointsProd[0] + k*nGPointsProd[1] + l*nGPointsProd[2];
				fMsg(ini,"parsedump", "%5d", (int) gridQuantity->val[p]);
			}
			fMsg(ini,"parsedump", "\n\n");
		}
	}
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
