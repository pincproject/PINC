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
	int *nGPointsProd = grid->nGPointsProd;
	// int *node = grid->node;

	msg(STATUS|ONCE, "Total grid points: \t %d, nGPointsProd = [%d , %d]",\
	 nGPointsProd[nDims], nGPointsProd[0], nGPointsProd[1]);


	//Get rank to put as grid values
	//rank 1 grid has all values 1 and so on
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	//Populate grid as 0.
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
			dump2DGrid(ini, gridQuantity);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	swapHalo(ini, gridQuantity);
	for(int r = 0; r < size; r++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank == r){
			dump2DGrid(ini, gridQuantity);
		}
	}

	return;
}


void dump2DGrid(dictionary *ini, GridQuantity *gridQuantity){

	int *node = gridQuantity->grid->node;

	msg(STATUS|ONCE, "Dumps 2D grid to parsefile, node[%d,%d]",node[0], node[1]);

	// int nDims = gridQuantity->grid->nDims;
	int *nGPoints = gridQuantity->grid->nGPoints;
	// int *nGPointsProd = gridQuantity->grid->nGPointsProd;


	int p = 0;

	if((node[0] == 0) & (node[1]==0)){

		fMsg(ini,"parsedump", "indexes: (%dx%d) \n \n", nGPoints[0], nGPoints[1]);

		for(int k = 0; k < nGPoints[1]; k++){
			for(int j = 0; j < nGPoints[0]; j++){
				fMsg(ini,"parsedump", "\t%d ", p);
				p++;
			}
		fMsg(ini,"parsedump", "\n ");
		}
	}

	fMsg(ini,"parsedump", "\nDump of 2D/1D grid: (%dx%d) and node[%d,%d] \n \n", nGPoints[0], nGPoints[1],  node[0], node[1]);

	p = 0;
	for(int k = 0; k < nGPoints[1]; k++){ //y-rows
		for(int j = 0; j < nGPoints[0]; j++){ //x-rows
			fMsg(ini,"parsedump", "\t%d", (int) gridQuantity->val[p]);
			p++;
		}
		fMsg(ini,"parsedump", "\n");
	}

	return;
}

void dumpHalo(dictionary *ini, GridQuantity *gridQuantity){

	//Load
	Grid *grid = gridQuantity->grid;
	int *nGhosts = grid->nGhosts;
	int *nGPoints = grid->nGPoints;
	int nDims = grid->nDims;
	double *halo = gridQuantity->halo;

	int nGhostPoints = 0;
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g +nDims])*nGPoints[g];
	}

	fMsg(ini , "parsedump", "\n******************************************************\n\n");
	fMsg(ini , "parsedump", "node[%d,%d]: ghostEdge = \n", grid->node[0],grid->node[1]);
	//Print GhostEdge
	for(int w = 0; w < nGhostPoints; w++){
		fMsg(ini, "parsedump" , "%d,",(int) halo[w]);
	}
	fMsg(ini , "parsedump", "\n\n******************************************************\n");


	return;
}
