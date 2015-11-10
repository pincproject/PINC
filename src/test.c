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
	int nDims = gridQuantity->grid->nDims;
	int *nGPoints = gridQuantity->grid->nGPoints;
	int *nGPointsProd = gridQuantity->grid->nGPointsProd;
	int totalGPoints = nGPointsProd[nDims];
	
	msg(STATUS|ONCE, "Total grid points: \t %d, nGPointsProd = [%d , %d]",\
	 totalGPoints, nGPointsProd[0], nGPointsProd[1]);

	//Populate grid as 0.
	for(int g = 0; g < totalGPoints; g++){
		gridQuantity->val[g]=0.;
	}

	int l = 0;	//Lower edge
	int h = 0;	//Higher edge
	int b = 1;	//Edge counter
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
	}

	
	dump2DGrid(ini, gridQuantity);
	dumpGhostVector(ini, gridQuantity);

	return;
}

void dump2DGrid(dictionary *ini, GridQuantity *gridQuantity){
	msg(STATUS|ONCE, "Dumps 2D grid to parsefile");

	int nDims = gridQuantity->grid->nDims;
	int *nGPoints = gridQuantity->grid->nGPoints;
	int *nGPointsProd = gridQuantity->grid->nGPointsProd;

	fMsg(ini,"parsedump", "Dump of 2D/1D indexes: (%dx%d) \n \n", nGPoints[0], nGPoints[1]);

	int p = 0;


	for(int k = 0; k < nGPoints[1]; k++){
		for(int j = 0; j < nGPoints[0]; j++){	
			fMsg(ini,"parsedump", "\t%d", p);
			p++;
		}
		fMsg(ini,"parsedump", "\n ");	
	}

	fMsg(ini,"parsedump", "\n Dump of 2D/1D grid: (%dx%d) \n \n", nGPoints[0], nGPoints[1]);

	p = 0;
	for(int k = 0; k < nGPoints[1]; k++){ //y-rows
		for(int j = 0; j < nGPoints[0]; j++){ //x-rows
			fMsg(ini,"parsedump", "\t%d", (int) gridQuantity->val[p]);	
			p++;
		}
		fMsg(ini,"parsedump", "\n ");	
	}
	
	return;
}

void dumpGhostVector(dictionary *ini, GridQuantity *gridQuantity){

	//Load
	Grid *grid = gridQuantity->grid;
	int *nGhosts = grid->nGhosts;
	int *nGPoints = grid->nGPoints;
	int nDims = grid->nDims;


	double *ghostEdge = getGhostEdge(ini, gridQuantity);

	int nGhostPoints = 0;
	for(int g = 0; g < nDims; g++){
		nGhostPoints += (nGhosts[g]+nGhosts[g +nDims])*nGPoints[g];
	}

	//Print GhostEdge
	for(int w = 0; w < nGhostPoints; w++){
		msg(STATUS, "ghostEdge[%d] = %d", w, (int) ghostEdge[w]);
	}
	
	return;
}