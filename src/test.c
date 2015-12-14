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
 * 		NOT TO BE INCLUDED IN FINAL PRODUCT
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include "test.h"

/******************************************************************************
 *	USEFUL CONSTANTS FOR TESTFUNCTIONS
 *****************************************************************************/
#define PI 3.14159265

/*****************************************************************************
 *					DEFINITIONS
 ****************************************************************************/

void testGridAndMGStructs(dictionary *ini, GridQuantity *gridQuantity, Multigrid *multigrid){
/*	Grid *grid = gridQuantity->grid;

	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}

	msg(STATUS|ONCE, "Performing a manual test of grid structs");

	//Writing information to parsedump (debugging)
	multigridParseDump(ini, multigrid);

	//Changing them from the multigrid struct
	for(int i = 0; i < 5; i++){
		multigrid->gridQuantities[0]->val[i] = 5.;
	}

	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES FROM MULTIGRID \n");

	//Check that both are changed
	//Writing information to parsedump (debugging)
	dumpGrid(ini, gridQuantity);
	multigridParseDump(ini, multigrid);


	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES BACK FROM GRIDQUANTITY \n");

	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}
	dumpGrid(ini,gridQuantity);
	multigridParseDump(ini, multigrid);
*/
	return;
}

void testGaussSeidel(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,
					MpiInfo *mpiInfo){
	msg(STATUS|ONCE, "Hello from GaussSeidel test");
	//Load Grid info
	GridQuantity *rho = multiRho->gridQuantities[0];
	GridQuantity *phi = multiPhi->gridQuantities[0];
	Grid *grid = rho->grid;
	int *nGPoints = grid->nGPoints;

	//Load GridQuantity
	double *rhoVal = rho->val;
	double *phiVal = phi->val;

	//Temp quick functions
	// double sin(double);

	//Variables
	// double angle;

	int ind = 0;
	for(int j = 0; j < nGPoints[0]; j++){
		// angle = 0;
		for (int k = 0; k<nGPoints[1]; k++) {
			// if(k>(0.5*nGPoints[1])) angle = 1.;
			rhoVal[ind] = (double) (j + 2*k);
			ind++;
		}
	}

	dumpGridIndexes(ini, rho);
	dumpGrid(ini,rho);

	for(int q = 0; q < 5; q++) multiRho->preSmooth(phi, rho);

	dumpGrid(ini, phi);

	for(int q = 0; q < 5; q++) multiRho->preSmooth(phi, rho);

	dumpGrid(ini, phi);


	return;
}

void testRestriction(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,
					MpiInfo *mpiInfo){

	msg(STATUS|ONCE, "Hello from Restriction test, not done!");
	//Load Grid info
	GridQuantity *rho = multiRho->gridQuantities[0];
	// GridQuantity *phi = multiPhi->gridQuantities[0];
	Grid *grid = rho->grid;
	int *nGPoints = grid->nGPoints;

	//Load GridQuantity
	double *rhoVal = rho->val;


	int ind = 0;
	for(int j = 0; j < nGPoints[0]; j++){
		for (int k = 0; k<nGPoints[1]; k++) {
			if(k < 3)rhoVal[ind] = (double) 1;
			ind++;
		}
	}

	//Load fine and coarse Grid
	GridQuantity *fRho = multiRho->gridQuantities[0];
	GridQuantity *cRho = multiRho->gridQuantities[1];

	// dumpGridIndexes(ini, fRho);
	// dumpGridIndexes(ini, cRho);
	dumpGrid(ini, fRho);

	multiRho->restrictor(multiRho->gridQuantities[0], multiRho->gridQuantities[1]);
	multiRho->prolongator(multiRho->gridQuantities[0], multiRho->gridQuantities[1]);

	dumpGrid(ini, fRho);


	// multiRho->prolongator(multiRho->gridQuantities[0], multiRho->gridQuantities[1]);

	return;
}

void testMultigrid(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,
					MpiInfo *mpiInfo){

	//Load
	GridQuantity *rho = multiRho->gridQuantities[0];
	Grid *grid = rho->grid;
	long int *nGPointsProd = grid->nGPointsProd;
	int nDims = grid->nDims;

	//Load GridQuantity
	double *rhoVal = rho->val;

	for(int ind = 0; ind < nGPointsProd[nDims]; ind++) rhoVal[ind] = 1.;

	dumpGrid(ini, rho);

	linearMGSolv(multiRho,multiPhi);


	return;
}

void testSwapHalo(dictionary *ini, GridQuantity *gridQuantity, MpiInfo *mpiInfo){

	//Load MPI info
	int mpiRank = mpiInfo->mpiRank;
	int mpiSize = mpiInfo->mpiSize;
	int *subdomain = mpiInfo->subdomain;


	//Load Grid info
	Grid *grid = gridQuantity->grid;
//	int nDims = grid->nDims;
	int *nGPoints = grid->nGPoints;

	//Load GridQuantity
	double *val = gridQuantity->val;

	int ind = 0;
	for(int j = 0; j < nGPoints[0]; j++){
		for (int k = 0; k<nGPoints[1]; k++) {
			val[ind] = (double) mpiRank;
			ind++;
		}
	}

	for(int rank = 0; rank < mpiSize; rank++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(mpiRank == rank){
			fMsg(ini, "parsedump", "rank = %d,\tsubdomain = [%d,%d]",
 					mpiRank, subdomain[0], subdomain[1]);
			dumpGrid(ini, gridQuantity);
		}
	}

	swapHalo(ini, gridQuantity, mpiInfo, 1);

	if(mpiRank == 0) fMsg(ini, "parsedump", "\n\nSwapping halos\n\n");

	for(int rank = 0; rank < mpiSize; rank++){
		MPI_Barrier(MPI_COMM_WORLD);
		if(mpiRank == rank){
			fMsg(ini, "parsedump", "rank = %d,\tsubdomain = [%d,%d]",
 					mpiRank, subdomain[0], subdomain[1]);
			dumpGrid(ini, gridQuantity);
		}
	}


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

	// swapHalo(ini, gridQuantity);

	dumpGrid(ini, gridQuantity);
}

void dumpGridIndexes(dictionary *ini, GridQuantity *gridQuantity){
	Grid *grid = gridQuantity->grid;
	int *nGPoints = grid->nGPoints;
	long int *nGPointsProd = grid->nGPointsProd;
	int nDims = grid->nDims;

	if(nDims == 3){
		fMsg(ini,"parsedump", "\nDump of 3D grid: (%dx%dx%d) \n \n",
					nGPoints[0], nGPoints[1], nGPoints[2]);
		//Cycles trough and prints the grid (not optimized)
		int p;
		for(int l = 0; l < nGPoints[2]; l++){
			fMsg(ini, "parsedump", "\t\t\t l = %d \n", l);
			for(int k = nGPoints[1] - 1; k > -1; k--){ //y-rows
				for(int j = 0; j < nGPoints[0]; j++){ //x-rows
					p = j*nGPointsProd[0] + k*nGPointsProd[1] + l*nGPointsProd[2];
					fMsg(ini,"parsedump", "%5d", p);
				}
				fMsg(ini,"parsedump", "\n\n");
			}
		}
	} else if(nDims==2) {
		fMsg(ini,"parsedump", "\t\t 2D grid: (%dx%d): \n \n",
					nGPoints[0], nGPoints[1]);
		int p;
		for(int k = nGPoints[1] - 1; k > -1; k--){ //y-rows
			for(int j = 0; j < nGPoints[0]; j++){ //x-rows
				p = j*nGPointsProd[0] + k*nGPointsProd[1];
				fMsg(ini,"parsedump", "%5d", p);
			}
			fMsg(ini,"parsedump", "\n\n");
		}

	}
}

void dumpGrid(dictionary *ini, GridQuantity *gridQuantity){

	Grid *grid = gridQuantity->grid;
	int *nGPoints = grid->nGPoints;
	long int *nGPointsProd = grid->nGPointsProd;
	int nDims = grid->nDims;

	msg(STATUS|ONCE, "Dumps grid to parsefile");
	if(nDims == 3){
		fMsg(ini,"parsedump", "\nDump of 3D grid: (%dx%dx%d) \n \n",
		 			nGPoints[0], nGPoints[1], nGPoints[2]);
		//Cycles trough and prints the grid (not optimized)
		int p;
		for(int l = 0; l < nGPoints[2]; l++){
			fMsg(ini, "parsedump", "\t\t\t l = %d \n", l);
			for(int k = nGPoints[1] - 1; k > -1; k--){ //y-rows
				for(int j = 0; j < nGPoints[0]; j++){ //x-rows
					p = j*nGPointsProd[0] + k*nGPointsProd[1] + l*nGPointsProd[2];
					fMsg(ini,"parsedump", "%5d", (int) gridQuantity->val[p]);
				}
				fMsg(ini,"parsedump", "\n\n");
			}
		}
	} else if(nDims==2) {
		fMsg(ini,"parsedump", "\t\t 2D grid: (%dx%d): \n \n",
		 			nGPoints[0], nGPoints[1]);
		int p;
		for(int k = nGPoints[1] - 1; k > -1; k--){ //y-rows
			for(int j = 0; j < nGPoints[0]; j++){ //x-rows
				p = j*nGPointsProd[0] + k*nGPointsProd[1];
				fMsg(ini,"parsedump", "%5d", (int) gridQuantity->val[p]);
			}
			fMsg(ini,"parsedump", "\n\n");
		}

	}

	return;
}
