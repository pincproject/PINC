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

	int nDims = multigrid->gridQuantities[0]->grid->nDims;

	if(!strcmp(preSmoothName,"gaussSeidel") && nDims == 2){
		multigrid->preSmooth = &gaussSeidel2D;

    } else {
    	msg(ERROR, "No Presmoothing algorithm specified");
    }

    if(!strcmp(postSmoothName,"gaussSeidel") && nDims == 2){
		multigrid->postSmooth = &gaussSeidel2D;
 	} else {
    	msg(ERROR, "No Postsmoothing algorithm specified");
    }

    if(!strcmp(coarseSolverName,"gaussSeidel") && nDims == 2){
		multigrid->coarseSolv = &gaussSeidel2D;
    } else {
    	msg(ERROR, "No coarse Grid Solver algorithm specified");
    }

    //Free!

//    printf("If this is the last I say, I'm bad and segfaults at freeing stuff \n");
    // free(cSolverName);
    // free(postSmoothName);
    // free(preSmoothName);

	return;
}

void setRestrictProlong(const dictionary *ini,Multigrid *multigrid){
	char *restrictor = iniparser_getstring((dictionary*)ini, "multigrid:restrictor", "\0");
	char *prolongator = iniparser_getstring((dictionary*)ini, "multigrid:prolongator", "\0");

	if(!strcmp(restrictor, "halfWeight")){
		multigrid->restrictor = &halfWeightRestrict;
	} else {
		msg(ERROR, "No restrict stencil specified");
	}
	if(!strcmp(prolongator, "bilinear")){
		multigrid->prolongator = &bilinearProlong;
	} else {
		msg(ERROR, "No prolongation stencil specified");
	}

}


GridQuantity **allocSubGrids(const dictionary *ini, GridQuantity *gridQuantity,
							const int nLevels){

	GridQuantity **gridQuantities = malloc(nLevels * sizeof(GridQuantity));
	int nDims;
	int nBoundaries;

	//Set first grid to point to f grid
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
		long int *nGPointsProd = longIntArrCumProd(nGPoints,nDims);

		//Assign to grid
		Grid *grid = malloc(sizeof(Grid));

		grid->nDims = nDims;
		grid->nGPoints = nGPoints;
		grid->nGPointsProd = nGPointsProd;
		grid->nGhosts = nGhosts;
		grid->nTGPoints = nTGPoints;
		//Creates GridQuantity
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
	int nMGCycles = iniparser_getint((dictionary *) ini, "multigrid:mgCycles", 0);
	int nPreSmooth = iniparser_getint((dictionary *) ini, "multigrid:nPreSmooth", 0);
	int nPostSmooth = iniparser_getint((dictionary *) ini, "multigrid:nPostSmooth", 0);
	int nCoarseSolve = iniparser_getint((dictionary *) ini, "multigrid:nCoarseSolve", 0);
	//Load data
	int nDims = gridQuantity->grid->nDims;
	int *nTGPoints = gridQuantity->grid->nTGPoints;

	//Sanity checks
	if(!nLevels) msg(ERROR, "Multi Grid levels is 0, direct solver not implemented yet \n");


	if(!nMGCycles) msg(ERROR, "MG cycles is 0 \n");


	// Sanity check (true grid points need to be a multiple of 2^(multigrid levels)
	for(int d = 0; d < nDims; d++){
		if(nTGPoints[d] % (int) pow(2, nLevels)){ //Sloppy and wrong
			msg(ERROR, "The number of True Grid Points needs to be a multiple of 2^nLevels");
		}
	}

	GridQuantity **gridQuantities = allocSubGrids(ini, gridQuantity, nLevels);

	//Store in multigrid struct
    Multigrid *multigrid = malloc(sizeof(Multigrid));

    multigrid->nLevels = nLevels;
    multigrid->nMGCycles = nMGCycles;
	multigrid->nPreSmooth = nPreSmooth;
	multigrid->nPostSmooth = nPostSmooth;
	multigrid->nCoarseSolve = nCoarseSolve;
    multigrid->gridQuantities = gridQuantities;

    //Setting the algorithms to be used, pointer functions
	setSolvers(ini, multigrid);
	setRestrictProlong(ini, multigrid);

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


void jacobian(GridQuantity *phi,const GridQuantity *rho, int nCycles){
	printf("Hello from Jacobian \n");
	return;
}

void gaussSeidel2D(GridQuantity *phi, const GridQuantity *rho, int nCycles){

	//Common variables
	Grid *grid = phi->grid;
	int *nTGPoints = grid->nTGPoints;
	long int *nGPointsProd = grid->nGPointsProd;



	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	for(int c = 0; c < nCycles;c++){
		long int ind = nGPointsProd[1] + 1;
		int j;
		//Red Pass
		// msg(STATUS, "Red indexes = ");
		for(int k = 1; k < nTGPoints[1] + 1; k ++){
			if(k%2) j = 1; else j = 2;
			for(; j < nTGPoints[0] + 1; j += 2){
				ind = j*nGPointsProd[0] + k*nGPointsProd[1];
				phiVal[ind] = (phiVal[ind+nGPointsProd[0]] + phiVal[ind-nGPointsProd[0]] +
							phiVal[ind+nGPointsProd[1]] + phiVal[ind-nGPointsProd[1]] -
							rhoVal[ind])*0.25;
				// msg(STATUS, "%d", ind);
			}
		}

		//Black Pass
		// msg(STATUS, "Black indexes = ");
		for(int k = 1; k < nTGPoints[1] + 1; k ++){
			if(k%2) j = 2; else j = 1;
			for(; j < nTGPoints[0] + 1; j += 2){
				ind = j*nGPointsProd[0] + k*nGPointsProd[1];
				phiVal[ind] = (phiVal[ind+nGPointsProd[0]] + phiVal[ind-nGPointsProd[0]] +
							phiVal[ind+nGPointsProd[1]] + phiVal[ind-nGPointsProd[1]] -
							rhoVal[ind])*0.25;
				// msg(STATUS, "%d", ind);
			}
		}
	}


	return;
}

void halfWeightRestrict(GridQuantity *f, GridQuantity *c){
	msg(STATUS, "Hello from restrictor");

	//Load f grid
	double *fVal = f->val;
	Grid *fGrid = f->grid;
	long int *fProd = fGrid->nGPointsProd;

	//Load c grid
	double *cVal = c->val;
	Grid *cGrid = c->grid;
	int *cNTGPoints = cGrid->nTGPoints;
	long int *cProd = cGrid->nGPointsProd;

	//Coarse and fine indexes
	long int cInd;
	long int fInd;

	int kF = 1;
	int jF = 1;

	//Cycle through the c grid
	//Spagetti code (To be improved)
	for(int kC=1; kC < cNTGPoints[1]+1; kC++){
		for(int jC=1; jC < cNTGPoints[0]+1;jC++){
			cInd = jC*cProd[0] + kC*cProd[1];
			fInd = jF*fProd[0] + kF*fProd[1];
			cVal[cInd] = 0.125*(4*fVal[fInd] + fVal[fInd+fProd[0]] + fVal[fInd-fProd[0]]
 						+ fVal[fInd+fProd[1]] + fVal[fInd-fProd[1]]);
			jF += 2;
		}
		kF += 2;
		jF = 1;
	}



	return;
}

void bilinearProlong(GridQuantity *f, GridQuantity *c){

	//Load f grid
	double *fVal = f->val;
	Grid *fGrid = f->grid;
	long int *fProd = fGrid->nGPointsProd;
	int *fNTGPoints = fGrid->nTGPoints;


	//Load c grid
	double *cVal = c->val;
	Grid *cGrid = c->grid;
	long int *cProd = cGrid->nGPointsProd;
	int *cNTGPoints = cGrid->nTGPoints;

	//Coarse and fine indexes
	long int cInd;
	long int fInd;

	int kF = 1;
	int jF = 1;

	//Cycle through the c grid, and copy in onto corresponding fGrid points
	//Spagetti code (To be improved)
	for(int kC=1; kC < cNTGPoints[1]+1; kC++){
		for(int jC=1; jC < cNTGPoints[0]+1;jC++){
			cInd = jC*cProd[0] + kC*cProd[1];
			fInd = jF*fProd[0] + kF*fProd[1];
			fVal[fInd] = cVal[cInd];
			jF += 2;
		}
		kF += 2;
		jF = 1;
	}

	//Odd numbered columns, interpolating vertically
	for(int kF=1; kF<fNTGPoints[1]+1;kF+=2){
		for(int jF=2; jF<fNTGPoints[0]+1;jF+=2){
				fInd = jF*fProd[0]+kF*fProd[1];
				fVal[fInd] = 0.5*(fVal[fInd+fProd[0]] + fVal[fInd-fProd[0]]);
		}
	}

	//Even numbered columns, interpolating horizontally
	for(int kF=2; kF<fNTGPoints[1]+1;kF+=2){
		for(int jF=1; jF<fNTGPoints[0]+1;jF+=2){
				fInd = jF*fProd[0]+kF*fProd[1];
				fVal[fInd] = 0.5*(fVal[fInd+fProd[1]] + fVal[fInd-fProd[1]]);
		}
	}



	return;
}



void linearMGSolv(Multigrid *multiRho, Multigrid *multiPhi){

	int nMGCycles = multiRho->nMGCycles;
	int nLevels = multiRho->nLevels;
	int nCoarseSolve = multiRho->nCoarseSolve;
	int nPreSmooth = multiRho->nPreSmooth;

	GridQuantity *rho = multiRho->gridQuantities[0];
	GridQuantity *phi = multiPhi->gridQuantities[0];

	// msg(STATUS, "%d", nLevels);

	for(int c = 0; c < nMGCycles; c++)
		for(int l = 0; l < nLevels; l++){
			if(l==nLevels-1) multiRho->coarseSolv(rho, phi, nCoarseSolve);	//Coarse grid
			else { //MG Rountine
				multiRho->preSmooth(rho,phi, nPreSmooth);
				//Defect calculation here
				//Restrict defect
				//Set inital guess
				//Call itself
			}
		}

	return;
}
