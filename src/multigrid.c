/**
 * @file		multigrid.c
 * @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
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

	if(!strcmp(preSmoothName,"gaussSeidel2D")){
		multigrid->preSmooth = &gaussSeidel2D;

    } else {
    	msg(ERROR, "No Presmoothing algorithm specified");
    }

    if(!strcmp(postSmoothName,"gaussSeidel2D")){
		multigrid->postSmooth = &gaussSeidel2D;
 	} else {
    	msg(ERROR, "No Postsmoothing algorithm specified");
    }

    if(!strcmp(coarseSolverName,"gaussSeidel2D")){
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


Grid **mgAllocSubGrids(const dictionary *ini, Grid *grid,
						const int nLevels){

	//Gather information on finest grid
	Grid **grids = malloc(nLevels * sizeof(*grids));

	int *trueSize = grid->trueSize;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;

	//Set first grid to point to f grid
	grids[0] = grid;

	//Cycle through subgrids
	for(int q = 1; q < nLevels; q++){

		//Allocate
		int *subTrueSize = malloc(rank*sizeof(subTrueSize));
		int *subSize = malloc(rank *sizeof(*subSize));

		//Set first entries
		subTrueSize[0] = trueSize[0];
		subSize[0] = subTrueSize[0];

		//The subgrid needs half the grid points
		for(int d = 1; d < rank; d++){
			if (!(trueSize[d]%(2*q)))	subTrueSize[d] = trueSize[d]/(2*q);
			else msg(ERROR, "Size of true grid in dim %d, is not divisible by 2*%d", d-1, q);
		}

		// Calculate the number of grid points (True points + ghost points)
		subSize[0] = trueSize[0];
		for(int d = 1 ; d < rank ; d ++){
			subSize[d] = subTrueSize[d];
			subSize[d] += nGhostLayers[d];
			subSize[d] += nGhostLayers[rank + d];
		}

		//Slice elements
		long int nSliceMax = 0;
		for(int d=0;d<rank;d++){
			long int nSlice = 1;
			for(int dd=0;dd<rank;dd++){
				if(dd!=d) nSlice *= subSize[dd];
			}
			if(nSlice>nSliceMax) nSliceMax = nSlice;
		}

		long int *subSizeProd = longIntArrCumProd(subSize,rank);

		double *val = malloc(subSizeProd[rank]*sizeof(*val));
		double *slice = malloc(nSliceMax*sizeof(*slice));

		//Assign to grid
		Grid *grid = malloc(sizeof(Grid));
		grid->rank = rank;
		grid->size = subSize;
		grid->sizeProd = subSizeProd;
		grid->nGhostLayers = nGhostLayers;
		grid->trueSize = subTrueSize;
		grid->val = val;
		grid->slice = slice;
		grid->h5 = 0;
		grids[q] = grid;
	}

	return grids;
}


/*************************************************
 *		DEFINITIONS
 ************************************************/

Multigrid *mgAlloc(const dictionary *ini, Grid *grid){

	//Multigrid
	int nLevels = iniparser_getint((dictionary *) ini, "multigrid:mgLevels", 0);
	int nMGCycles = iniparser_getint((dictionary *) ini, "multigrid:mgCycles", 0);
	int nPreSmooth = iniparser_getint((dictionary *) ini, "multigrid:nPreSmooth", 0);
	int nPostSmooth = iniparser_getint((dictionary *) ini, "multigrid:nPostSmooth", 0);
	int nCoarseSolve = iniparser_getint((dictionary *) ini, "multigrid:nCoarseSolve", 0);
	//Load data
	int nDims = grid->rank-1;
	int *trueSize = grid->trueSize;


	//Sanity checks
	if(!nLevels) msg(ERROR, "Multi Grid levels is 0, direct solver not implemented yet \n");


	if(!nMGCycles) msg(ERROR, "MG cycles is 0 \n");


	// Sanity check (true grid points need to be a multiple of 2^(multigrid levels)
	for(int d = 0; d < nDims; d++){
		if(trueSize[d+1] % (int) pow(2, nLevels)){ //Sloppy and wrong
			msg(ERROR, "The number of True Grid Points needs to be a multiple of 2^nLevels");
		}
	}

	Grid **grids = mgAllocSubGrids(ini, grid, nLevels);

	//Store in multigrid struct
    Multigrid *multigrid = malloc(sizeof(Multigrid));

    multigrid->nLevels = nLevels;
    multigrid->nMGCycles = nMGCycles;
	multigrid->nPreSmooth = nPreSmooth;
	multigrid->nPostSmooth = nPostSmooth;
	multigrid->nCoarseSolve = nCoarseSolve;
    multigrid->grids = grids;

    //Setting the algorithms to be used, pointer functions
	setSolvers(ini, multigrid);
	setRestrictProlong(ini, multigrid);

  	return multigrid;

}

void mgFree(Multigrid *multigrid){
	//
	// Grid **grids = multigrid->grids;
	// int nLevels = multigrid->nLevels;

	// for(int n = 1; n < nLevels; n++)
	// {
	// 	gFree(grids[n]);
	// }
	// free(multigrid);

	msg(STATUS, "Freeing of multigrid not complete");
	return;
}


void jacobian(Grid *phi,const Grid *rho, int nCycles){
	printf("Hello from Jacobian \n");
	return;
}

void gaussSeidel2D(Grid *phi, const Grid *rho, int nCycles){

	//Common variables
	int *trueSize = phi->trueSize;
	long int *sizeProd = phi->sizeProd;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	for(int c = 0; c < nCycles;c++){
		long int ind = sizeProd[1] + 1;
		int j;
		//Red Pass
		// msg(STATUS, "Red indexes = ");
		for(int k = 1; k < trueSize[1] + 1; k ++){
			if(k%2) j = 1; else j = 2;
			for(; j < trueSize[0] + 1; j += 2){
				ind = j*sizeProd[0] + k*sizeProd[1];
				phiVal[ind] = (phiVal[ind+sizeProd[0]] + phiVal[ind-sizeProd[0]] +
							phiVal[ind+sizeProd[1]] + phiVal[ind-sizeProd[1]] -
							rhoVal[ind])*0.25;
				// msg(STATUS, "%d", ind);
			}
		}

		//Black Pass
		// msg(STATUS, "Black indexes = ");
		for(int k = 1; k < trueSize[1] + 1; k ++){
			if(k%2) j = 2; else j = 1;
			for(; j < trueSize[0] + 1; j += 2){
				ind = j*sizeProd[0] + k*sizeProd[1];
				phiVal[ind] = (phiVal[ind+sizeProd[0]] + phiVal[ind-sizeProd[0]] +
							phiVal[ind+sizeProd[1]] + phiVal[ind-sizeProd[1]] -
							rhoVal[ind])*0.25;
				// msg(STATUS, "%d", ind);
			}
		}
	}


	return;
}

void halfWeightRestrict(const Grid *fine, Grid *coarse){
	// msg(STATUS, "Hello from restrictor");
	//
	// //Load f grid
	// double *fVal = fine->val;
	// long int *fProd = fine->sizeProd;
	//
	// //Load c grid
	// double *cVal = coarse->val;
	// int *ctrueSize = coarse->trueSize;
	// long int *cProd = coarse->sizeProd;
	//
	// //Coarse and fine indexes
	// long int cInd;
	// long int fInd;
	//
	// int kF = 1;
	// int jF = 1;
	//
	// //Cycle through the c grid
	// //Spagetti code (To be improved)
	// for(int kC=1; kC < ctrueSize[1]+1; kC++){
	// 	for(int jC=1; jC < ctrueSize[0]+1;jC++){
	// 		cInd = jC*cProd[0] + kC*cProd[1];
	// 		fInd = jF*fProd[0] + kF*fProd[1];
	// 		cVal[cInd] = 0.125*(4*fVal[fInd] + fVal[fInd+fProd[0]] + fVal[fInd-fProd[0]]
 // 						+ fVal[fInd+fProd[1]] + fVal[fInd-fProd[1]]);
	// 		jF += 2;
	// 	}
	// 	kF += 2;
	// 	jF = 1;
	// }
	//
	//

	return;
}

void bilinearProlong(Grid *fine, const Grid *coarse){

	// //Load f grid
	// double *fVal = f->val;
	// Grid *fGrid = f->grid;
	// long int *fProd = fGrid->sizeProd;
	// int *ftrueSize = fGrid->trueSize;
	//
	//
	// //Load c grid
	// double *cVal = c->val;
	// Grid *cGrid = c->grid;
	// long int *cProd = cGrid->sizeProd;
	// int *ctrueSize = cGrid->trueSize;
	//
	// //Coarse and fine indexes
	// long int cInd;
	// long int fInd;
	//
	// int kF = 1;
	// int jF = 1;
	//
	// //Cycle through the c grid, and copy in onto corresponding fGrid points
	// //Spagetti code (To be improved)
	// for(int kC=1; kC < ctrueSize[1]+1; kC++){
	// 	for(int jC=1; jC < ctrueSize[0]+1;jC++){
	// 		cInd = jC*cProd[0] + kC*cProd[1];
	// 		fInd = jF*fProd[0] + kF*fProd[1];
	// 		fVal[fInd] = cVal[cInd];
	// 		jF += 2;
	// 	}
	// 	kF += 2;
	// 	jF = 1;
	// }
	//
	// //Odd numbered columns, interpolating vertically
	// for(int kF=1; kF<ftrueSize[1]+1;kF+=2){
	// 	for(int jF=2; jF<ftrueSize[0]+1;jF+=2){
	// 			fInd = jF*fProd[0]+kF*fProd[1];
	// 			fVal[fInd] = 0.5*(fVal[fInd+fProd[0]] + fVal[fInd-fProd[0]]);
	// 	}
	// }
	//
	// //Even numbered columns, interpolating horizontally
	// for(int kF=2; kF<ftrueSize[1]+1;kF+=2){
	// 	for(int jF=1; jF<ftrueSize[0]+1;jF+=2){
	// 			fInd = jF*fProd[0]+kF*fProd[1];
	// 			fVal[fInd] = 0.5*(fVal[fInd+fProd[1]] + fVal[fInd-fProd[1]]);
	// 	}
	// }



	return;
}



void linearMGSolv(Multigrid *multiRho, Multigrid *multiPhi){

	int nMGCycles = multiRho->nMGCycles;
	int nLevels = multiRho->nLevels;
	int nCoarseSolve = multiRho->nCoarseSolve;
	int nPreSmooth = multiRho->nPreSmooth;

	Grid *rho = multiRho->grids[0];
	Grid *phi = multiPhi->grids[0];

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
