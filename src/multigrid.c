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

	char *preSmoothName = iniparser_getstring((dictionary*)ini, "modules:preSmooth", "\0");
    char *postSmoothName = iniparser_getstring((dictionary*)ini, "modules:postSmooth", "\0");
    char *coarseSolverName = iniparser_getstring((dictionary*)ini, "modules:coarseSolv", "\0");

	int nDims = multigrid->grids[0]->rank-1;

	if(!strcmp(preSmoothName,"gaussSeidel")){
		if(nDims == 2)	multigrid->preSmooth = &gaussSeidel2D;
		// else if(nDims == 3) multigrid->preSmooth = &gaussSeidel3D;

    } else {
    	msg(ERROR, "No Presmoothing algorithm specified");
    }

    if(!strcmp(postSmoothName,"gaussSeidel")){
		multigrid->postSmooth = &gaussSeidel2D;
 	} else {
    	msg(ERROR, "No Postsmoothing algorithm specified");
    }

    if(!strcmp(coarseSolverName,"gaussSeidel")){
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

	int rank = multigrid->grids[0]->rank;

	if(!strcmp(restrictor, "halfWeight")){
		if(rank == 3)	multigrid->restrictor = &halfWeightRestrict2D;
		else if(rank == 4) multigrid->restrictor = &halfWeightRestrict3D;
		else msg(ERROR, "No restricting algorithm for D%d", rank-1);
	} else {
		msg(ERROR, "No restrict stencil specified");
	}
	if(!strcmp(prolongator, "bilinear")){
		if(rank == 3)	multigrid->prolongator = &bilinearProlong2D;
		else if(rank==4)	multigrid->prolongator = &bilinearProlong3D;
		else msg(ERROR, "No restricting algorithm for D%d", rank-1);
	} else {
		msg(ERROR, "No prolongation stencil specified");
	}

}


Grid **mgAllocSubGrids(const dictionary *ini, Grid *grid,
						const int nLevels){

	//Gather information on finest grid
	Grid **grids = malloc((nLevels+1) * sizeof(Grid));

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


void jacobian(Grid *phi,const Grid *rho, const int nCycles){

	//Common variables
	int rank = phi->rank;
	long int *sizeProd = phi->sizeProd;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;


	// //Debug stuff
	// int ind = 3*sizeProd[1] + 3*sizeProd[2];
	//
	// msg(STATUS, "grid point #%d at point: [3,3] = %f", ind, rhoVal[ind]);


	double *tempVal = malloc (sizeProd[rank]*sizeof(*tempVal));
	for(int c = 0; c < nCycles; c++){
		// Index of neighboring nodes
		int gj = sizeProd[1];
		int gjj= -sizeProd[1];
		int gk = sizeProd[2];
		int gkk= -sizeProd[2];



		for(int g = 0; g < sizeProd[rank]; g++){
			tempVal[g] = 0.25*(	phiVal[gj] + phiVal[gjj] +
								phiVal[gk] + phiVal[gkk] + rhoVal[g]);

			// //Debug Stuff
			// if(g == ind){
			// 	msg(STATUS)
			// 	msg(STATUS, "i+: %f", rhoVal[gj]);
			// 	msg(STATUS, "i-: %f", rhoVal[gjj]);
			// 	msg(STATUS, "j+: %f", rhoVal[gk]);
			// 	msg(STATUS, "j-: %f", rhoVal[gkk]);
			// }

			gj++;
			gjj++;
			gk++;
			gkk++;
		}

		for(int q = 0; q < sizeProd[rank]; q++) phiVal[q] = tempVal[q];
	}

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

void halfWeightRestrict3D(const Grid *fine, Grid *coarse){

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;
	int *cSize = coarse->size;

	//Indexes
	long int c = cSizeProd[1] + cSizeProd[2] + cSizeProd[3];

	long int f = fSizeProd[1] + fSizeProd[2] + fSizeProd[3];
	long int fj  = f + fSizeProd[1];
	long int fjj = f - fSizeProd[1];
	long int fk  = f + fSizeProd[2];
	long int fkk = f - fSizeProd[2];
	long int fl  = f + fSizeProd[3];
 	long int fll = f - fSizeProd[3];

	int cKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2];
	int fKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2] + fSizeProd[2];
	int cLEdgeInc = (nGhostLayers[3] + nGhostLayers[rank + 3])*cSizeProd[2];
	int fLEdgeInc = (nGhostLayers[3] + nGhostLayers[rank + 3])*fSizeProd[2] + fSizeProd[3];

	int coeff = 1./12.;
	//
	// msg(STATUS, "Coarse SizeProds = [%d, %d, %d, %d]", cSizeProd[0], cSizeProd[1], cSizeProd[2], cSizeProd[3]);
	// msg(STATUS, "Fine SizeProds = [%d, %d, %d, %d]", fSizeProd[0], fSizeProd[1], fSizeProd[2], fSizeProd[3]);
	// // msg(STATUS, "nSize: %d", cSizeProd[rank]);

	//Cycle Coarse grid
	for(int l = nGhostLayers[3]; l<cSize[3]-nGhostLayers[rank+3]; l++){
		for(int k = nGhostLayers[2]; k < cSize[2]-nGhostLayers[rank + 2]; k++){
			for(int j = nGhostLayers[1]; j < cSize[1]-nGhostLayers[rank+1]; j++){
				// msg(STATUS, "c=%d, f = [%d] (%d %d %d %d %d %d)", c, f , fj, fjj, fk, fkk, fl, fll);
				cVal[c] = coeff*(6*fVal[f] + fVal[fj] + fVal[fjj] + fVal[fk] + fVal[fkk] + fVal[fl] + fVal[fll]);
				c++;
				f  +=2;
				fj +=2;
				fjj+=2;
				fk +=2;
				fkk+=2;
				fl +=2;
				fll+=2;
			}
			c  += cKEdgeInc;
			f  += fKEdgeInc;
			fj += fKEdgeInc;
			fjj+= fKEdgeInc;
			fk += fKEdgeInc;
			fkk+= fKEdgeInc;
			fl += fKEdgeInc;
			fll+= fKEdgeInc;
		}
		c  += cLEdgeInc;
		f  += fLEdgeInc;
		fj += fLEdgeInc;
		fjj+= fLEdgeInc;
		fk += fLEdgeInc;
		fkk+= fLEdgeInc;
		fl += fLEdgeInc;
		fll+= fLEdgeInc;
	}

	return;
}

void halfWeightRestrict2D(const Grid *fine, Grid *coarse){
	msg(STATUS, "Hello from restrictor");

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;
	int *cSize = coarse->size;

	//Indexes
	long int c = cSizeProd[2] + cSizeProd[1];

	long int f = fSizeProd[2] + fSizeProd[1];
	long int fj = f + fSizeProd[1];
	long int fjj = f - fSizeProd[1];
	long int fk = f + fSizeProd[2];
	long int fkk = f - fSizeProd[2];

	int cKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2];
	int fKEdgeInc = cKEdgeInc + fSizeProd[2];

	//Cycle Coarse grid
	for(int k = nGhostLayers[2]; k < cSize[2]-nGhostLayers[rank + 2]; k++){
		for(int j = nGhostLayers[1]; j < cSize[1]-nGhostLayers[rank+1]; j++){
			// msg(STATUS, "c=%d, f = [%d] (%d %d %d %d)", c, f , fj, fjj, fk, fkk);
			cVal[c] = 0.125*(4*fVal[f] + fVal[fj] + fVal[fjj] + fVal[fk] + fVal[fkk]);
			c++;
			f  +=2;
			fj +=2;
			fjj+=2;
			fk +=2;
			fkk+=2;
		}
		c  += cKEdgeInc;
		f  += fKEdgeInc;
		fj += fKEdgeInc;
		fjj+= fKEdgeInc;
		fk += fKEdgeInc;
		fkk+= fKEdgeInc;
	}

	return;
}

void bilinearProlong3D(Grid *fine, const Grid *coarse){
	return;
}


void bilinearProlong2D(Grid *fine, const Grid *coarse){

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int *fSize = fine->sizeProd;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;

	//Indexes
	long int f = fSizeProd[2] + fSizeProd[1];

	long int c = cSizeProd[2] + cSizeProd[1];

	int cKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2];
	int fKEdgeInc = cKEdgeInc + fSizeProd[2];

	//Direct insertion c->f
	for(int k = nGhostLayers[2]; k < fSize[2]-nGhostLayers[rank + 2]; k++){
		for(int j = nGhostLayers[1]; j < fSize[1]-nGhostLayers[rank+1]; j++){

		}
	}
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
