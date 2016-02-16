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
		else if(nDims == 3) multigrid->preSmooth = &gaussSeidel3D;
		else msg(ERROR, "No Presmoothing algorithm set for dimensions %d", nDims);
    } else if (!strcmp(preSmoothName, "jacobian2D")){
		multigrid->preSmooth = &jacobian2D;
 	} else {
    	msg(ERROR, "No Presmoothing algorithm specified");
    }

    if(!strcmp(postSmoothName,"gaussSeidel")){
		if(nDims == 2)	multigrid->postSmooth = &gaussSeidel2D;
		else if(nDims == 3) multigrid->postSmooth = &gaussSeidel3D;
		else msg(ERROR, "No postsmoothing algorithm set for dimensions %d", nDims);
	} else if (!strcmp(postSmoothName, "jacobian2D")){
		multigrid->postSmooth = &jacobian2D;
 	} else {
    	msg(ERROR, "No Postsmoothing algorithm specified");
    }

    if(!strcmp(coarseSolverName,"gaussSeidel")){
		if(nDims == 2)	multigrid->coarseSolv = &gaussSeidel2D;
		else if(nDims == 3) multigrid->coarseSolv = &gaussSeidel3D;
		else msg(ERROR, "No coarsesolver algorithm set for dimensions %d", nDims);
	} else if (!strcmp(coarseSolverName, "jacobian2D")){
		multigrid->coarseSolv = &jacobian2D;
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
	bndType *bnd = grid->bnd;
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
		for(int d = 1; d < rank; d++)	subTrueSize[d] = trueSize[d]/(2*q);

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

		long int *subSizeProd = malloc(rank*sizeof(*subSizeProd));
		ailCumProd(subSize, subSizeProd, rank);

		//Alloc slice and val
		double *val = malloc(subSizeProd[rank]*sizeof(*val));
		double *slice = malloc(nSliceMax*sizeof(*slice));

		//Ghost layer vector
		int *subNGhostLayers = malloc(rank*2*sizeof(*subNGhostLayers));
		for(int d = 0; d < 2*rank; d++)	subNGhostLayers[d] = nGhostLayers[d];

		//Copying boundaries
		bndType *subBnd = malloc(rank*2*sizeof(*subBnd));
		for(int d = 0; d < 2*rank; d++)	subBnd[d] = bnd[d];


		//Assign to grid
		Grid *grid = malloc(sizeof(Grid));
		grid->rank = rank;
		grid->size = subSize;
		grid->sizeProd = subSizeProd;
		grid->nGhostLayers = subNGhostLayers;
		grid->bnd = subBnd;
		grid->trueSize = subTrueSize;
		grid->val = val;
		grid->slice = slice;
		grid->h5 = 0;
		grids[q] = grid;
	}

	return grids;
}

/*************************************************
 *		Inline functions
 ************************************************/

 inline static void loopRedBlack2D(double *rhoVal,double *phiVal,long int *sizeProd, int *trueSize, int kEdgeInc,
 				long int g, long int gj, long int gjj, long int gk, long int gkk){

 	gj = g + sizeProd[1];
 	gjj= g - sizeProd[1];
 	gk = g + sizeProd[2];
 	gkk= g - sizeProd[2];

 	for(int k = 1; k < trueSize[2]; k +=2){
 		for(int j = 1; j < trueSize[1]; j += 2){
 			phiVal[g] = 0.25*(	phiVal[gj] + phiVal[gjj] +
 								phiVal[gk] + phiVal[gkk] + rhoVal[g]);
 			g	+=2;
 			gj	+=2;
 			gjj	+=2;
 			gk	+=2;
 			gkk	+=2;
 		}
 		g	+=kEdgeInc;
 		gj	+=kEdgeInc;
 		gjj	+=kEdgeInc;
 		gk	+=kEdgeInc;
 		gkk	+=kEdgeInc;
 	}

 	return;
 }

 inline static void loopRedBlack3D(double *rhoVal,double *phiVal,long int *sizeProd, int *trueSize, int kEdgeInc, int lEdgeInc,
 				long int g, long int gj, long int gjj, long int gk, long int gkk, long int gl, long int gll){

 	gj = g + sizeProd[1];
 	gjj= g - sizeProd[1];
 	gk = g + sizeProd[2];
 	gkk= g - sizeProd[2];
 	gl = g + sizeProd[3];
 	gll= g - sizeProd[3];

 	for(int l = 0; l<trueSize[3]; l+=2){
 		for(int k = 0; k < trueSize[2]; k+=2){
 			for(int j = 0; j < trueSize[1]; j+=2){
 				// msg(STATUS, "g=%d", g);
 				phiVal[g] = 0.125*(phiVal[gj] + phiVal[gjj] +
 								phiVal[gk] + phiVal[gkk] +
 								phiVal[gl] + phiVal[gll] + rhoVal[g]);
 				g	+=2;
 				gj	+=2;
 				gjj	+=2;
 				gk	+=2;
 				gkk	+=2;
 				gl	+=2;
 				gll	+=2;
 			}
 		g	+=kEdgeInc;
 		gj	+=kEdgeInc;
 		gjj	+=kEdgeInc;
 		gk	+=kEdgeInc;
 		gkk	+=kEdgeInc;
 		gl	+=kEdgeInc;
 		gll	+=kEdgeInc;
 		}
 	g	+=lEdgeInc;
 	gj	+=lEdgeInc;
 	gjj	+=lEdgeInc;
 	gk	+=lEdgeInc;
 	gkk	+=lEdgeInc;
 	gl	+=lEdgeInc;
 	gll	+=lEdgeInc;
 	}

 	return;
 }



/*************************************************
 *		DEFINITIONS
 ************************************************/
/*************************************************
 * 		ALLOCATIONS
 * 		DESTRUCTORS
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
	if(nLevels<2) msg(ERROR, "Multi Grid levels is 0, or 1, direct solver not implemented yet \n");


	if(!nMGCycles) msg(ERROR, "MG cycles is 0 \n");


	// Sanity check (true grid points need to be a multiple of 2^(multigrid levels)
	for(int d = 0; d < nDims; d++){
		if(trueSize[d+1] % (int) 2*nLevels){ //Sloppy and wrong
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

	Grid **grids = multigrid->grids;
	int nLevels = multigrid->nLevels;

	for(int n = 1; n < nLevels; n++)
	{
		gFree(grids[n]);
	}
	free(multigrid);

	return;
}

/******************************************************
 *		Iterative Solvers
 *****************************************************/


void jacobian2D(Grid *phi,const Grid *rho, const int nCycles, const  MpiInfo *mpiInfo){

	//Common variables
	int rank = phi->rank;
	long int *sizeProd = phi->sizeProd;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	//Temporary value
	double *tempVal = malloc (sizeProd[rank]*sizeof(*tempVal));

	for(int c = 0; c < nCycles; c++){
		// Index of neighboring nodes
		int gj = sizeProd[1];
		int gjj= -sizeProd[1];
		int gk = sizeProd[2];
		int gkk= -sizeProd[2];

		for(long int g = 0; g < sizeProd[rank]; g++){
			tempVal[g] = 0.25*(	phiVal[gj] + phiVal[gjj] +
								phiVal[gk] + phiVal[gkk] - rhoVal[g]);

			gj++;
			gjj++;
			gk++;
			gkk++;
		}

		for(long int q = 0; q < sizeProd[rank]; q++) phiVal[q] = tempVal[q];

		for(int d = 1; d < rank; d++) gSwapHaloDim(phi, mpiInfo, d);

	}

	return;
}

void jacobian3D(Grid *phi,const Grid *rho, const int nCycles, const  MpiInfo *mpiInfo){

	//Common variables
	int rank = phi->rank;
	long int *sizeProd = phi->sizeProd;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	//Temporary value
	double *tempVal = malloc (sizeProd[rank]*sizeof(*tempVal));
	double coeff = 1./6;

	for(int c = 0; c < nCycles; c++){
		// Index of neighboring nodes
		int gj = sizeProd[1];
		int gjj= -sizeProd[1];
		int gk = sizeProd[2];
		int gkk= -sizeProd[2];
		int gl = sizeProd[3];
		int gll= -sizeProd[3];

		for(long int g = 0; g < sizeProd[rank]; g++){
			tempVal[g] = coeff*(phiVal[gj] + phiVal[gjj] +
								phiVal[gk] + phiVal[gkk] +
								phiVal[gl] + phiVal[gll] +
								- rhoVal[g]);

			gj++;
			gjj++;
			gk++;
			gkk++;
			gl++;
			gll++;
		}

		for(long int q = 0; q < sizeProd[rank]; q++) phiVal[q] = tempVal[q];

		for(int d = 1; d < rank; d++) gSwapHaloDim(phi, mpiInfo, d);

	}

	return;
}



void gaussSeidel2D(Grid *phi, const Grid *rho, int nCycles, const MpiInfo *mpiInfo){

	//Common variables
	int *trueSize = phi->trueSize;
	int *nGhostLayers = phi->nGhostLayers;
	long int *sizeProd = phi->sizeProd;
	int rank = phi->rank;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	//Indexes
	long int g;
	long int gj;
	long int gjj;
	long int gk;
	long int gkk;

	for(int c = 0; c < nCycles;c++){

		//Increments
		int kEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2] + sizeProd[2];

		/**************************
		 *	Red Pass
		 *************************/
		//Odd numbered rows
		g = nGhostLayers[1] + sizeProd[2];
		loopRedBlack2D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, g, gj, gjj, gk, gkk);

		//Even numbered columns
		g = nGhostLayers[1] + 1 + 2*sizeProd[2];
		loopRedBlack2D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, g, gj, gjj, gk, gkk);

		for(int d = 1; d < rank; d++) gSwapHaloDim(phi, mpiInfo, d);
		// gBnd(rho,mpiInfo);
		// gBnd(phi,mpiInfo);


		/***********************************
		 *	Black pass
		 **********************************/
		//Odd numbered rows
		g = nGhostLayers[1] + 1 + sizeProd[2];
		loopRedBlack2D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, g, gj, gjj, gk, gkk);

		//Even numbered columns
		g = nGhostLayers[1] + 2*sizeProd[2];
		loopRedBlack2D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, g, gj, gjj, gk, gkk);


		for(int d = 1; d < rank; d++) gSwapHaloDim(phi, mpiInfo, d);
		// gBnd(rho,mpiInfo);
		// gBnd(phi,mpiInfo);

	}


	return;
}


void gaussSeidel3D(Grid *phi, const Grid *rho, int nCycles, const MpiInfo *mpiInfo){

	//Common variables
	int *trueSize = phi->trueSize;
	int *size = phi->size;
	int *nGhostLayers = phi->nGhostLayers;
	long int *sizeProd = phi->sizeProd;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	//Indexes
	long int g;
	int gj = sizeProd[1];
	int gk = sizeProd[2];
	int gl = sizeProd[3];

	double coeff = 1./6.;


	for(int c = 0; c < nCycles; c++){

		/*********************
		 *	Red Pass
		 ********************/
		g = sizeProd[3]*nGhostLayers[3];
		for(int l = 0; l < trueSize[3];l++){
			for(int k = 0; k < size[2]; k++){
				for(int j = 0; j < size[1]; j+=2){
					phiVal[g] = coeff*(	phiVal[g+gj] + phiVal[g-gj] +
										phiVal[g+gk] + phiVal[g-gk] +
										phiVal[g+gl] + phiVal[g-gl] - rhoVal[g]);
					g	+=2;
				}

				// g+=(-1 + 2*(k%2))*(-1 + 2*(l%2));

				if(l%2){
					if(k%2)	g+=1; else g-=1;
				} else {
					if(k%2) g-=1; else g+=1;
				}

			}
			if(l%2) g-=1; else g+=1;
			// g -= -1 + 2*(l%2);
		}

		gSwapHalo(phi, mpiInfo);

		/*********************
		 *	Black pass
		 ********************/
		 g = sizeProd[1] + sizeProd[3]*nGhostLayers[3];
		 for(int l = 0; l < trueSize[3];l++){
		 	for(int k = 0; k < size[2]; k++){
		 		for(int j = 0; j < size[1]; j+=2){
		 			phiVal[g] = coeff*(	phiVal[g+gj] + phiVal[g-gj] +
		 								phiVal[g+gk] + phiVal[g-gk] +
		 								phiVal[g+gl] + phiVal[g-gl] - rhoVal[g]);

		 			g	+=2;
		 		}
					if(l%2){
						if(k%2)	g-=1; else g+=1;
					} else {
						if(k%2) g+=1; else g-=1;
					}
				// g+=(-1 + 2*(k%2))*(-1 + 2*(l%2));


		 	}
				if(l%2) g+=1; else g-=1;
			// g -= -1 + 2*(l%2);

		 }

		gSwapHalo(phi, mpiInfo);
	}

	return;
}



void gaussSeidel3DNew(Grid *phi, const Grid *rho, int nCycles, const MpiInfo *mpiInfo){

	//Common variables
	int *trueSize = phi->trueSize;
	int *nGhostLayers = phi->nGhostLayers;
	long int *sizeProd = phi->sizeProd;
	int rank = phi->rank;

	//Seperate values
	double *phiVal = phi->val;
	double *rhoVal = rho->val;

	//Indexes
	long int g;
	long int gj;
	long int gjj;
	long int gk;
	long int gkk;
	long int gl;
	long int gll;


	int kEdgeInc = (nGhostLayers[1] + nGhostLayers[rank+1]) + sizeProd[2];
	int lEdgeInc = (nGhostLayers[2] + nGhostLayers[rank+2])*sizeProd[2] + sizeProd[3];


	for(int c = 0; c < nCycles;c++){

		/**************************
		 *	Red Pass
		 *************************/
		//Odd layers - Odd Rows
		g = nGhostLayers[1]*sizeProd[1] + nGhostLayers[2]*sizeProd[2] + nGhostLayers[3]*sizeProd[3];
		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
						g, gj, gjj, gk, gkk, gl, gll);

		//Odd layers - Even Rows
		g = (nGhostLayers[1]+1)*sizeProd[1] + (nGhostLayers[2]+1)*sizeProd[2] + nGhostLayers[3]*sizeProd[3];
		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
						g, gj, gjj, gk, gkk, gl, gll);

		//Even layers - Odd Rows
		g = (nGhostLayers[1])*sizeProd[1] + (nGhostLayers[2])*sizeProd[2] + (nGhostLayers[3]+1)*sizeProd[3];
		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
						g, gj, gjj, gk, gkk, gl, gll);

		//Even layers - Even Rows
		g = (nGhostLayers[1] + 1)*sizeProd[1] + (nGhostLayers[2]+1)*sizeProd[2] + (nGhostLayers[3]+1)*sizeProd[3];
		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
						g, gj, gjj, gk, gkk, gl, gll);

		gSwapHalo(phi, mpiInfo);

		/***********************************
		 *	Black pass
		 **********************************/
		 //Odd layers - Odd Rows
 		g = (nGhostLayers[1]*sizeProd[1]+1) + nGhostLayers[2]*sizeProd[2] + nGhostLayers[3]*sizeProd[3];
 		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
 						g, gj, gjj, gk, gkk, gl, gll);

 		//Odd layers - Even Rows
 		g = (nGhostLayers[1])*sizeProd[1] + (nGhostLayers[2]+1)*sizeProd[2] + nGhostLayers[3]*sizeProd[3];
 		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
 						g, gj, gjj, gk, gkk, gl, gll);

 		//Even layers - Odd Rows
 		g = (nGhostLayers[1])*sizeProd[1] + (nGhostLayers[2])*sizeProd[2] + (nGhostLayers[3]+1)*sizeProd[3];
 		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
 						g, gj, gjj, gk, gkk, gl, gll);

		//Even layers - Even Rows
 		g = (nGhostLayers[1]+1)*sizeProd[1] + (nGhostLayers[2]+1)*sizeProd[2] + (nGhostLayers[3]+1)*sizeProd[3];
 		loopRedBlack3D(rhoVal, phiVal, sizeProd, trueSize, kEdgeInc, lEdgeInc,
 						g, gj, gjj, gk, gkk, gl, gll);

		gSwapHalo(phi, mpiInfo);
	}


	return;
}

/***********************************************************
 *			RESTRICTORS/PROLONGATORS
 **********************************************************/

void halfWeightRestrict3D(const Grid *fine, Grid *coarse){

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;
	int *cTrueSize = coarse->trueSize;


	//Indexes
	long int c = cSizeProd[1]*nGhostLayers[1] + cSizeProd[2]*nGhostLayers[2] + cSizeProd[3]*nGhostLayers[3];

	long int f = fSizeProd[1]*nGhostLayers[1] + fSizeProd[2]*nGhostLayers[2] + fSizeProd[3]*nGhostLayers[3];
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

	double coeff = 1./12.;


	//Cycle Coarse grid
	for(int l = 0; l<cTrueSize[3]; l++){
		for(int k = 0; k < cTrueSize[2]; k++){
			for(int j = 0; j < cTrueSize[1]; j++){
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

void bilinearProlong3D(Grid *fine, const Grid *coarse,const  MpiInfo *mpiInfo){

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int *fSize = fine->size;
	int *fTrueSize =fine->trueSize;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;
	int *cTrueSize = coarse->trueSize;

	//Help Indexes
	long int f = fSizeProd[1] + fSizeProd[2] + fSizeProd[3];
	long int c = cSizeProd[1] + cSizeProd[2] + cSizeProd[3];
	long int fNext;
	long int fPrev;

	//Edge jumps
	int cKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2];
	int fKEdgeInc = cKEdgeInc + fSizeProd[2];
	int cLEdgeInc = (nGhostLayers[3] + nGhostLayers[rank + 3])*cSizeProd[2];
	int fLEdgeInc = (nGhostLayers[3] + nGhostLayers[rank + 3])*fSizeProd[2] + fSizeProd[3];

	//Direct insertion c->f
	for(int l = 0; l < cTrueSize[3]; l++){
		for(int k = 0; k < cTrueSize[2]; k++){
			for(int j = 0; j < cTrueSize[1]; j++){
				fVal[f] = cVal[c];
				c++;
				f+=2;
			}
			c+= cKEdgeInc;
			f+= fKEdgeInc;
		}
		c+= cLEdgeInc;
		f+= fLEdgeInc;
	}

	//Filling ghostlayer
	gSwapHaloDim(fine, mpiInfo, 3);

	//Interpolation 3rd Dim
	f = fSizeProd[1] + fSizeProd[2] + 2*fSizeProd[3];
	fNext = f + fSizeProd[3];
	fPrev = f - fSizeProd[3];

	for(int l = 0; l < fTrueSize[3]; l+=2){
		for(int k = 0; k < fSize[2]; k+=2){
			for(int j = 0; j < fSize[1]; j+=2){
				fVal[f] = 0.5*(fVal[fPrev]+fVal[fNext]);
				f +=2;
				fNext +=2;
				fPrev +=2;
			}
			f		+=fSizeProd[2];
			fNext 	+=fSizeProd[2];
			fPrev 	+=fSizeProd[2];
		}
		f		+=fSizeProd[3];
		fNext 	+=fSizeProd[3];
		fPrev 	+=fSizeProd[3];
	}

	gSwapHaloDim(fine, mpiInfo, 2);

	//Interpolation 2nd Dim
	f = fSizeProd[1] + 2*fSizeProd[2] + fSizeProd[3];
	fNext = f + fSizeProd[2];
	fPrev = f - fSizeProd[2];

	for(int l = 0; l < fTrueSize[3]; l++){
		for(int k = 0; k < fSize[2]; k+=2){
			for(int j = 0; j < fSize[1]; j+=2){
				fVal[f] = 0.5*(fVal[fPrev]+fVal[fNext]);
				f +=2;
				fNext +=2;
				fPrev +=2;
			}
			f		+=fSizeProd[2];
			fNext 	+=fSizeProd[2];
			fPrev 	+=fSizeProd[2];
		}
	}

	gSwapHaloDim(fine, mpiInfo, 1);

	//Interpolation 2nd Dim
	f = 2*fSizeProd[1] + fSizeProd[2] + fSizeProd[3];
	fNext = f + fSizeProd[1];
	fPrev = f - fSizeProd[1];

	for(int l = 0; l < fTrueSize[3]; l++){
		for(int k = 0; k < fTrueSize[2]; k++){
			for(int j = 0; j < fSize[1]; j+=2){
				fVal[f] = 0.5*(fVal[fPrev]+fVal[fNext]);
				f +=2;
				fNext +=2;
				fPrev +=2;
			}
		}
		f		+=2*fSizeProd[2];
		fNext 	+=2*fSizeProd[2];
		fPrev 	+=2*fSizeProd[2];
	}


	return;
}


void bilinearProlong2D(Grid *fine, const Grid *coarse, const MpiInfo *mpiInfo){

	//Load fine grid
	double *fVal = fine->val;
	long int *fSizeProd = fine->sizeProd;
	int *fSize = fine->size;
	int rank = fine->rank;
	int *nGhostLayers = fine->nGhostLayers;

	//Load coarse grid
	double *cVal = coarse->val;
	long int *cSizeProd = coarse->sizeProd;
	int *cSize = coarse->size;

	//Help Indexes
	long int f = fSizeProd[2] + fSizeProd[1];
	long int c = cSizeProd[2] + cSizeProd[1];
	long int fNext;
	long int fPrev;

	int cKEdgeInc = nGhostLayers[2] + nGhostLayers[rank + 2];
	int fKEdgeInc = cKEdgeInc + fSizeProd[2];

	//Direct insertion c->f
	for(int k = nGhostLayers[2]; k < cSize[2]-nGhostLayers[rank + 2]; k++){
		for(int j = nGhostLayers[1]; j < cSize[1]-nGhostLayers[rank+1]; j++){
			fVal[f] = cVal[c];
			c++;
			f+=2;
		}
		c+= cKEdgeInc;
		f+= fKEdgeInc;
	}

	//Filling ghost cells
	gSwapHaloDim(fine, mpiInfo, 2);

 	f= fSizeProd[1];
	fNext = f + fSizeProd[2];
	fPrev = f - fSizeProd[2];
	//Odd numbered columns, interpolating vertically
	for(int k = 0; k < fSize[2]; k+=2){
		for(int j = 0; j < fSize[1]; j+=2){
			fVal[f] += 0.5*(fVal[fPrev]+fVal[fNext]);
			f 		+= 2;
			fNext 	+= 2;
			fPrev 	+= 2;
		}
		f		+= fSizeProd[2];
		fNext 	+= fSizeProd[2];
		fPrev 	+= fSizeProd[2];
	}

	//Filling ghost cells
	gSwapHaloDim(fine, mpiInfo, 1);

	//Even numbered columns, interpolating horizontally
	f = 0;
	fNext = f + fSizeProd[1];
	fPrev = f - fSizeProd[1];

	for(int k = 0; k < fSize[2]; k+=1){
		for(int j = 0; j < fSize[1]; j+=2){
			fVal[f] += 0.5*(fVal[fPrev]+fVal[fNext]);
			f 		+= 2;
			fNext 	+= 2;
			fPrev 	+= 2;
		}
	}

	return;
}

/*******************************************************
 *			VARIOUS COMPUTATIONS (RESIDUAL)
 ******************************************************/

void mgResidual(Grid *res, const Grid *rho, const Grid *phi,const MpiInfo *mpiInfo){

	//Load
	long int *sizeProd = res->sizeProd;
	int rank = res->rank;
	double *resVal = res->val;
	double *rhoVal = rho->val;

	gFinDiff2nd2D(res, phi, mpiInfo);

	for (long int g = 0; g < sizeProd[rank]; g++) resVal[g] -= rhoVal[g];

	for(int d = 1; d < rank; d++) 	gSwapHaloDim(res, mpiInfo, d);

	return;
}


double mgResMass3D(Grid *grid){

	//Load
	int rank = grid->rank;
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	int *nGhostLayers = grid->nGhostLayers;
	double *val = grid->val;

	long int ind;
	double mass = 0;

	for(int l = nGhostLayers[3]; l < size[3]-nGhostLayers[rank+3]; l++){
		for(int k = nGhostLayers[2]; k < size[2]-nGhostLayers[rank+2]; k++){
			for(int j = nGhostLayers[1]; j < size[1]-nGhostLayers[rank+1]; j++){
				ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				mass += val[ind]*val[ind];
			}
		}
	}



	return mass;
}

/*****************************************************
 *			MG CYCLES
 ****************************************************/

void inline static mgVRecursive(int level, int targetLvl, Multigrid *mgRho, Multigrid *mgPhi,
 									Multigrid *mgRes, const MpiInfo *mpiInfo){

	//Solve and return at coarsest level
	if(level == targetLvl){
		gSwapHalo(mgPhi->grids[level], mpiInfo);
		mgRho->coarseSolv(mgPhi->grids[level], mgRho->grids[level], mgRho->nCoarseSolve, mpiInfo);
		mgRho->prolongator(mgRes->grids[level-1], mgPhi->grids[level], mpiInfo);
		return;
	}
	// msg(STATUS, "rank = %d, lvl = %d, targetlvl = %d", mpiInfo->mpiRank, level, targetLvl);


	//Gathering info
	int nPreSmooth = mgRho->nPreSmooth;
	int nPostSmooth= mgRho->nPostSmooth;
	Grid *phi = mgPhi->grids[level];
	Grid *rho = mgRho->grids[level];
	Grid *res = mgRes->grids[level];

	// msg(STATUS, "Got here");
	//Prepare to go down
	gSwapHalo(rho, mpiInfo);
	gBnd(rho,mpiInfo);
	mgRho->preSmooth(phi, rho, nPreSmooth, mpiInfo);
	// msg(STATUS, "Passed bnd?");

	mgResidual(res, rho, phi, mpiInfo);

	gSwapHalo(res, mpiInfo);
	mgRho->restrictor(res, mgRho->grids[level + 1]);

	//Go coarser and solve
	mgVRecursive(level + 1, targetLvl, mgRho, mgPhi, mgRes, mpiInfo);

	//Prepare to go up
	gAddTo( phi, res );

	gSwapHalo(phi,mpiInfo);
	gBnd(phi,mpiInfo);
	mgRho->postSmooth(phi, rho, nPostSmooth, mpiInfo);

	if(level > 0){
		mgRho->prolongator(mgRes->grids[level-1], phi, mpiInfo);
	}
	return;
}

void mgVRegular(int level,int targetLvl, Multigrid *mgRho, Multigrid *mgPhi,
 									Multigrid *mgRes, const MpiInfo *mpiInfo){

	//Gathering info
	int nPreSmooth = mgRho->nPreSmooth;
	int nPostSmooth= mgRho->nPostSmooth;
	int nCoarseSolv= mgRho->nCoarseSolve;


	//Down to coarsest level
	for(int level = 0; level <targetLvl; level ++){
		//Load grids
		Grid *phi = mgPhi->grids[level];
		Grid *rho = mgRho->grids[level];
		Grid *res = mgRes->grids[level];

		gSwapHalo(phi,mpiInfo);
		gBnd(phi,mpiInfo);
		mgRho->preSmooth(phi, rho, nPreSmooth, mpiInfo);
		mgResidual(res, rho, phi, mpiInfo);
		mgRho->restrictor(res, mgRho->grids[level + 1]);
	}

	//Solve at coarsest
	mgRho->coarseSolv(mgPhi->grids[targetLvl], mgRho->grids[targetLvl], nCoarseSolv, mpiInfo);
	mgRho->prolongator(mgRes->grids[targetLvl-1], mgPhi->grids[targetLvl], mpiInfo);

	//Up to finest
	for(int level = targetLvl-1; level >-1; level --){
		Grid *phi = mgPhi->grids[level];
		Grid *rho = mgRho->grids[level];
		Grid *res = mgRho->grids[level];


		//Prepare to go up
		gAddTo( phi, res );
		mgRho->postSmooth(phi, rho, nPostSmooth, mpiInfo);
		if(level > 0)	mgRho->prolongator(mgRes->grids[level-1], phi, mpiInfo);
	}

	return;

}


void linearMGSolv(Multigrid *mgRho, Multigrid *mgPhi, Multigrid *mgRes, const MpiInfo *mpiInfo){

	int nMGCycles = mgRho->nMGCycles;
	int targetLvl = mgRho->nLevels-1;

	gZero(mgPhi->grids[0]);
	for(int c = 0; c < nMGCycles; c++){
		mgVRecursive(0,targetLvl, mgRho, mgPhi, mgRes, mpiInfo);
		// mgVRegular(0, targetLvl, mgRho, mgPhi, mgRes, mpiInfo);
	}

	return;
}
