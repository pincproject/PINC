/**
 * @file		grid.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for grid.c
 * @date		17.12.15
 */

#include "pinc.h"
#include "test.h"
#include "iniparser.h"
#include <math.h>

// Loads the 5x4x3 grid which is used as example several places, and sets
// elements equal to linear index p
Grid *loadGrid543(int nValues){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"grid:trueSize","5,4,3");
	iniparser_set(ini,"grid:stepSize","1,1,1");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");

	Grid *grid = gAlloc(ini,nValues);
	long int nElements = grid->sizeProd[grid->rank];
	for(int p=0;p<nElements;p++) grid->val[p] = p;

	iniparser_freedict(ini);
	return grid;

}

static int testGAlloc(){

	Grid *grid = loadGrid543(2);
	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	int *nGhostLayers = grid->nGhostLayers;
	double *stepSize = grid->stepSize;

	int expectedSize[] = {2,5,4,3};
	long int expectedSizeProd[] = {1,2,10,40,120};
	int expectedNGhostLayers[] = {0,0,0,0,0,0};
	double expectedStepSize[] = {1,1,1,1};

	utAssert(grid->rank==4,"wrong rank assigned by gAlloc");
	utAssert(aiEq(size        ,expectedSize        ,4),"wrong size assigned by gAlloc");
	utAssert(alEq(sizeProd    ,expectedSizeProd    ,5),"wrong sizeProd assigned by gAlloc");
	utAssert(aiEq(nGhostLayers,expectedNGhostLayers,6),"wrong nGhostLayers assigned by gAlloc");
	utAssert(adEq(stepSize    ,expectedStepSize    ,4,0),"wrong stepSize assigned by gAlloc");

	return 0;

}

static int testSwapHalo(){
	//Not written because several comp nodes are not run in test mode,
	//instead done in main.local.c


	return 0;
}

static int testFinDiff1st(){
	//Tests F(x,y,z) = x^2 - z -> d/dx F = 2x, d/dy F = 0 and d/dz = -1
	dictionary *ini = iniGetDummy();

	int testResult = 0;

	iniparser_set(ini, "grid:trueSize", "12,12,12");
	iniparser_set(ini, "grid:stepSize", "1,1,1");
	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1,1,1");

	Grid *phi = gAlloc(ini, 1);
	Grid *E = gAlloc(ini, 3);

	iniparser_freedict(ini);

	//Load
	int *size = phi->size;
	int *trueSize = phi->trueSize;
	long int *sizeProd = phi->sizeProd;
	long int *fieldSizeProd = E->sizeProd;
	int *nGhostLayers = phi->nGhostLayers;

	double *phiVal = phi->val;
	double *eVal = E->val;

	int ind;
	for(int j = 0; j < size[1]; j++){
		for (int k = 0; k<size[2]; k++) {
			for(int l = 0; l<size[3]; l++){
				ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				phiVal[ind] = (double) (j*j-l);
			}
		}
	}

	gFinDiff1st(phi, E);

	//Test d/dx
	double answer;
	//Testing the inner grid
	for(int j = nGhostLayers[1]; j < trueSize[1]; j++){
		for(int k = nGhostLayers[2]; k < trueSize[2]; k++){
			for(int l = nGhostLayers[3]; l < trueSize[3]; l++){
				answer = (double) 2*j;
				ind = j*fieldSizeProd[1] + k*fieldSizeProd[2] + l*fieldSizeProd[3];
				if((eVal[ind]-answer < -0.001) || (eVal[ind] -answer > 0.001)) testResult = 1;
			}
		}
	}
	utAssert(testResult==0, "Failed centered Finite Difference: x Dim");
	testResult = 0;
	//Test d/dy
	for(int j = nGhostLayers[1]; j < trueSize[1]; j++){
		for(int k = nGhostLayers[2]; k < trueSize[2]; k++){
			for(int l = nGhostLayers[3]; l < trueSize[3]; l++){
				answer = (double) 0;
				ind = 1 + j*fieldSizeProd[1] + k*fieldSizeProd[2] + l*fieldSizeProd[3];
				if((eVal[ind]-answer < -0.001) || (eVal[ind] -answer > 0.001)) testResult = 1;
			}
		}
	}
	utAssert(testResult==0, "Failed centered Finite Difference: y Dim");
	testResult = 0;
	//Test d/dy
	for(int j = nGhostLayers[1]; j < trueSize[1]; j++){
		for(int k = nGhostLayers[2]; k < trueSize[2]; k++){
			for(int l = nGhostLayers[3]; l < trueSize[3]; l++){
				answer = (double) -1;
				ind = 2 + j*fieldSizeProd[1] + k*fieldSizeProd[2] + l*fieldSizeProd[3];
				if((eVal[ind]-answer < -0.001) || (eVal[ind] -answer > 0.001)){
					testResult = 1;
				}
			}
		}
	}
	utAssert(testResult==0, "Failed centered Finite Difference: z Dim");

	gFree(phi);
	gFree(E);

	return 0;
}


static int testFinDiff2nd2D(){
	/*
	 * Tests d^2/dx^2 by setting the grid to f(x,y)=(3*x^2 - 0.5y^2) ignoring ghost nodes
	 * The laplacan(f(x,y)) should then be 5 over the inner grid, where the operation
	 * is applied
	 */
	dictionary *ini = iniGetDummy();

	int testResult = 0;

	iniparser_set(ini, "grid:trueSize", "12,12");
	iniparser_set(ini, "grid:stepSize", "1,1");
	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1");

	Grid *rho = gAlloc(ini, 1);
	Grid *phi = gAlloc(ini, 1);

	iniparser_freedict(ini);


	int *size = rho->size;
	int *trueSize = rho->trueSize;
	long int *sizeProd = rho->sizeProd;
	int *nGhostLayers = rho->nGhostLayers;

	double *rhoVal = rho->val;
	double *phiVal = phi->val;

	int ind;
	for(int j = 0; j < size[1]; j++){
		for (int k = 0; k<size[2]; k++) {
			ind = j*sizeProd[1] + k*sizeProd[2];// + l*sizeProd[3];
			rhoVal[ind] = (double) (3*j*j - 0.5*k*k);
		}
	}

	gFinDiff2nd2D(phi,rho);

	//Testing the inner grid
	double ans = 5.;
	for(int j = nGhostLayers[1]; j < trueSize[1]; j++){
		for(int k = nGhostLayers[2]; k < trueSize[2]; k++){
			ind = j*sizeProd[1] + k*sizeProd[2];
			if((phiVal[ind]-ans > 0.01)|| (phiVal[ind]-ans < -0.01)) testResult = 1;
		}
	}

	utAssert(testResult==0, "Failed laplacian 2D");

	testResult = 0;

	gFree(rho);
	gFree(phi);

	return 0;
}

static int testgFinDiff2nd3D(){

	/*
	 * 3D Test laplacian(x - 2y^2 + 4z^3)
	 */

	dictionary *ini = iniGetDummy();

 	int testResult = 0;

	iniparser_set(ini, "grid:trueSize", "12,12,12");
	iniparser_set(ini, "grid:stepSize", "1,1,1");
	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1,1,1");

	Grid *rho = gAlloc(ini, 1);
	Grid *phi = gAlloc(ini, 1);

	iniparser_freedict(ini);

	int *size = rho->size;
	int *trueSize = rho->trueSize;
	long int *sizeProd = rho->sizeProd;
	int *nGhostLayers = rho->nGhostLayers;

	double *rhoVal = rho->val;
	double *phiVal = phi->val;

	int ind;

	for(int j = 0; j < size[1]; j++){
		for (int k = 0; k<size[2]; k++) {
			for(int l = 0; l < size[3]; l++){
				ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				rhoVal[ind] = (double) (j - 2.*k*k + 4*l*l*l);
			}
	 	}
	}

	gFinDiff2nd3D(phi,rho);

	double answer;
	//Testing the inner grid
	for(int j = nGhostLayers[1]; j < trueSize[1]; j++){
		 for(int k = nGhostLayers[2]; k < trueSize[2]; k++){
			for(int l = nGhostLayers[3]; l < trueSize[2]; l++){
				ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
				answer = (double) -4 + 24*l;
   			 	if((phiVal[ind]-answer > 0.001) || (phiVal[ind]-answer < -0.001))
 					testResult = 1;
			}
		}
	}

	utAssert(testResult==0, "Failed laplacian 3D");

	gFree(rho);
	gFree(phi);

	return;
}

static int testGCreateNeighborhood(){

	dictionary *ini = iniGetDummy();

	iniparser_set(ini,"grid:trueSize","10,11,12");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");
	iniparser_set(ini,"grid:thresholds","1.5,1.5,1.0,-1.5,-1.5,-1.0");
	iniparser_set(ini,"grid:nEmigrantsAlloc","1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27");

	Grid *grid = gAlloc(ini,1);
	MpiInfo *mpiInfo = gAllocMpi(ini);

	aiSet(mpiInfo->nSubdomains,3,5,4,3);
	aiSet(mpiInfo->subdomain,3,3,1,1);
	aiCumProd(mpiInfo->nSubdomains,mpiInfo->nSubdomainsProd,3);

	mpiInfo->mpiRank = aiDotProd(mpiInfo->subdomain,mpiInfo->nSubdomainsProd,3);
	mpiInfo->mpiSize = mpiInfo->nSubdomainsProd[mpiInfo->nDims];

	gCreateNeighborhood(ini,mpiInfo,grid);

	utAssert(mpiInfo->nNeighbors==27, "Wrong number of neighbours: %i",mpiInfo->nNeighbors);
	utAssert(mpiInfo->neighborhoodCenter==13, "Wrong center of tha hood computed: %i",mpiInfo->neighborhoodCenter);

	double *result = malloc(6*sizeof(*result));
	adSet(result,6,1.5,1.5,1.0,8.5,9.5,11.0);
	utAssert(adEq(mpiInfo->thresholds,result,6,pow(10,-15)), "Wrong thresholds assigned (using negative values)");

	long int *resultl = malloc(27*sizeof(*resultl));
	alSet(resultl,27,1,2,3,4,5,6,7,8,9,10,11,12,13,0,15,16,17,18,19,20,21,22,23,24,25,26,27);
	utAssert(alEq(mpiInfo->nEmigrantsAlloc,resultl,27), "Wrong number of migrants allocated for (full specification)");

	gDestroyNeighborhood(mpiInfo);
	iniparser_set(ini,"grid:nEmigrantsAlloc","1,2,3");
	gCreateNeighborhood(ini,mpiInfo,grid);

	alSet(resultl,27,1,2,1,2,3,2,1,2,1,2,3,2,3,0,3,2,3,2,1,2,1,2,3,2,1,2,1);
	utAssert(alEq(mpiInfo->nEmigrantsAlloc,resultl,27), "Wrong number of migrants allocated for (smart specification)");

	gDestroyNeighborhood(mpiInfo);
	iniparser_set(ini,"grid:nEmigrantsAlloc","4");
	gCreateNeighborhood(ini,mpiInfo,grid);

	alSetAll(resultl,27,4);
	resultl[13]=0;
	utAssert(alEq(mpiInfo->nEmigrantsAlloc,resultl,27), "Wrong number of migrants allocated for (all equal specification)");

	return 0;
}

// All tests for grid.c is contained in this function
void testGrid(){

	utRun(&testGValDebug);
	utRun(&testSwapHalo);
	utRun(&testFinDiff1st);
	utRun(&testFinDiff2nd2D);
	utRun(&testgFinDiff2nd3D);
	utRun(&testGAlloc);
	utRun(&testGCreateNeighborhood);

}
