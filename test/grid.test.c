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


int testGValDebug(){

	// dictionary *ini = iniOpenDummy();
	// dictionary *ini = iniGetDummy();

	// double debye = iniparser_getdouble(ini,"grid:debye",0);
	// printf("debye=%f\n",debye);

	return 0;
	// Grid *grid = gAlloc(ini,3);
	// MpiInfo *mpiInfo = gAllocMpi(ini);
	// free(ini);
	//
	// gValDebug(grid,mpiInfo);
	//
	// // double *val = grid->val;
	//
	// for(int i=0;i<60;i++){
	// 	// printf("%f\n",val[i]);
	// }
	//
	// return 0;

}

static int testSwapHalo(){
	//Not written because several comp nodes are not run in test mode
	//
	// dictionary *ini = iniGetDummy();
	//
	// iniparser_set(ini, "grid:trueSize", "6,6");
	// iniparser_set(ini, "grid:stepSize", "1,1");
	// iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1");
	//
	// Grid *grid = gAlloc(ini, 1);
	//
	// //Load
	// int *size = val->size;
	// int *trueSize = val->trueSize;
	// long int *sizeProd = val->sizeProd;
	// int *nGhostLayers = val->nGhostLayers;
	//
	// //Fill grid
	// for(int q = 0; q < sizeProd[rank]; q++) val[q] = 0.0;
	//
	// //Fill boundary layers
	//

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

	return 0;
}

// All tests for grid.c is contained in this function
void testGrid(){
	utRun(&testGValDebug);
	utRun(&testSwapHalo);
	utRun(&testFinDiff1st);
	utRun(&testFinDiff2nd2D);
	utRun(&testgFinDiff2nd3D);

}
