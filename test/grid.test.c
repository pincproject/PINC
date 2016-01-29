/**
 * @file		grid.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
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
	utAssert(aEq((char*)size        ,(char*)expectedSize        ,4*sizeof(*size)        ),"wrong size assigned by gAlloc");
	utAssert(aEq((char*)sizeProd    ,(char*)expectedSizeProd    ,5*sizeof(*sizeProd)    ),"wrong sizeProd assigned by gAlloc");
	utAssert(aEq((char*)nGhostLayers,(char*)expectedNGhostLayers,6*sizeof(*nGhostLayers)),"wrong nGhostLayers assigned by gAlloc");
	utAssert(aEq((char*)stepSize    ,(char*)expectedStepSize    ,4*sizeof(*stepSize)    ),"wrong stepSize assigned by gAlloc");

	return 0;

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
	free(mpiInfo->nSubdomainsProd);
	mpiInfo->nSubdomainsProd = aiCumProd(mpiInfo->nSubdomains,3);

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
	utRun(&testGAlloc);
	utRun(&testGCreateNeighborhood);
}
