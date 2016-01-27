/**
 * @file		multigrid.test.c
 * @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for multigrid.c
 * @date		17.12.15
 */

#include "pinc.h"
#include "test.h"
#include "multigrid.h"
#include "iniparser.h"


static int testStructs(){
	/*
	 * Checks that val is accessible and changeable from
	 * multigrid and grid
	*/

	dictionary *ini = iniGetDummy();

	iniparser_set(ini, "grid:trueSize", "12,12,12");
	iniparser_set(ini, "grid:stepSize", "1,1,1");
	iniparser_set(ini, "grid:nGhostLayers", "0,0,0,0,0,0");


	Grid *grid = gAlloc(ini, 1);

	Multigrid *multigrid = mgAlloc(ini, grid);
	iniparser_freedict(ini);


	long int *sizeProd = grid->sizeProd;
	int rank = grid->rank;

	for(int i = 0; i < sizeProd[rank]; i++){
		grid->val[i] = 1.;
	}

	//Accesible from MG
	utAssert(multigrid->grids[0]->val[10]-1. < 0.001, "Data not accesible from MG");

	//Changing them from the multigrid struct
	for(int i = 0; i < sizeProd[rank]; i++){
		multigrid->grids[0]->val[i] = 5.;
	}

	utAssert(grid->val[10]-5. < 0.001, "Grid data not changed from MG");

	gFree(grid);
	// mgFree(multigrid);


	return 0;
}

static int testGaussSeidel(){

	return 0;
}

static int testRestrictor(){
	/*
	 * Set up a predefined fine grid, then checks the restrictor against a
	 * precalculated results
	 */
	//TBF(inished)

	dictionary *ini = iniGetDummy();

	iniparser_set(ini, "grid:trueSize", "4,4");
	iniparser_set(ini, "grid:stepSize", "1,1");
	iniparser_set(ini, "grid:nGhostLayers", "1,1,1,1");
	iniparser_set(ini, "multigrid:restrictor", "halfWeight");
	iniparser_set(ini, "multigrid:mgLevels", "1");
	iniparser_set(ini, "msgfiles:parsefile", "parsedump.txt");

	Grid *fine = gAlloc(ini, 1);
	Multigrid *multigrid = mgAlloc(ini, fine);
	Grid *coarse = multigrid->grids[1];

	iniparser_freedict(ini);

	// multigrid->restrictor(fine, coarse);

	//Populate and dump fine grid indexes
	long int *fSizeProd = fine->sizeProd;
	int rank = fine->rank;
	double *fVal = fine->val;

	for(long int g = 0; g < fSizeProd[rank]; g++){
		fVal[g] = (double) g;
	}

	dumpWholeGrid(ini, fine);

	//Populate and dump coarse grid indexes
	long int *cSizeProd = coarse->sizeProd;
	double *cVal = coarse->val;

	for(long int g = 0; g < 5; g++){
		// cVal[g] = (double) g;
	}

	// dumpWholeGrid(ini, coarse);



	return 0;
}

// All tests for grid.c is contained in this function
void testMultigrid(){
	utRun(&testStructs);
	utRun(&testGaussSeidel);
	utRun(&testRestrictor);
}
