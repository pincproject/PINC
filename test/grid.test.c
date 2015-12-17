/**
 * @file		grid.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for grid.c
 * @date		17.12.15
 */

#include "pinc.h"
#include "test.h"
#include "iniparser.h"

dictionary *iniOpenDummy(){
	return iniparser_load("test/test.ini");

	dictionary *ini = malloc(sizeof(*ini));
	iniparser_set(ini,"grid:trueSize","3,3,3");

}

int testGValDebug(){

	dictionary *ini = iniOpenDummy();
	Grid *grid = gAlloc(ini,3);
	MpiInfo *mpiInfo = gAllocMpi(ini);
	free(ini);

	gValDebug(grid,mpiInfo);

	// double *val = grid->val;

	for(int i=0;i<60;i++){
		// printf("%f\n",val[i]);
	}

	return 0;

}

// All tests for grid.c is contained in this function
void testGrid(){
	utRun(&testGValDebug);
}
