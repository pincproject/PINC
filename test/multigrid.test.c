/**
 * @file		multigrid.test.c
 * @author		Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for multigrid.c
 * @date		17.12.15
 */

#include "pinc.h"
#include "test.h"
#include "iniparser.h"

static int testStructs(){
	// iniparser_load("test/test.ini");
	//
	//
	// dictionary *ini = malloc(sizeof(*ini));

	// Grid *grid = gAlloc(ini, 1);
	// for(int i = 0; i < 5; i++){
	// 	grid->val[i] = 1.;
	// }
	//
	// //Changing them from the multigrid struct
	// for(int i = 0; i < 5; i++){
	// 	multigrid->grids[0]->val[i] = 5.;
	// }
	//
	//
	// for(int i = 0; i < 5; i++){
	// 	grid->val[i] = 1.;
	// }

	int testCase = 5;
	utAssert(testCase == 6, "Failed test: 5 is not equal to 6");
	return 0;
}

// All tests for grid.c is contained in this function
void testMultigrid(){
	utRun(&testStructs);
}
