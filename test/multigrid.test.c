/**
 * @file		multigrid.test.c
 * @author		Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for multigrid.c
 * @date		17.12.15
 */

#include "pinc.h"
#include "test.h"

static int testStructs(){

	int testCase = 5;
	utAssert(testCase == 6, "Hello");
	return 0;
}

// All tests for grid.c is contained in this function
void testMultigrid(){
	utRun(&testStructs);
	printf("Hello\n");
}
