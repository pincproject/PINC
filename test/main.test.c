/**
 * @file	    main.test.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 * @copyright   University of Oslo, Norway
 * @brief	    Unit testing main routine
 * @date        16.12.15
 */

#include "test.h"

int main(int argc, char *argv[]){

	testAux();
	testIo();
	testGrid();
	testPopulation();
	testPusher();
	testMultigrid();
	utSummary();

	return 0;
}
