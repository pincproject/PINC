/**
 * @file	    main.test.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 * @copyright   University of Oslo, Norway
 * @brief	    Unit testing main routine
 * @date        16.12.15
 */

#include "test.h"
#include "iniparser.h"
#include <mpi.h>

int main(int argc, char *argv[]){

	MPI_Init(&argc,&argv);

	testAux();
	testIo();
	testGrid();
	testPopulation();
	testPusher();
	testMultigrid();
	utSummary();

	MPI_Finalize();

	return 0;
}
