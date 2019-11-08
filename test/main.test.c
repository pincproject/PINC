/**
 * @file	    main.test.c
 * @brief	    Unit testing main routine
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 */

#include "test.h"
#include "iniparser.h"
#include <mpi.h>

int main(int argc, char *argv[]){

	MPI_Init(&argc,&argv);
	iniSetDummy(argc,argv);

/* 	testAux();
	testIo();
	testGrid();
	testPopulation();
	testPusher();
	testMultigrid(); */
	testObject();
	utSummary();

	MPI_Finalize();

	return 0;
}
