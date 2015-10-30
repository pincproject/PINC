/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 * Replaces old DiP3D main.c file by Wojciech Jacek Miloch.
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include "multigrid.h"


int main(int argc, char *argv[]){

	MPI_Init(&argc,&argv);
	msg(STATUS,"PINC started.");

	/*
	 * READ INPUT FILE
	 */
	 
	dictionary *ini = iniOpen(argc,argv);
//	ini_complete_time(ini);
//	ini_complete_grid(ini);

	Population *pop = allocPopulation(ini);

	Grid *grid = allocGrid(ini, 1);

	msg(STATUS, "#Grid values per grid point = %d", grid->nValues);

	/*
	 * SUCCESSFUL EXIT
	 */



	freePopulation(pop);
	freeGrid(grid);

	iniparser_freedict(ini);
	msg(STATUS,"PINC completed successfully!");

	MPI_Finalize();
	return 0;

}
