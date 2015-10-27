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

	/*
	 * SUCCESSFUL EXIT
	 */

	freePopulation(pop);

	iniparser_freedict(ini);
	msg(STATUS,"PINC completed successfully!");

	MPI_Finalize();
	return 0;

}


