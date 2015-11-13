/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 * Replaces old DiP3D main.c file by Wojciech Jacek Miloch.
 */

#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include "multigrid.h"
#include "test.h"

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	msg(STATUS,"PINC started.");

	/*
	 * INITIALIZE PINC VARIABLES
	 */

	dictionary *ini = iniOpen(argc,argv);

	Grid *grid = allocGrid(ini);
	GridQuantity *gridQuantity = allocGridQuantity(ini, grid, 1);
	Multigrid *multigrid = allocMultigrid(ini, gridQuantity);


	/*
	 * 	Test Area
	 */

/*	testGridAndMGStructs(ini, gridQuantity, multigrid);
*/	testBoundarySendRecieve(ini, gridQuantity, multigrid);

	/*
	 * FINALIZE PINC VARIABLES
	 */
	freeGrid(grid);
	freeGridQuantity(gridQuantity);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
