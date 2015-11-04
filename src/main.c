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
	gridParseDump(ini, grid);

	GridQuantity *gridQuantity = allocGridQuantity(ini, grid, 1);

	Multigrid *multigrid = allocMultigrid(ini, gridQuantity);

/*
	msg(STATUS, "#Grid values per grid point = %d", grid->nValues);
*/
	/*
	 * TEST ZONE
	 */
	 multigrid->preSmooth();


	/*
	 * FINALIZE PINC VARIABLES
	 */
	freeGrid(grid);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
