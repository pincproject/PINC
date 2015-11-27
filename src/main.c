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
#include <hdf5.h>

#include "pinc.h"
#include "iniparser.h"
#include "multigrid.h"
#include "test.h"

int main(int argc, char *argv[]){


	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);

	msg(STATUS|ONCE,"PINC started.");
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * INITIALIZE PINC VARIABLES
	 */

	dictionary *ini = iniOpen(argc,argv);

	Population *pop = allocPopulation(ini);
	Grid *grid = allocGrid(ini);
	GridQuantity *gridQuantity = allocGridQuantity(ini, grid, 1);
	Multigrid *multigrid = allocMultigrid(ini, gridQuantity);
	MpiInfo *mpiInfo = allocMpiInfo(ini);

	/*
	 * 	Test Area
	 */

	/*
	 * FINALIZE PINC VARIABLES
	 */

	// testGetSlice(ini, gridQuantity);
	testSwapHalo(ini, gridQuantity, mpiInfo);


	freeMpiInfo(mpiInfo);
	freePopulation(pop);
	freeGrid(grid);
	freeGridQuantity(gridQuantity);
//	freeMultigrid(multigrid);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
