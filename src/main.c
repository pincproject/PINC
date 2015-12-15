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
	Grid *grid = allocGrid(ini);
	GridQuantity *rho = allocGridQuantity(ini, grid, 1);
	GridQuantity *phi = allocGridQuantity(ini, grid, 1);
	Multigrid *multiRho = allocMultigrid(ini, rho);
	Multigrid *multiPhi = allocMultigrid(ini, phi);
	MpiInfo *mpiInfo = allocMpiInfo(ini);

	/*
	 *	TEST AREA
	 */


 // 	testGetSlice(ini, gridQuantity);
 // 	testSwapHalo(ini, gridQuantity, mpiInfo);
 // 	testGaussSeidel(ini, multiRho, multiPhi, mpiInfo);
	// testRestriction(ini, multiRho, multiPhi, mpiInfo);
	// testMultigrid(ini, multiRho, multiPhi,mpiInfo);
	testDerivatives(ini, rho, phi);

	/*
	 * FINALIZE PINC VARIABLES
	 */
	freeGrid(grid);
	freeGridQuantity(rho);
	freeGridQuantity(phi);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
