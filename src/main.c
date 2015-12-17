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
#include "pusher.h"
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
	Grid *rho = gAlloc(ini, 1);
	Grid *phi = gAlloc(ini, 1);


	Multigrid *multiRho = mgAlloc(ini, rho);
	// Multigrid *multiPhi = mgAlloc(ini, phi);

	MpiInfo *mpiInfo = gAllocMpi(ini);

	//
	// int nDims = grid->nDims;
	// double *denorm = malloc(sizeof(*denorm));
	// double *dimens = malloc(sizeof(*dimens));
	//
	// for(int d = 0; d < nDims;d++) denorm[d] = 1.;
	// for(int d = 0; d < nDims;d++) dimens[d] = 1.;
	//
	// createGridQuantityH5(ini,rho,mpiInfo, denorm, dimens, "rho");
	// createGridQuantityH5(ini,phi,mpiInfo,denorm,dimens, "phi");

	/*
	 *	TEST AREA
	 */


 // 	testGetSlice(ini, gridQuantity);
 // 	testSwapHalo(ini, gridQuantity, mpiInfo);
 // 	testGaussSeidel(ini, multiRho, multiPhi, mpiInfo);
 // 	testRestriction(ini, multiRho, multiPhi, mpiInfo);
 // 	testMultigrid(ini, multiRho, multiPhi,mpiInfo);
 // 	testDerivatives(ini, rho, phi);

	/*
	 * FINALIZE PINC VARIABLES
	 */
	gFree(rho);
	gFree(phi);
	gFreeMpi(mpiInfo);
	mgFree(multiRho);

	//
	// closeGridQuantityH5(rho);
	// closeGridQuantityH5(phi);
	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
