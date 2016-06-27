/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 */

#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <hdf5.h>

#include "core.h"
#include "iniparser.h"
#include "pusher.h"
#include "multigrid.h"

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv); // No stdout before this, needs MPI
	msg(STATUS|ONCE,"PINC started.");
	MPI_Barrier(MPI_COMM_WORLD);

	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);

	Population *pop = pAlloc(ini);
	Grid *rho = gAlloc(ini,1);
	MpiInfo *mpiInfo = gAllocMpi(ini);

	pPosUniform(ini,pop,mpiInfo,rngSync);
	pPosPerturb(ini,pop,mpiInfo);
	puMigrate(pop,mpiInfo,rho);

	puDistr3D1(pop,rho);	// Possible source
	gHaloOp(addSlice,rho,mpiInfo);	// Possible source

	double denorm[3];
	adSetAll(denorm,3,1.0);
	gOpenH5(ini,rho,mpiInfo,denorm,denorm,"rho");
	gWriteH5(rho,mpiInfo,0.0);
	gCloseH5(rho);

	gFreeMpi(mpiInfo);
	gFree(rho);
	pFree(pop);

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
