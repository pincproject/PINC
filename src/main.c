/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
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

	// -2. Sanitize

	// -1. Allocate all datatypes

	// 0. Specify all function pointers from ini

	// 1. Specify phase space distribution
	// 2. rho: puDistr3D1(); (distribute)
	// 3. phi: linearMG();
	// 4. E: gFinDiff1rd();

	// Accelerate half-step
	// gMul(E,0.5);
	// puAcc3D1(pop,E); // (accelerate)
	// gMul(E,2);

	/*
	 *	TEST AREA
	 */

	// for(long int n=1;n<=N;n++){
	//
	// 	// Everything in here is function pointers
	//
	// 	// Move
	// 	move();
	// 	migrate();			// Including boundaries and safety testing
	//
	// 	// Weighting
	// 	distribute();
	// 	interactAdd();		// Fix boundaries for rho
	//
	// 	// Field solver
	// 	solver();			// Including boundaries
	// 	finDiff();
	// 	swapHalo();
	//
	// 	imposeExternal();	// To add external field
	// 	potentialEnergy();	// Calculate potential energy for step n
	//
	// 	// Accelerate
	// 	accelerate();		// Including total kinetic energy for step n
	//
	// 	// Diagnostics
	// 	if(n%a==0) savePop();
	// 	if(n%b==0) saveGrid();
	// 	if(n%c==0) saveVelocityDistr();
	//
	// }

	/*
	 * FINALIZE PINC VARIABLES
	 */

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
