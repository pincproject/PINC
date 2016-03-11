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

void regularRoutine(dictionary *ini){

	//Random number seeds
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	/*
	 * INITIALIZE PINC VARIABLES
	 */

	//MPI struct
	MpiInfo *mpiInfo = gAllocMpi(ini);

	 //Setting up particles.
	Population *pop = pAlloc(ini);

	//Allocating grids
	Grid *E = gAlloc(ini, 3);
	Grid *rho = gAlloc(ini, 1);
	Grid *res = gAlloc(ini, 1);
	Grid *phi = gAlloc(ini, 1);

	//Alloc multigrids
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);

	//Alloc h5 files
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	pCreateH5(ini, pop, "pop");
	gCreateH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gCreateH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gCreateH5(ini, E, mpiInfo, denorm, dimen, "E");

	free(denorm);
	free(dimen);


	/***************************************************************
	 *		ACTUAL simulation stuff
	 **************************************************************/
	int nTimesteps = iniparser_getint(ini, "time:nTimesteps", 0);

	//Inital conditions
	pPosUniform(ini, pop, mpiInfo, rng);



	//Get initial E-field
	puDistr3D1(pop, rho);
	gHaloOp(setSlice, rho, mpiInfo);
	mgSolver(mgVRegular, mgRho, mgPhi, mgRes, mpiInfo);
	gFinDiff1st(phi, E);

	//Half step
	gMul(E, 0.5);
	puAcc3D1(pop, E);
	gMul(E, 2.0);

	 //Time loop
	for(int t = 0; t < nTimesteps; t++){

		//Move particles
		puMove(pop);
		puMigrate(pop, mpiInfo, E);

		//Compute E field
		puDistr3D1(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo);
		mgSolver(mgVRegular, mgRho, mgPhi, mgRes, mpiInfo);
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo);

		//Apply external E
		// gAddTo(Ext);

		//Accelerate
		puAcc3D1(pop, E);		// Including total kinetic energy for step n

		//Write h5 files
		pWriteH5(pop, mpiInfo, (double) t, (double) t);
	}


	 /*
	 * FINALIZE PINC VARIABLES
	 */
	 gFreeMpi(mpiInfo);

	//  //Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(phi);
	gCloseH5(E);

	//Free memory
	mgFree(mgRho);
	mgFree(mgPhi);
	mgFree(mgRes);
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);


	 /*
		* FINALIZE THIRD PARTY LIBRARIES
		*/
	// iniClose(ini); 	//No iniClose??
	gsl_rng_free(rng);

	return;
}

void mgRoutine(dictionary *ini){

	//Mpi
	MpiInfo *mpiInfo = gAllocMpi(ini);

	//Grids
	Grid *phi = gAlloc(ini, 1);
	Grid *rho = gAlloc(ini, 1);
	Grid *res= gAlloc(ini, 1);

	//Multilevel grids
	Multigrid *mgPhi = mgAlloc(ini, phi);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);

	//Fill rho with lazy heaviside
	int *subdomain = mpiInfo->subdomain;
	int rank = rho->rank;
	long int *sizeProd = rho->sizeProd;
	double *val = rho->val;

	double charge = (double) (subdomain[0]>0); // -1 + 2*(subdomain[0]>0);
	for(int g = 0; g < sizeProd[rank]; g++) val[g] = charge;

	//Prep to store grids
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	gCreateH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gCreateH5(ini, phi, mpiInfo, denorm, dimen, "phi");

	free(denorm);
	free(dimen);

	//Run solver
	mgSolver(mgVRegular, mgRho, mgPhi, mgRes, mpiInfo);

	//Compute residual and mass
	mgResidual(res,rho, phi, mpiInfo);
	double mass = mgResMass3D(res,mpiInfo);
	if(mpiInfo->mpiRank == 0)	msg(STATUS, "The residual mass is %f ",mass);

	gWriteH5(rho,mpiInfo,0.);
	gWriteH5(phi,mpiInfo,0.);
	gWriteH5(res,mpiInfo,0.);

	gCloseH5(phi);
	gCloseH5(rho);
	gCloseH5(res);






	return;
}


int main(int argc, char *argv[]){


	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv);
	msg(STATUS|ONCE,"PINC started.");
	MPI_Barrier(MPI_COMM_WORLD);

	//Choose routine from ini file
	char *routine = iniparser_getstring(ini, "main:routine", "\0");

	if(!strcmp(routine, "regular"))				regularRoutine(ini);
	if(!strcmp(routine, "mgRoutine"))			mgRoutine(ini);


	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}

/*****************************************************************
 *			Blueprint
 ****************************************************************/

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
