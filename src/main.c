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

	pOpenH5(ini, pop, "pop");
	gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gOpenH5(ini, E, mpiInfo, denorm, dimen, "E");

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);


	/***************************************************************
	 *		ACTUAL simulation stuff
	 **************************************************************/
	int nTimeSteps = iniparser_getint(ini, "time:nTimeSteps", 0);

	//Inital conditions
	pPosUniform(ini, pop, mpiInfo, rng);
	// TBD: INITIAL VELOCITY

	//Get initial E-field
	puDistr3D1(pop, rho);
	gHaloOp(setSlice, rho, mpiInfo);
	mgSolver(mgVRegular, mgRho, mgPhi, mgRes, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo);

	//Half step
	gMul(E, 0.5);
	puAcc3D1(pop, E);
	gMul(E, 2.0);

	// Time loop
	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	for(int n = 1; n < nTimeSteps; n++){

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
		puAcc3D1KE(pop, E);		// Includes kinetic energy for step n

		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		//Write h5 files
		pWriteH5(pop, mpiInfo, (double) n, (double) n);
		pWriteEnergy(history,pop,(double)n);
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
	xyCloseH5(history);

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

/**************************************************************
 *			TEMP, to reading of h5 files are ready
 *************************************************************/

void debugFillHeaviside(Grid *grid, MpiInfo *mpiInfo){

	//Load
	int *size = grid->size;
	int *trueSize = grid->trueSize;
	long int *sizeProd = grid->sizeProd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	double *val = grid->val;

	//Hard code try
	long int ind = 0;
	if(nSubdomains[2]==1){
		for(int j = 1; j < size[1]-1; j++){
			for (int k = 1; k<size[2]-1; k++) {
				for(int l = 1; l < size[3]-1; l++){
					ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
					if(l < trueSize[3]/2) val[ind] = -1;
					else if (l == trueSize[3]/2. || l == trueSize[3]) val[ind] = 0.;
					else val[ind] = 1.;
				}
			}
		}
	} else {
		for(int j = 1; j < size[1]-1; j++){
			for (int k = 1; k<size[2]-1; k++) {
				for(int l = 1; l < size[3]-1; l++){
					ind = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
					if(subdomain[2]<nSubdomains[2]/2) val[ind] = -1.;
					else val[ind] = 1.;
				}
			}
		}
		//Set in 0 at between domains (use sendSlice since it is reset every time it is used)
		double *slice = grid->sendSlice;
		for(int l = 0; l < size[3]; l++)	slice[l] = 0;
		setSlice(slice, grid, 3, size[3]-2);
	}

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

	//Prep to store grids
	int rank = rho->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gOpenH5(ini, res, mpiInfo, denorm, dimen, "res");

	free(denorm);
	free(dimen);

	// Timer *t = tAlloc(rank);

	double tol = 50;
	double err = 10001;

	//Compute stuff
	debugFillHeaviside(rho, mpiInfo);

	Timer *t = tAlloc();

	fMsg(ini, "mgLog", "New run \n\n");

	while(err>tol){
		//Run solver
		tStart(t);
		mgSolver(mgVRegular, mgRho, mgPhi, mgRes, mpiInfo);
		tStop(t);
		//Compute residual and mass
		mgResidual(res,rho, phi, mpiInfo);
		err = mgResMass3D(res,mpiInfo);
		// msg(STATUS|ONCE, "The error mass (e^2) is %f, time %llu", err, t->total);
		tMsg(t->total, "Hello:");
		fMsg(ini, "mgLog", "Hello \n");

	}

	gWriteH5(rho,mpiInfo,0.);
	gWriteH5(phi,mpiInfo,0.);
	gWriteH5(res,mpiInfo,0.);


	gCloseH5(phi);
	gCloseH5(rho);
	gCloseH5(res);

	tFree(t);

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
