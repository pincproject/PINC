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

#include "core.h"
#include "pusher.h"
#include "multigrid.h"

/*
 * MAIN ROUTINES (RUN MODES PERHAPS A BETTER NAME?)
 */
void regularRoutine(dictionary *ini);
void debugFillHeaviside(Grid *grid, MpiInfo *mpiInfo);
void mgRoutine(dictionary *ini);

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE PINC
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv); // No printing before this
	msg(STATUS|ONCE, "PINC started.");    // Needs MPI
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * CHOOSE PINC MAIN ROUTINE (RUN MODE PERHAPS A BETTER NAME?)
	 */
	char *routine = iniGetStr(ini,"main:routine");
	if(!strcmp(routine, "regular"))		regularRoutine(ini);
	if(!strcmp(routine, "mgRoutine"))	mgRoutine(ini);
	free(routine);

	/*
	 * FINALIZE PINC
	 */
	iniClose(ini);
	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}

void regularRoutine(dictionary *ini){

	// Random number seeds
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
//	gsl_rng_set(rng,2);

	/*
	 * INITIALIZE PINC VARIABLES
	 */


	// MPI struct
	MpiInfo *mpiInfo = gAllocMpi(ini);

	 // Setting up particles.
	Population *pop = pAlloc(ini);

	// Allocating grids
	Grid *E   = gAlloc(ini, 3);
	Grid *rho = gAlloc(ini, 1);
	Grid *res = gAlloc(ini, 1);
	Grid *phi = gAlloc(ini, 1);

	// Creating a neighbourhood in the rho Grid variable to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	// Alloc multigrids
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);

	//Set mgSolver
	MgAlgo mgAlgo = getMgAlgo(ini);



	// Alloc h5 files
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	pOpenH5(ini, pop, "pop");
	gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);

	/***************************************************************
	 *		ACTUAL simulation stuff
	 **************************************************************/
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");


	// Initalize particles
//	pPosUniform(ini, pop, mpiInfo, rng);
	pPosLattice(ini, pop, mpiInfo);
	pVelZero(pop);

	// Perturb particles
	pPosPerturb(ini, pop, mpiInfo);

	puExtractEmigrants3D(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	// Get initial charge density
	puDistr3D1(pop, rho);
	gHaloOp(addSlice, rho, mpiInfo, 1);

	// Get initial E-field
	mgSolver(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, 0);	// ?? What does this do?

	// Advance velocities half a step
	gMul(E, 0.5);
	puAcc3D1KE(pop, E);		// Includes kinetic energy for step n
	gMul(E, 2.0);

	// Time loop
	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS|ONCE,"Computing time-step %i",n);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,1.0);

		// Move particles
		puMove(pop);

		// Migrate particles (periodic boundaries)
		puExtractEmigrants3D(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		puDistr3D1(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo, 1);

		// Compute electric potential phi
		gZero(phi);
		gZero(res);
		mgSolver(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		// gMul(E, normE);
		gHaloOp(setSlice, E, mpiInfo, 0);

		// Apply external E
		// gAddTo(Ext);

		// Accelerate particle and compute kinetic energy for step n
		puAcc3D1KE(pop, E);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		//Write h5 files
		gWriteH5(E, mpiInfo, (double) n);
		gWriteH5(rho, mpiInfo, (double) n);
		gWriteH5(phi, mpiInfo, (double) n);
		// pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		pWriteEnergy(history,pop,(double)n);

	}


	 /*
	 * FINALIZE PINC VARIABLES
	 */
	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(phi);
	gCloseH5(E);
	xyCloseH5(history);

	// Free memory
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
	gsl_rng_free(rng);

}



void mgRoutine(dictionary *ini){

	//Mpi
	MpiInfo *mpiInfo = gAllocMpi(ini);

	//Rand Seed
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	//Grids
	Grid *phi = gAlloc(ini, 1);
	Grid *rho = gAlloc(ini, 1);
	Grid *res= gAlloc(ini, 1);
	Grid *sol = gAlloc(ini, 1);
	Grid *E = gAlloc(ini, 3);
	Grid *error =gAlloc(ini, 1);

	//Multilevel grids
	Multigrid *mgPhi = mgAlloc(ini, phi);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);

	MgAlgo mgAlgo = getMgAlgo(ini);

	//Sets the boudary slices
	gSetBndSlices(phi, mpiInfo);
	mgRestrictBnd(mgPhi);

	// gBnd(phi, mpiInfo);
	// gBnd(mgPhi->grids[1], mpiInfo);

	// dumpWholeGrid(ini, phi);
	// dumpWholeGrid(ini, mgPhi->grids[1]);
	//
	// return;


	int rank = rho->rank;
	Timer *t = tAlloc(rank);

	// double tol = 77000.;
	// double err = tol+1.;

	//Compute stuff
	fillHeaviside(rho, mpiInfo);
	fillHeaviSol(sol, mpiInfo);
	// fillPointCharge(rho, mpiInfo);
	// fillPolynomial(rho, mpiInfo);
	// fillPointSol(sol, mpiInfo);
	// fillExp(sol, mpiInfo);
	// fillSin(rho, mpiInfo);
	// fillSinSol(sol, mpiInfo);
	// fillCst(rho, mpiInfo);
	// fillRng(rho, mpiInfo, rng);

	// gHaloOp(setSlice, sol, mpiInfo);
	// gFinDiff2nd3D(rho, sol);

	msg(STATUS|ONCE, "mgLevels = %d", mgRho->nLevels);
	gNeutralizeGrid(rho, mpiInfo);

	double tol = 10000;
	double errSquared = 10001;
	double resSquared = 10001;
	double counter = tol+1;

	while(counter>tol){
		counter--;
		// Run solver
		gZero(res);
		tStart(t);
		mgSolver(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		// for(int n = 0; n < mgRho->nMGCycles; n++){
		// // // 	// mgGS3D(phi, rho, mgRho->nPreSmooth, mpiInfo);
		// // // 	// mgGS3D(phi, rho, mgRho->nPostSmooth, mpiInfo);
		// // // 	// mgJacob3D(phi, rho, mgRho->nPreSmooth, mpiInfo);
		// // // 	// mgJacob3D(phi, rho, mgRho->nPostSmooth, mpiInfo);
		// 	mgRho->preSmooth(phi, rho, mgRho->nPreSmooth, mpiInfo);
		// 	mgRho->postSmooth(phi, rho, mgRho->nPreSmooth, mpiInfo);
		// // //
		// }

		tStop(t);
		//
		// //Compute error
		mgCompError(phi,sol,error);
		errSquared = mgSumTrueSquared(error, mpiInfo);
		resSquared = mgSumTrueSquared(res, mpiInfo);


		// Compute residual and mass
		// gZero(res);
		// gHaloOp(setSlice, rho, mpiInfo, 0);
		// gHaloOp(setSlice, phi,mpiInfo, 0);
		// gBnd(phi, mpiInfo);
		// mgResidual(res,rho, phi, mpiInfo);
		// gHaloOp(setSlice, res, mpiInfo, 0);
		// err = mgResMass3D(res,mpiInfo);
		msg(STATUS|ONCE, "Error squared (e^2) = %f", errSquared);
		msg(STATUS|ONCE, "Residual squared (res^2) = %f", resSquared);
			// The res squared (res^2) = %f", errSquared, resSquared);
	}


	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");


	/*********************************************************************
	 *			STORE GRIDS
	 ********************************************************************/

	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	//(Re)Compute E, error and residual
	gHaloOp(setSlice, phi, mpiInfo, 0);
	gNeutralizeGrid(phi, mpiInfo);
	gBnd(phi, mpiInfo);
	gFinDiff1st(phi, E);
	mgCompError(phi,sol,error);
	mgResidual(res,rho, phi, mpiInfo);


	gOpenH5(ini, E, mpiInfo, denorm, dimen, "E_0");
	gWriteH5(E, mpiInfo, 0.);
	gCloseH5(E);

	gOpenH5(ini, sol, mpiInfo, denorm, dimen, "sol_0");
	gWriteH5(sol, mpiInfo, 0.);
	gCloseH5(sol);

	gOpenH5(ini, error, mpiInfo, denorm, dimen, "error_0");
	gWriteH5(error, mpiInfo, 0.);
	gCloseH5(error);


	//Saving lvl of grids
	char fName[12];
	for(int lvl = 0; lvl <mgRho->nLevels; lvl ++){

		rho = mgRho->grids[lvl];
		phi = mgPhi->grids[lvl];
		res = mgRes->grids[lvl];

		sprintf(fName, "rho_%d", lvl);
		gOpenH5(ini, rho, mpiInfo, denorm, dimen, fName);
		sprintf(fName, "phi_%d", lvl);
		gOpenH5(ini,  phi, mpiInfo, denorm, dimen, fName);
		sprintf(fName, "res_%d", lvl);
		gOpenH5(ini, res, mpiInfo, denorm, dimen, fName);


		gWriteH5(rho,mpiInfo,0.);
		gWriteH5(phi,mpiInfo,0.);
		gWriteH5(res,mpiInfo,0.);

		gCloseH5(phi);
		gCloseH5(rho);
		gCloseH5(res);

	}


	free(denorm);
	free(dimen);


	tFree(t);
	gFreeMpi(mpiInfo);

	gsl_rng_free(rng);

	return;
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
