/**
 * @file	    main.c
 * @brief	    PINC main routine.
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 *
 * Main routine for PINC (Particle-IN-Cell).
 */

#include "core.h"
#include "pusher.h"
#include "multigrid.h"
#include "spectral.h"
#include "object.h"

/**
 * @brief Regular PIC simulation
 */
void regular(dictionary *ini);
funPtr regular_set(dictionary *ini){ return regular; }

/**
 * @brief Do a pure Poisson solve
 */
void poisson(dictionary *ini);
funPtr poisson_set(dictionary *ini){ return poisson; }

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE PINC
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv); // No printing before this
	msg(STATUS, "PINC %s started.", VERSION);    // Needs MPI
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * CHOOSE PINC RUN MODE
	 */
	void (*run)() = select(ini,"methods:mode",	regular_set,
                                                poisson_set,
												mgMode_set,
												mgModeErrorScaling_set,
												sMode_set);
	run(ini);

	/*
	 * FINALIZE PINC
	 */
	iniClose(ini);
	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}

void regular(dictionary *ini){

	/*
	 * SELECT METHODS
	 */
	void (*acc)()   			= select(ini,	"methods:acc",
												puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set);

	void (*distr)() 			= select(ini,	"methods:distr",
												puDistr3D1_set,
												puDistrND1_set,
												puDistrND0_set);

	void (*extractEmigrants)()	= select(ini,	"methods:migrate",
												puExtractEmigrants3D_set,
												puExtractEmigrantsND_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set,
												sSolver_set);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	void *solver = solverAlloc(ini, rho, phi);
	// Object *obj = oAlloc(ini);

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */
	pOpenH5(ini, pop, units, "pop");
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");
  // oOpenH5(ini, obj, mpiInfo, units, 1, "test");
  // oReadH5(obj, mpiInfo);

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	 * INITIAL CONDITIONS
	 */

	// Initalize particles
	// pPosUniform(ini, pop, mpiInfo, rngSync);
	pPosLattice(ini, pop, mpiInfo);
	pVelZero(pop);
	// pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// Perturb particles
	pPosPerturb(ini, pop, mpiInfo);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);

	puMigrate(pop, mpiInfo, rho);

	/*
	 * INITIALIZATION (E.g. half-step)
	 */

	// Get initial charge density
	distr(pop, rho);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

	// Get initial E-field
	solve(solver, rho, phi, mpiInfo);
	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);


	// Advance velocities half a step
	gMul(E, 0.5);
	acc(pop, E);
	gMul(E, 2.0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i",n);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);

		// Move particles
		puMove(pop);
		// oRayTrace(pop, obj);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		distr(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		// gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		// solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

		// mgSolve(solver, rho, phi, mpiInfo);
		// sSolve(solver, rho, phi, mpiInfo);

		solve(solver, rho, phi, mpiInfo);

		gHaloOp(setSlice, phi, mpiInfo, TOHALO); // Needed by sSolve but not mgSolve

		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);

		gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);

		// Accelerate particle and compute kinetic energy for step n
		acc(pop, E);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		//Write h5 files
		// gWriteH5(E, mpiInfo, (double) n);
		// gWriteH5(rho, mpiInfo, (double) n);
		gWriteH5(phi, mpiInfo, (double) n);
		// pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		pWriteEnergy(history,pop,(double)n);

	}

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */
	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(phi);
	gCloseH5(E);
	// oCloseH5(obj);
	xyCloseH5(history);

	// Free memory
	// sFree(solver);
	// mgFreeSolver(solver);
	solverFree(solver);
	gFree(rho);
	gFree(phi);
	gFree(E);
	pFree(pop);
	uFree(units);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}

void poisson(dictionary *ini){

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set,
												sSolver_set);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	void *solver = solverAlloc(ini, rho, phi);

	// Setting Boundary slices 
	gSetBndSlices(phi, mpiInfo);

	/*
	 * PREPARE FILES FOR WRITING
	 */
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");

    /* gReadH5(rho, mpiInfo, 0.); */
    gFillSin(rho, 1, mpiInfo, 0);
    gWriteH5(rho, mpiInfo, 0.);
    /* adPrint(rho->val, rho->sizeProd[rho->rank]); */
    adPrint(rho->val, 10);

	Timer *t = tAlloc(mpiInfo->mpiRank);

    tStart(t);

    gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
    solve(solver, rho, phi, mpiInfo);
    gHaloOp(setSlice, phi, mpiInfo, TOHALO); // Needed by sSolve but not mgSolve

    tStop(t);

    adPrint(phi->val, 10);
    gWriteH5(phi, mpiInfo, 0.);

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */
	gFreeMpi(mpiInfo);

	// Close h5 files
	gCloseH5(rho);
	gCloseH5(phi);

	// Free memory
	solverFree(solver);
	gFree(rho);
	gFree(phi);
	uFree(units);
}
