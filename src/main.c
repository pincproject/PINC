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
#include "collisions.h"

void regular(dictionary *ini);
funPtr regular_set(dictionary *ini){ return regular; }

//delete BorisTest in final Version
void BorisTestMode(dictionary *ini);
funPtr BorisTestMode_set(dictionary *ini){ return BorisTestMode; }


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
												mccTestMode_set,
												BorisTestMode_set,
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
		// gWriteH5(phi, mpiInfo, (double) n);
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


// delete below here in final version

void BorisTestMode(dictionary *ini){


	/*
	 * SELECT METHODS
	 */
	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set,
												puBoris3D1_set,
												puBoris3D1KE_set,
												puBoris3D1KETEST_set);

	void (*distr)() = select(ini,"methods:distr",	puDistr3D1_set,
													puDistrND1_set,
													puDistrND0_set);

	void (*solverInterface)() = select(ini,"methods:poisson", mgSolver_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
																puExtractEmigrantsND_set);

	// char *str;
	//
	// str = iniGetStr("methods:acc");
	// void (*acc)() = NULL;
	// if(!strcmp(str,"puAcc3D1")) acc = puAcc3D1_set();
	// if(!strcmp(str,"puAcc3D1KE")) acc = puAcc3D1KE_set();
	// if(acc==NULL) msg(ERROR,"methods:acc=%s is an invalid option")



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
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	void *solver = solverAlloc(ini, rho, phi);



	// Object *obj = oAlloc(ini);

	// using Boris algo
	int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	//Set mgSolve
	//MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1


	/*
	 * PREPARE FILES FOR WRITING
	 */

	//---------------------------REMOVE v ?
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;
	//----------------------------------------------------

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

	free(denorm);
	free(dimen);

	/*
	 * INITIAL CONDITIONS
	 */

	// Initalize particles
	//pPosUniform(ini, pop, mpiInfo, rngSync);
	//pPosLattice(ini, pop, mpiInfo);
	//pVelZero(pop);
	//double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	//msg( STATUS, "velthermal1 = %f, velthermal2 = %f", velThermal[0], velThermal[1]);
	//pVelConstant(ini, pop, velThermal[0], velThermal[1]); //constant values for vel.
	//pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// Manually initialize a single particle
	if(mpiInfo->mpiRank==0){
		double pos[3] = {8., 8., 8.};
		double vel[3] = {0.02, 0., 1.};
		pNew(pop, 0, pos, vel);
		//double pos1[3] = {17., 17., 16.};
		//double vel1[3] = {0.1, 0., 0.1};
		//pNew(pop, 1, pos1, vel1); //second particle
	}

	// Perturb particles
	//pPosPerturb(ini, pop, mpiInfo);

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
	puAddEext(ini, pop, E);

	gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	// adScale(T, 3*nSpecies, 0.5);
	// adScale(S, 3*nSpecies, 0.5);

	gZero(E);
	acc(pop, E, T, S);

	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	double x_min = 20;
	double x_max = 0;
	double y_min = 20;
	double y_max = 0;

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(rank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i ",n);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);

		// Move particles
		adPrint(pop->pos, 3);
		x_min = pop->pos[0]<x_min ? pop->pos[0] : x_min;
		x_max = pop->pos[0]>x_max ? pop->pos[0] : x_max;
		y_min = pop->pos[1]<y_min ? pop->pos[1] : y_min;
		y_max = pop->pos[1]>y_max ? pop->pos[1] : y_max;
		puMove(pop);
		// oRayTrace(pop, obj);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);
		//msg(STATUS, "HEEEEEEEEEEERRRRRRRRRRREEEEEEEEEEEE");
		// Compute charge density
		distr(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		solve(solver, rho, phi, mpiInfo);

		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);

		gAssertNeutralGrid(E, mpiInfo);


		// Apply external E
		// gAddTo(Ext);
		gZero(E); //remove interparticle forces
		puAddEext(ini, pop, E);

		//puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n

		acc(pop, E, T, S);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		//Write h5 files
		//gWriteH5(E, mpiInfo, (double) n);
		gWriteH5(rho, mpiInfo, (double) n);
		gWriteH5(phi, mpiInfo, (double) n);
		pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		pWriteEnergy(history,pop,(double)n);

	}

	msg(STATUS, "x in [%f, %f]", x_min, x_max);
	msg(STATUS, "y in [%f, %f]", y_min, y_max);

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	 * FINALIZE PINC VARIABLES
	 */
	free(S);
	free(T);

	gFreeMpi(mpiInfo);

	// Close h5 files
	pCloseH5(pop);
	//gCloseH5(rho);
	gCloseH5(phi);
	//gCloseH5(E);
	// oCloseH5(obj);
	xyCloseH5(history);

	// Free memory
	//mgFree(mgRho);
	//mgFree(mgPhi);
	//mgFree(mgRes);
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);
	solverFree(solver);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}
