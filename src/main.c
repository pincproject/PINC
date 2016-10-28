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
#include "object.h"

void regular(dictionary *ini);
funPtr regular_set(dictionary *ini){ return regular; }

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE PINC
	 */
	MPI_Init(&argc,&argv);
	dictionary *ini = iniOpen(argc,argv); // No printing before this
	msg(STATUS, "PINC %s started.", VERSION);    // Needs MPI
	MPI_Barrier(MPI_COMM_WORLD);
	parseIndirectInput(ini);

	/*
	 * CHOOSE PINC RUN MODE
	 */
	void (*run)() = select(ini,"methods:mode",	regular_set,
												mgMode_set,
												mgModeErrorScaling_set,
												puModeParticle_set,
												puModeInterp_set);
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
	void (*acc)()   = select(ini,"methods:acc",	puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set);

	void (*distr)() = select(ini,"methods:distr",	puDistr3D1_set,
													puDistrND1_set,
													puDistrND0_set);

	void (*solve)() = select(ini,"methods:poisson", mgSolve_set);

	void (*extractEmigrants)() = select(ini,"methods:migrate",	puExtractEmigrants3D_set,
																puExtractEmigrantsND_set);

	// char *str;
	//
	// str = iniGetStr("methods:acc");
	// void (*acc)() = NULL;
	// if(!strcmp(str,"puAcc3D1")) acc = puAcc3D1_set();
	// if(!strcmp(str,"puAcc3D1KE")) acc = puAcc3D1KE_set();
	// if(acc==NULL) msg(ERROR,"methods:acc=%s is an invalid option")

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *deltaRho = gAlloc(ini, SCALAR);
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);
	Object *obj = oAlloc(ini);

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

	// Setting Boundary slices
	gSetBndSlices(phi, mpiInfo);

	//Set mgSolve
	MgAlgo mgAlgo = getMgAlgo(ini);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1


	/*
	 * PREPARE FILES FOR WRITING
	 */
	int rank = phi->rank;
	double *denorm = malloc((rank-1)*sizeof(*denorm));
	double *dimen = malloc((rank-1)*sizeof(*dimen));

	for(int d = 1; d < rank;d++) denorm[d-1] = 1.;
	for(int d = 1; d < rank;d++) dimen[d-1] = 1.;

	pOpenH5(ini, pop, "pop");
	gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");

    oOpenH5(ini, obj, mpiInfo, denorm, dimen, "test");

    oReadH5(obj, mpiInfo);
    
    //Compute capacitance matrix
    oComputeCapacitanceMatrix(obj, ini, mpiInfo);

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);

	/*
	 * INITIAL CONDITIONS
	 */

	//  oComputeCapacitanceMatrix(obj,phi,deltaRho);

	// Initalize particles
	pPosUniform(ini, pop, mpiInfo, rngSync);
	pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// pPosLattice(ini, pop, mpiInfo);
	// pVelZero(pop);

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
	solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
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

	Timer *t = tAlloc(rank);

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
		// oRayTrace(pop, obj, deltaRho);
		puMove(pop);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);

		// Compute charge density
		distr(pop, rho);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

		gAssertNeutralGrid(phi, mpiInfo);

		// Second run with solver to account for charges
		// oApplyCapacitanceMatrix(obj,phi,deltaRho);
		// gAdd(rho, deltaRho);
		// solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);


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
	oCloseH5(obj);
	xyCloseH5(history);

	// Free memory
	mgFree(mgRho);
	mgFree(mgPhi);
	mgFree(mgRes);
	gFree(rho);
	gFree(deltaRho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);
	oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}
