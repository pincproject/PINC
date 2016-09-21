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

void regular(dictionary *ini);
funPtr regular_set(dictionary *ini){ return regular; }

void mgRun(dictionary *ini);
funPtr mgRun_set(dictionary *ini){ return mgRun; }

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
	void (*run)() = select(ini,"methods:mode",regular_set,mgRun_set);
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

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *Et  = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);

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

	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	free(denorm);
	free(dimen);

	/*
	 * INITIAL CONDITIONS
	 */

	double left = E->nGhostLayers[1];
	int *L = gGetGlobalSize(ini);
	double midway = L[0]/2.0;


	double stepSize = E->stepSize[1];

	double pos[] = {0.5*midway+left};
	double vel[] = {0.0};
	pNew(pop,0,pos,vel);

	double slope =  1.0/midway;
	msg(STATUS,"midway: %f, L: %i",midway,L[0]);

	free(L);

	for(int g=1; g<Et->trueSize[1]+1; g++){
		phi->val[g] = -0.5*stepSize*slope*pow((double)g-midway-left,2);
		Et->val[g] = slope*(g-midway-left);
	}
	adPrint(Et->val, Et->size[1]);
	gNormalizeE(ini, Et);
	adPrint(Et->val, Et->size[1]);
	gNormalizePhi(ini, phi);
	// adPrint(phi->val, phi->size[1]);
	gHaloOp(setSlice,phi,mpiInfo,TOHALO);
	// adPrint(phi->val, phi->size[1]);
	gFinDiff1st(phi, E);
	gMul(E,-1.0);
	// adPrint(E->val, E->size[1]);
	gHaloOp(setSlice,E,mpiInfo,TOHALO);
	adPrint(E->val, E->size[1]);

	pWriteH5(pop, mpiInfo, 0.0, 0.0);

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

		// msg(STATUS,"Computing time-step %i",n);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,1.0);

		tStart(t);

		// Move particles
		// msg(STATUS,"n: %i, pos: %.2f",n,pop->pos[0]-1);
		puMove(pop);
		puPeriodic(pop,E);

		pWriteH5(pop, mpiInfo, (double) n, (double)n-0.5);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);

		acc(pop, E);

		tStop(t);

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
	xyCloseH5(history);

	// Free memory
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	gFree(Et);
	pFree(pop);

	gsl_rng_free(rngSync);

}
