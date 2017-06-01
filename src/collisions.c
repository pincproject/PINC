/**
* @file		multigrid.c
* @brief		Collisional module, Montecarlo Collision Method.
* @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
*
* Main module collecting functions conserning collisions
*
* MCC details......
*
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_randist.h>

#include "core.h"
#include "collisions.h"
#include "multigrid.h"
#include "pusher.h"

/******************************************************************************
* 				Local functions
*****************************************************************************/

int mccTest(int one, int two){
	//bogus test to check if function calls and includes work
	return one+two;
}

/******************************************************************************
* 				Global functions
*****************************************************************************/


/******************************************************************************
* 				Local functions
*****************************************************************************/

static void mccSanity(dictionary *ini, const char* name, int nSpecies){

	int countSpecies = iniGetInt(ini,"population:nSpecies");
	if(countSpecies!=nSpecies){
		msg(ERROR,"%s only supports 2 species",name);
	}
}

void mccReadcrossect(const char* filename){
	FILE* fp = fopen(filename, "r" );

	fclose(fp);
}

void mccCollideElectron(Population *pop,  double Pnull, const gsl_rng *rng){


	msg(STATUS,"colliding Electrons");
	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;
	//double *vel = &pop.vel[i*pop.nDims];
	double R = gsl_rng_uniform(rng);

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	double vx, vy, vz;

	//for(int s=0; s<nSpecies; s++){
	int errorcounter = 0;

	double Ekin = 0;
	long int pStart = pop->iStart[0];		// *nDims??
	long int pStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (1-Pnull)*(pop->iStop[0]-pop->iStart[0])/nDims; // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((pStop-pStart)*nDims)/ (double)(NparticleColl));
	if ((mccStepSize-1)*(NparticleColl)/nDims > (pStop-pStart)){ //errorcheck, remove
		msg(WARNING,"particle collisions out of bounds in mccCollideElectron");
	}
	msg(STATUS,"number of particles in array = %i", (pStop-pStart)/nDims);
	msg(STATUS,"number of particles colliding = %i", (NparticleColl));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize/nDims));

	long int mccStop = pStart + mccStepSize*NparticleColl; //enshure non out of bounds
	for(long int p=pStart;p<mccStop;p+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pnull*(pStop-pStart);p++)
		errorcounter += 1;
		//msg(STATUS,"errorcounter = %i", (errorcounter));
		if (errorcounter > NparticleColl){ //errorcheck, remove
			msg(WARNING,"errorcounter = %i should be the same as number of coll\
			particles = %i", (errorcounter),NparticleColl/nDims);
			msg(WARNING,"particle collisions out of bounds in mccCollideElectron");
			//errorcounter = 0;
		}
		R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
		q = p + (R*(((double)mccStepSize/(double)nDims)));
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );
		vx = vel[q];
		vy = vel[q+1];
		vz = vel[q+2];
		double testy;


		// need n,n+1,n+2 for x,y,j ?
		Ekin = 0.5*(vx*vx + vy*vy + vz*vz)*mass[0]; // values stored in population for speed??
		R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
		angleChi = acos( (2+Ekin-2*pow((1+Ekin),R) )/(Ekin));
		anglePhi = 2*PI*R; // set different R?
		angleTheta = acos(vx);
		A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
		vel[q] = vx*cos(angleChi)+A*(vy*vy + vz*vz); //Vx
		vel[q+1] = vy*cos(angleChi)+A*vz-A*vx*vy; //Vy
		vel[q+2] = vz*cos(angleChi)-A*vy-A*vx*vz; //Vz

	}
	msg(STATUS,"errorcounter = %i should be the same as number of coll particles\
	= %i", (errorcounter),NparticleColl);

	//free?
}
 //1590430

void mccCollideIon(const dictionary *ini, Population *pop, double Pnull, const gsl_rng *rng){

	msg(STATUS,"colliding Ions");

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *mass = pop->mass;
	//double *vel = &pop.vel[i*pop.nDims];
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	double timeStep = iniGetDouble(ini,"time:timeStep");
	double stepSize = iniGetDouble(ini,"grid:stepSize");
	double velTh = (timeStep/stepSize)*velThermal[1]; //Ions !?

	double R = gsl_rng_uniform(rng);
	double fractR = gsl_rng_uniform(rng);

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	double fraction = 0.7; //defines fraction of particles colliding w el/ch-ex

	double vxMW, vyMW, vzMW;

	//for(int s=0; s<nSpecies; s++){
	int errorcounter = 0;

	double Ekin = 0;

	long int pStart = pop->iStart[1];			// *nDims??
	long int pStop = pop->iStop[1];      // make shure theese are actually ions!

	long int NparticleColl = (1-Pnull)*(pop->iStop[1]-pop->iStart[1])/nDims; // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((pStop-pStart)*nDims)/ (double)(NparticleColl));
	if ((mccStepSize-1)*(NparticleColl)/nDims > (pStop-pStart)){ //errorcheck, remove
		msg(WARNING,"particle collisions out of bounds in mccCollideIon");
		msg(WARNING,"%i is bigger than array = %i",(mccStepSize-1)*(NparticleColl),(pStop-pStart) );
	}
	msg(STATUS,"number of particles in array = %i", (pStop-pStart)/nDims);
	msg(STATUS,"number of particles colliding = %i", (NparticleColl));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize/nDims));

	long int mccStop = pStart + mccStepSize*NparticleColl; //enshure non out of bounds
	for(long int p=pStart;p<mccStop;p+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pnull*(pStop-pStart);p++)
		errorcounter += 1;
		if (errorcounter > NparticleColl){ //errorcheck, remove
			msg(WARNING,"errorcounter = %i should be the same as number of coll\
			particles = %i", (errorcounter),NparticleColl/nDims);
			msg(WARNING,"particle collisions out of bounds in mccCollideIon");
			//errorcounter = 0;
		}
		fractR = gsl_rng_uniform(rng); //decides type of coll.
		R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
		q = (R*(((double)mccStepSize/(double)nDims))) + p;
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );

		vxMW = gsl_ran_gaussian_ziggurat(rng,velTh); //maxwellian dist?
		vyMW = gsl_ran_gaussian_ziggurat(rng,velTh);
		vzMW = gsl_ran_gaussian_ziggurat(rng,velTh);

		double vxTran = vel[q]-vxMW;   //simple transfer, should be picked from maxwellian
		double vyTran = vel[q+1]-vyMW;  //vel-velMaxwellian
		double vzTran = vel[q+2]-vzMW;  // this will for now break conservation laws
		if(fractR < fraction){ // if test slow, but enshures randomness...
			// elastic:
			//msg(STATUS,"elastic");
			// need n,n+1,n+2 for x,y,j ?
			Ekin = 0.5*(vxTran*vxTran + vyTran*vyTran + vzTran*vzTran)*mass[1]; // values stored in population for speed??
			R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
			angleChi = acos( (2+Ekin-2*pow((1+Ekin),R) )/(Ekin));
			anglePhi = 2*PI*R; // set different R?
			angleTheta = acos(vxTran);
			A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
			vel[q] = vxTran*cos(angleChi)+A*(vyTran*vyTran + vzTran*vzTran) + vxMW; //Vx
			vel[q+1] = vyTran*cos(angleChi)+A*vzTran-A*vxTran*vyTran + vyMW; //Vy
			vel[q+2] = vzTran*cos(angleChi)-A*vyTran-A*vxTran*vzTran + vzMW; //Vz
		}else{
			// Charge exchange:
			//msg(STATUS,"ch-exchange");
			angleChi = acos(sqrt(1-R)); // for M_ion = M_neutral
			anglePhi = 2*PI*R; // set different R?
			angleTheta = acos(vxTran); // transfomed vx or not?
			A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
			vel[q] = vxMW*cos(angleChi)+A*(vyMW*vyMW + vzMW*vzMW); //Vx
			vel[q+1] = vyMW*cos(angleChi)+A*vzMW-A*vxMW*vyMW; //Vy
			vel[q+2] = vzMW*cos(angleChi)-A*vyMW-A*vxMW*vzMW; //Vz
		}
	}

	msg(STATUS,"errorcounter = %i should be the same as number of coll particles\
	= %i", (errorcounter),NparticleColl);

	//free?
}



/*************************************************
*		Inline functions
************************************************/



/*************************************************
*		DEFINITIONS
************************************************/



/*************************************************
* 		ALLOCATIONS
* 		DESTRUCTORS
************************************************/


/*******************************************************
*			VARIOUS COMPUTATIONS (RESIDUAL)
******************************************************/



/*************************************************
*		RUNS
************************************************/

funPtr mccTestMode_set(dictionary *ini){
	//test sanity here!
	mccSanity(ini,"mccTestMode",2);
	return mccTestMode;
}

void mccTestMode(dictionary *ini){
	int errorvar = 0;
	msg(STATUS, "start mcc Test Mode");

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



	/*
	* INITIALIZE PINC VARIABLES
	*
	*done in "main" for "regular mode"
	*/
	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini);
	Grid *E   = gAlloc(ini, VECTOR);
	Grid *rho = gAlloc(ini, SCALAR);
	Grid *res = gAlloc(ini, SCALAR);
	Grid *phi = gAlloc(ini, SCALAR);
	Multigrid *mgRho = mgAlloc(ini, rho);
	Multigrid *mgRes = mgAlloc(ini, res);
	Multigrid *mgPhi = mgAlloc(ini, phi);


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
	// oOpenH5(ini, obj, mpiInfo, denorm, dimen, "test");
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
	* compute initial half-step
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

		puMove(pop);
		// oRayTrace(pop, obj);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		pPosAssertInLocalFrame(pop, rho);


		//-----------------------------
		// Compute charge density
		distr(pop, rho);
		//msg(STATUS,"HEEEEEEEEEEERRRRRRRRRRREEEEEEEEEEEE");
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		//---------------------------
		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);

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


		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		/*
		* mcc specific variables
		*/
		double mccTimestep = iniGetDouble(ini,"time:timeStep");
		double frequency = iniGetDouble(ini,"collisions:collisionFrequency");
		double PnullElectron = 1-exp(-(frequency*mccTimestep));
		//msg(STATUS, "freq is %f",frequency);
		msg(STATUS, "Pnull for Electrons is %f",PnullElectron);


		mccCollideElectron(pop, PnullElectron, rng); //race conditions?????
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		mccCollideIon(ini, pop, PnullElectron, rng);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary


	}
	msg(STATUS, "Test returned %d", errorvar);

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
	mgFree(mgRho);
	mgFree(mgPhi);
	mgFree(mgRes);
	gFree(rho);
	gFree(phi);
	gFree(res);
	gFree(E);
	pFree(pop);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}
