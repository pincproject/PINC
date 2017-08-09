/**
* @file		 	collisions.c
* @brief		Collisional module, Montecarlo Collision Method.
* @author		Steffen Mattias Brask <steffen.brask@fys.uio.no>
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
	// check error, see vahedi + surendra p 181
	// check for v_ion + v_neutral > maxVel (moves beyond a cell in one timestep)
	int countSpecies = iniGetInt(ini,"population:nSpecies");
	if(countSpecies!=nSpecies){
		msg(ERROR,"%s only supports 2 species",name);
	}
}

void mccReadcrossect(const char* filename){
	// read crossections from file?
	FILE* fp = fopen(filename, "r" );

	fclose(fp);
}
double mccGetMaxVel(const Population *pop, int species){

	double *vel = pop->vel;
	double MaxVelocity = 0;
	double NewVelocity = 0;
	int nDims = pop->nDims;

	long int iStart = pop->iStart[species];
	long int iStop  = pop->iStop[species];
	for(int i=iStart; i<iStop; i+=nDims){
		NewVelocity = sqrt(vel[i]*vel[i]+vel[i+1]*vel[i+1]+vel[i+2]*vel[i+2]);
		//msg(STATUS,"NewVel = %f ",NewVelocity);
		if(NewVelocity>MaxVelocity){
			MaxVelocity=NewVelocity;
		}
	}
	return MaxVelocity;
}

double mccSigmaCEX(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...
	double temp = 0.7*pow(10,4); // hardcoded order of---
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
double mccSigmaIonElastic(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...
	double temp = 0.3*pow(10,4); // hardcoded order of---
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
double mccSigmaElectronElastic(double eps, double DebyeLength){
	// need a function for calculating real cross sects from data...
	double temp = 0.5*pow(10,5); // hardcoded order of---
	//temp = temp/(DebyeLength*DebyeLength); // convert to PINC dimensions
	return temp;
}
void mccGetPnullIon(double m, double thermalVel, double dt, double nt, double DebyeLength,double *Pnull, double *nullFreq, Population *pop){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	// for now lets do
	double max_v = mccGetMaxVel(pop,1); //2.71828*thermalVel; // e*thermalVel, needs to be max_velocity function
	//msg(STATUS,"maxVelocity =  %f", max_v);
	*nullFreq = (mccSigmaCEX(0.0, DebyeLength)+mccSigmaIonElastic(0.0, DebyeLength))*max_v*nt; //nullFreq = max_eps(sigma_T *v)*nt (or max frequency of colls?)
	//msg(STATUS,"nullFreq =  %f", *nullFreq );
	*Pnull = 1-exp(-(*nullFreq*dt));
	//*Pnull = 0.01;
	//msg(STATUS,"getPnull =  %f", *Pnull);
	//return Pnull, nullFreq;
}
void mccGetPnullElectron(double m, double thermalVel, double dt, double nt, double DebyeLength,double *Pnull, double *nullFreq, Population *pop){
	// determine max local number desity / max_x(n_t(x_i)) (for target species, neutrals)
	// Get å functional form of sigma_T (total crossect as funct. of energy)
	// determine the speed, v(eps_i) = sqrt(2*eps_i *m_s)
	// determine, max_eps(sigma_T *v)

	// for now lets do
	double max_v = mccGetMaxVel(pop,0);//2.71828*thermalVel; // e*thermalVel, needs to be max_velocity
	//msg(STATUS,"maxVelocity =  %f", max_v);
	*nullFreq = mccSigmaElectronElastic(0.0, DebyeLength)*max_v*nt; //nullFreq = max_eps(sigma_T *v)*nt (or max frequency of colls?)
	//msg(STATUS,"nullFreq =  %f", *nullFreq );
	*Pnull = 1-exp(-(*nullFreq*dt));
	//*Pnull = 0.01;
	//msg(STATUS,"getPnull =  %f", *Pnull);
	//return Pnull, nullFreq;
}
double mccGetMyCollFreq(double (*sigma)(double, double), double m, double vx,double vy,double vz, double dt, double nt, double DebyeLength){
	// collision freq given a cross section function
	double v = sqrt(vx*vx + vy*vy + vz*vz);
	double eps_i = 0.5*m*v*v;
	double sigma_T = sigma(eps_i, DebyeLength);
	double MyFreq = v*sigma_T*nt;
	//msg(STATUS,"Myfreq =  %f", MyFreq);
	return MyFreq;
}

double mccEnergyDiffIonElastic(double Ekin, double theta, double mass1, double mass2){
	double temp;
	//msg(STATUS, "Ekin = %.32f", Ekin);
	//msg(STATUS, "theta = %f", theta);
	//msg(STATUS, "mass1 = %f", mass1);
	//msg(STATUS, "mass2 = %f", mass2);
	temp = Ekin*((2*mass1*mass2)/(mass1+mass2))*(1.0-cos(theta));
	//msg(STATUS, "temp = %.32f", temp);
	return temp;
}

void mccCollideElectron(Population *pop,  double *Pnull, double *nullFreqElectron, const gsl_rng *rng, double dt,double nt, double DebyeLength){

	//msg(STATUS,"colliding Electrons");

	int nDims = pop->nDims;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double R = gsl_rng_uniform(rng);
	double Rp = gsl_rng_uniform(rng);

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	double vx, vy, vz;

	//int errorcounter = 0;
	double Ekin = 0;
	long int iStart = pop->iStart[0];		// *nDims??
	long int iStop = pop->iStop[0];      // make shure theese are actually electrons!
	long int NparticleColl = (*Pnull)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	//if ((mccStepSize-1)*(NparticleColl) > (iStop-iStart)){ //errorcheck, remove
	//	msg(WARNING,"particle collisions out of bounds in mccCollideElectron");
	//}
	msg(STATUS,"number of particles in array = %i", (iStop-iStart));
	msg(STATUS,"number of particles colliding = %i", (NparticleColl));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds
	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pnull*(pStop-pStart);p++)
		//errorcounter += 1;
		//msg(STATUS,"errorcounter = %i", (errorcounter));
		//if (errorcounter > NparticleColl){ //errorcheck, remove
		//	msg(WARNING,"errorcounter = %i should be the same as number of coll\
			particles = %i", (errorcounter),NparticleColl);
		//	msg(WARNING,"particle collisions out of bounds\
			in mccCollideElectron");
			//errorcounter = 0;
		//}
		R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
		Rp = gsl_rng_uniform(rng); // separate rand num. for prob.
		q = (i + (R*mccStepSize))*nDims;
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );
		vx = vel[q];
		vy = vel[q+1];
		vz = vel[q+2];
		// need n,n+1,n+2 for x,y,j ?

		double MyCollFreq = mccGetMyCollFreq(mccSigmaElectronElastic, mass[0],vx,vy,vz,dt,nt,DebyeLength); //prob of coll for particle i

		//msg(STATUS,"is Rp = %f < MyCollFreq/nullFreqElectron = %f", Rp, (MyCollFreq/ *nullFreqElectron));

		if (Rp<(MyCollFreq/ *nullFreqElectron)){
		// Pnull gives the max possible colls. but we will "never" reach this
		// number, so we need to test for each particles probability.

			//if (( (MyCollFreq)/ *nullFreqElectron)>1){
				//put in sanity later, need real cross sects
				//msg(WARNING,"MyCollFreqElectron)/ *nullFreqElectron > 1");
			//}
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

	}
	//msg(STATUS,"errorcounter = %i should be the same as number of\
	coll particles = %i", (errorcounter),NparticleColl);

	//free?
}
 //1590430

void mccCollideIon(const dictionary *ini, Population *pop, double *Pnull,
	double *nullFreqIon, const gsl_rng *rng, double timeStep, double nt, double DebyeLength){

	msg(STATUS,"colliding Ions");

	//int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *vel = pop->vel;
	double *mass = pop->mass;
	double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	//double timeStep = iniGetDouble(ini,"time:timeStep");
	double stepSize = iniGetDouble(ini,"grid:stepSize");
	double velTh = (timeStep/stepSize)*NvelThermal; //neutrals
	//free(velThermal);

	double R, R1, Rp;

	double angleChi = 0; //scattering angle in x, y
	double anglePhi = 0; //scattering angle in j
	double angleTheta = 0;
	double A = 0; //scaling Factor used for precomputing (speedup)
	long int q = 0;
	//double fraction = 0.7; //defines fraction of particles colliding w el/ch-ex
	double vxMW, vyMW, vzMW;
	//int errorcounter = 0;
	double Ekin = 0;
	double EkinAfter = 0;
	double EnergyDiff = 1;

	long int iStart = pop->iStart[1];			// *nDims??
	long int iStop = pop->iStop[1];      // make shure theese are actually ions!
	long int NparticleColl = (*Pnull)*(pop->iStop[0]-pop->iStart[0]); // 1 dim ?, Number of particles too collide
	long int mccStepSize = floor((double)((iStop-iStart))\
	/ (double)(NparticleColl));
	//if ((mccStepSize-1)*(NparticleColl)/nDims > (iStop-iStart)){ //errorcheck, remove
	//	msg(WARNING,"particle collisions out of bounds in mccCollideIon");
	//	msg(WARNING,"%i is bigger than array = %i",(mccStepSize-1)\
	//	*(NparticleColl),(iStop-iStart) );
	//}
	msg(STATUS,"number of particles in array = %i", (iStop-iStart));
	msg(STATUS,"number of particles colliding = %i", (NparticleColl));
	msg(STATUS,"Particles per box to pick one from = %i", (mccStepSize));

	long int mccStop = iStart + mccStepSize*NparticleColl; //enshure non out of bounds
	for(long int i=iStart;i<mccStop;i+=mccStepSize){  //check this statement #segfault long int p=pStart;p<pStart+Pnull*(pStop-pStart);p++)
		//errorcounter += 1;
		//if (errorcounter > NparticleColl){ //errorcheck, remove
		//	msg(WARNING,"errorcounter = %i should be the same as number of coll\
		//	particles = %i", (errorcounter),NparticleColl/nDims);
		//	msg(WARNING,"particle collisions out of bounds in mccCollideIon");
			//errorcounter = 0;
		//}
		Rp = gsl_rng_uniform(rng); //decides type of coll.
		R = gsl_rng_uniform(rng); // New random number per particle. maybe print?
		q = (i + (R*mccStepSize))*nDims;
		//msg(STATUS,"R=%g, p = %i, q = %i", R,p,q );

		vxMW = gsl_ran_gaussian_ziggurat(rng,velTh); //maxwellian dist?
		vyMW = gsl_ran_gaussian_ziggurat(rng,velTh);
		vzMW = gsl_ran_gaussian_ziggurat(rng,velTh);

		double vxTran = vel[q]-vxMW;   //simple transfer, should be picked from maxwellian
		double vyTran = vel[q+1]-vyMW;  //vel-velMaxwellian
		double vzTran = vel[q+2]-vzMW;  // this will for now break conservation laws???

		double MyCollFreq1 = mccGetMyCollFreq(mccSigmaIonElastic, mass[1],vxTran,vyTran,vzTran,timeStep,nt,DebyeLength); //prob of coll for particle i
		double MyCollFreq2 = mccGetMyCollFreq(mccSigmaCEX, mass[1],vxTran,vyTran,vzTran,timeStep,nt,DebyeLength); //duplicate untill real cross sections are implemented

		if (Rp<( (MyCollFreq1+MyCollFreq2)/ *nullFreqIon)){ // MyCollFreq1+MyCollFreq2 defines total coll freq.

			if (( (MyCollFreq1+MyCollFreq2)/ *nullFreqIon)>1){
				//put in sanity later, need real cross sects
				msg(WARNING,"(MyCollFreq1+MyCollFreq2)/ *nullFreqIon > 1");
			}
			if(Rp < MyCollFreq1/ *nullFreqIon){ // if test slow, but enshures randomness...
				// elastic:

				//msg(STATUS,"elastic");


				//msg(ERROR,"Elastic collision"); // temporary test
				Ekin = 0.5*(vel[q]*vel[q] + vel[q+1]*vel[q+1] + vel[q+2]*vel[q+2])*mass[0];
				R1 = gsl_rng_uniform(rng);
				angleChi = acos(sqrt(1-R)); // for M_ion = M_neutral
				anglePhi = 2*PI*R1; // set different R?
				angleTheta = 2*angleChi;
				//EnergyDiff = mccEnergyDiffIonElastic(Ekin, angleTheta, mass[0], mass[1]);
				A = (sin(angleChi)*sin(anglePhi))/sin(angleTheta);
				vel[q] = vxTran*cos(angleChi)+A*(vyTran*vyTran + vzTran*vzTran) \
				+ vxMW; //Vx
				vel[q+1] = vyTran*cos(angleChi)+A*vzTran-A*vxTran*vyTran + vyMW; //Vy
				vel[q+2] = vzTran*cos(angleChi)-A*vyTran-A*vxTran*vzTran + vzMW; //Vz
				EkinAfter = 0.5*(vel[q]*vel[q] + vel[q+1]*vel[q+1] + vel[q+2]*vel[q+2])*mass[0];
				//msg(STATUS, "---------- energydiff = %.32f", (Ekin-EkinAfter) );
				//msg(STATUS, "analytical energydiff = %.32f",EnergyDiff );
			}else{
				// Charge exchange:
				// flip ion and neutral. Neutral is new ion.
				//msg(STATUS,"ch-ex");
				//errorcounter += 1;
				vel[q] = vxMW; //Vx
				vel[q+1] = vyMW; //Vy
				vel[q+2] = vzMW; //Vz
			}
		}
	}

	//msg(STATUS,"errorcounter = %i should be the same as number of\
	coll particles = %i", (errorcounter),NparticleColl);

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
	puAccND0KE_set,
	puBoris3D1_set,
	puBoris3D1KE_set,
	puBoris3D1KETEST_set);

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

	/*
	* mcc specific variables
	*/
	// make struct???
	double PnullElectron = 0;
	double nullFreqElectron = 0;
	double PnullIon = 0;
	double nullFreqIon = 0;
	double *mass = pop->mass;
	int nSpecies = pop->nSpecies;
	double nt = iniGetDouble(ini,"collisions:numberDensityNeutrals"); //constant for now
	double DebyeLength = iniGetDouble(ini,"grid:debye"); // for converting crosssection dimensions
	double mccTimestep = iniGetDouble(ini,"time:timeStep");
	//double frequency = iniGetDouble(ini,"collisions:collisionFrequency"); // for static Pnull...
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);

	// using Boris algo
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));
	puGet3DRotationParameters(ini, T, S);

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

	//pOpenH5(ini, pop, "pop");
	//gOpenH5(ini, rho, mpiInfo, denorm, dimen, "rho");
	gOpenH5(ini, phi, mpiInfo, denorm, dimen, "phi");
	//gOpenH5(ini, E,   mpiInfo, denorm, dimen, "E");
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

	//msg(STATUS, "constants = %f, %f", velThermal[0], velThermal[1]);
	// Initalize particles
	//pPosLattice(ini, pop, mpiInfo);
	pPosUniform(ini, pop, mpiInfo, rngSync);
	//pVelZero(pop);
	//pVelConstant(ini, pop, velThermal[0], velThermal[1]); //constant values for vel.
	pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");

	// Perturb particles
	// pPosPerturb(ini, pop, mpiInfo);

	// compute Pnull once outside loop. More correct is every dt
	mccGetPnullElectron(mass[0], velThermal[0], mccTimestep, nt, DebyeLength, &PnullElectron, &nullFreqElectron, pop);
	mccGetPnullIon(mass[1], velThermal[1], mccTimestep, nt, DebyeLength, &PnullIon, &nullFreqIon, pop);
	//msg(STATUS, "freq is %f",frequency);

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

	// add External E
	puAddEext(ini, pop, E);


	// Advance velocities half a step
	gMul(E, 0.5);
	for(int d = 0; d < nSpecies*3;d++) T[d] *= 0.5;
	for(int d = 0; d < nSpecies*3;d++) S[d] *= 0.5;
	acc(pop, E, T, S);
	gMul(E, 2.0);
	//msg(STATUS, "S = %f", T);
	for(int d = 0; d < nSpecies*3;d++) T[d] *= 2.0;
	for(int d = 0; d < nSpecies*3;d++) S[d] *= 2.0;

	//Write initial h5 files
	//gWriteH5(E, mpiInfo, 0.0);
	//gWriteH5(rho, mpiInfo, 0.0);
	gWriteH5(phi, mpiInfo, 0.0);
	//pWriteH5(pop, mpiInfo, 0.0, 0.5);
	pWriteEnergy(history,pop,0.0);



	int collsOnOff;
	collsOnOff = iniGetInt(ini,"collisions:collisionsOnOff");
	if(collsOnOff == 0){
		PnullElectron = 0;
		PnullIon = 0;
		msg(STATUS, "----COLLISIONS ARE TURNED OFF----");
	}

	double debugVel = 0;

	msg(STATUS, "Pnull for Electrons is %f",PnullElectron);
	msg(STATUS, "Pnull for Ions is %f",PnullIon);

	/*
	* TIME LOOP
	*/

	Timer *t = tAlloc(rank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		msg(STATUS,"Computing time-step %i",n);
		debugVel = mccGetMaxVel(pop,0);
		//msg(STATUS,"maxVelocity = %f",debugVel);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);


		// Move particles
		//msg(STATUS, "moving particles");
		puMove(pop);
		// oRayTrace(pop, obj);

		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		mccCollideElectron(pop, &PnullElectron, &nullFreqElectron, rng, mccTimestep, nt,DebyeLength); //race conditions?????
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary
		mccCollideIon(ini, pop, &PnullIon, &nullFreqIon, rng, mccTimestep, nt, DebyeLength);
		MPI_Barrier(MPI_COMM_WORLD);	// Temporary, shouldn't be necessary

		// Migrate particles (periodic boundaries)
		//msg(STATUS, "extracting emigrants");
		extractEmigrants(pop, mpiInfo);
		//msg(STATUS, "Migrating particles");
		puMigrate(pop, mpiInfo, rho);

		// Check that no particle resides out-of-bounds (just for debugging)
		//msg(STATUS, "checking particles out of bounds");
		pPosAssertInLocalFrame(pop, rho);


		//-----------------------------
		// Compute charge density
		//msg(STATUS, "computing charge density");
		distr(pop, rho);
		//msg(STATUS, "MPI exchanging density");
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);

		//---------------------------
		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(rho, mpiInfo);

		// Compute electric potential phi
		//msg(STATUS, "Compute electric potential phi");
		solve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(phi, mpiInfo);

		// Compute E-field
		//msg(STATUS, "compute E field");
		gFinDiff1st(phi, E);
		//msg(STATUS, "MPI exchanging E");
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		//msg(STATUS, "gMul( E, -1)");
		gMul(E, -1.);

		//msg(STATUS, "gAssertNeutralGrid");
		gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		// gAddTo(Ext);
		puAddEext(ini, pop, E);

		// Accelerate particle and compute kinetic energy for step n
		//msg(STATUS, "Accelerate particles");
		acc(pop, E, T, S);
		tStop(t);

		// Sum energy for all species
		//msg(STATUS, "sum kinetic energy");
		pSumKinEnergy(pop);
		//msg(STATUS, "compute pot energy");
		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		//msg(STATUS, "writing to file");
		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);
		//Write h5 files
		//gWriteH5(E, mpiInfo, (double) n);
		//gWriteH5(rho, mpiInfo, (double) n);
		gWriteH5(phi, mpiInfo, (double) n);
		//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		pWriteEnergy(history,pop,(double)n);

		//free(nt);
		//free(mccTimestep);
		//free(frequency);
		//free(velThermal);
		msg(STATUS, "   -    ");


	}
	//msg(STATUS, "Test returned %d", errorvar);

	if(mpiInfo->mpiRank==0) tMsg(t->total, "Time spent: ");

	/*
	* FINALIZE PINC VARIABLES
	*/
	gFreeMpi(mpiInfo);

	// Close h5 files
	//pCloseH5(pop);
	//gCloseH5(rho);
	gCloseH5(phi);
	//gCloseH5(E);
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
	free(S);
	free(T);
	// oFree(obj);

	gsl_rng_free(rngSync);
	gsl_rng_free(rng);

}
