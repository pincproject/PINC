/**
 * @file		pusher.test.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 * @copyright	University of Oslo, Norway
 * @brief		Unit tests for pusher.c
 * @date		17.12.15
 */

#include "test.h"
#include "pinc.h"
#include "pusher.h"
#include <math.h>

/*
 * Test the acceleration of particles in constant E-field in x-direction. Tests
 * three particles with varying q and m to test specie-specific normalization
 * used in pusher. Since E is constant it is not a good test for the
 * interpolation algorithm.
 */
static int testConstE(){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","10,10,10");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","1,1,-1");
	iniparser_set(ini,"population:m","1,2,1");
	iniparser_set(ini,"grid:trueSize","256,256,256");
	iniparser_set(ini,"grid:stepSize","1,1,1");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");
	iniparser_set(ini,"time:timestep","1");

	// Assign population
	Population *pop = pAlloc(ini);
	double *pos = pop->pos;

	double posV[] = {100,100,100};
	double velV[] = {0,0,0};

	pNew(pop,0,posV,velV);
	pNew(pop,1,posV,velV);
	pNew(pop,2,posV,velV);

	// Assign grid quantities
	Grid *E = gAlloc(ini,3);
	double val[] = {1,0,0};
	gSet(E,val);
	gNormalizeE(ini,E);	// Normalization should yield 1 too. Thus testing this as well.

	// Accelerate half-step
	gMul(E,0.5);
	puAcc3D1(pop,E);
	gMul(E,2);

	int N = 5;
	for(int n=1;n<=N;n++){

		puMove(pop);
		puAcc3D1(pop,E);

		double ana;

		ana = 100 + 0.5*pow(n,2);
		utAssert( fabs(ana-pos[0])  < pow(10,-15),
			"Numeric and analytical results does not match. Dir: 0, Specie: 0, step: %i, ana: %f, num: %f",n,ana,pos[0]);

		ana = 100 + 0.25*pow(n,2);
		utAssert( fabs(ana-pos[30]) < pow(10,-15),
			"Numeric and analytical results does not match. Dir: 0, Specie: 1, step: %i, ana: %f, num: %f",n,ana,pos[30]);

		ana = 100 - 0.5*pow(n,2);
		utAssert( fabs(ana-pos[60]) < pow(10,-15),
			"Numeric and analytical results does not match. Dir: 0, Specie: 2, step: %i, ana: %f, num: %f",n,ana,pos[60]);

	}

	return 0;

}

/*
 * Tests puAcc3D1 and with it puInterp3D1. Only one specie, specie-specific
 * renormalization not tested. See testConstE().
 */
static int testPuAcc3D1(){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","10,10,10");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","1,1,-1");
	iniparser_set(ini,"population:m","1,2,1");
	iniparser_set(ini,"time:timeStep","1");
	iniparser_set(ini,"grid:stepSize","1,1,1");
	iniparser_set(ini,"grid:trueSize","5,4,3");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");

	Grid *grid = gAlloc(ini,3);
	for(int p=0;p<grid->sizeProd[grid->rank];p++) grid->val[p] = p;

	Population *pop = pAlloc(ini);
	double *vel = pop->vel;

	double velV[] = {100,100,100}; // Non-zero to test that v+=dv and not v=dv

	double posV[] = {2.5,1.5,0.5};
	pNew(pop,0,posV,velV);

	posV[0] = 0.1;
	posV[1] = 0.2;
	posV[2] = 0.3;
	pNew(pop,0,posV,velV);

	puAcc3D1(pop,grid);

	// Specie 0, particle 0, center in cell
	utAssert( fabs( vel[0]-160 ) < pow(10,-13), "Centered interpolation failed, x-component");
	utAssert( fabs( vel[1]-161 ) < pow(10,-13), "Centered interpolation failed, y-component");
	utAssert( fabs( vel[2]-162 ) < pow(10,-13), "Centered interpolation failed, z-component");

	// Specie 0, particle 1, non-centered
	utAssert( fabs( vel[3]-121.3 ) < pow(10,-13), "Non-centered interpolation failed");

	return 0;
}

static int testPuDistr3D1(){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","10,10,10");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","1,1,-1");
	iniparser_set(ini,"population:m","1,2,1");
	iniparser_set(ini,"time:timeStep","3");
	iniparser_set(ini,"grid:stepSize","2,2,2");
	iniparser_set(ini,"grid:trueSize","5,4,3");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");

	// Compute normalization factor
	int nDims;
	double timeStep = iniparser_getdouble(ini,"time:timeStep",0.0);
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",&nDims);
	double cellVolume = adProd(stepSize,nDims);
	double norm = pow(timeStep,2)/cellVolume;


	Grid *rho = gAlloc(ini,1);
	gZero(rho);
	double *val = rho->val;

	Population *pop = pAlloc(ini);
	double velV[] = {100,100,100};

	// One centere specie 0 particle
	double posV[] = {0.5,0.5,0.5};
	pNew(pop,0,posV,velV);

	// One non-centered specie 0 particle
	posV[0] = 2.2;
	posV[1] = 0.3;
	posV[2] = 0.4;
	pNew(pop,0,posV,velV);

	// Two specie 0 particles
	posV[0] = 0.5;
	posV[1] = 2.5;
	posV[2] = 0.5;
	pNew(pop,0,posV,velV);
	posV[0] = 1.2;
	pNew(pop,0,posV,velV);

	puDistr3D1(pop,rho);

	utAssert( fabs( val[0] -norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[1] -norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[5] -norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[6] -norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[20]-norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[21]-norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[25]-norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");
	utAssert( fabs( val[26]-norm*0.125 ) < pow(10,-13), "Distribution of one centered specie 0 particle failed");

	utAssert( fabs( val[2] -norm*0.336 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[3] -norm*0.084 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[7] -norm*0.144 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[8] -norm*0.036 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[22]-norm*0.224 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[23]-norm*0.056 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[27]-norm*0.096 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");
	utAssert( fabs( val[28]-norm*0.024 ) < pow(10,-13), "Distribution of one non-centered specie 0 particle failed");

	utAssert( fabs( val[10]-norm*0.125 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[15]-norm*0.125 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[30]-norm*0.125 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[35]-norm*0.125 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[11]-norm*0.325 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[16]-norm*0.325 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[31]-norm*0.325 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[36]-norm*0.325 ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[12]-norm*0.05  ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[17]-norm*0.05  ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[32]-norm*0.05  ) < pow(10,-13), "Distribution of two specie 0 particles failed");
	utAssert( fabs( val[37]-norm*0.05  ) < pow(10,-13), "Distribution of two specie 0 particles failed");

	return 0;

}

// Multi-specie test utilizing specie renormalization
static int testPuDistr3D1renorm(){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","10,10,10");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","-1,1,2");
	iniparser_set(ini,"population:m","10,1,10");
	iniparser_set(ini,"time:timeStep","3");
	iniparser_set(ini,"grid:stepSize","2,2,2");
	iniparser_set(ini,"grid:trueSize","5,4,3");
	iniparser_set(ini,"grid:nGhostLayers","0,0,0,0,0,0");

	// Compute normalization factor
	int nDims, nSpecies;
	double timeStep = iniparser_getdouble(ini,"time:timeStep",0.0);
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",&nDims);
	double *q = iniGetDoubleArr(ini,"population:q",&nSpecies);
	double *m = iniGetDoubleArr(ini,"population:m",&nSpecies);
	double cellVolume = adProd(stepSize,nDims);
	double norm = (q[0]/m[0])*pow(timeStep,2)/cellVolume;


	Grid *rho = gAlloc(ini,1);
	gZero(rho);
	double *val = rho->val;

	Population *pop = pAlloc(ini);
	double velV[] = {100,100,100};

	double posV[] = {0.2,0.2,0.5};
	pNew(pop,0,posV,velV);

	posV[0] = 0.8;
	pNew(pop,1,posV,velV);

	posV[0] = 0.2;
	posV[1] = 0.8;
	pNew(pop,2,posV,velV);

	puDistr3D1(pop,rho);

	utAssert( fabs( val[0] +norm*0.08 ) < pow(10,-13), "Distribution of multiple species failed 1");
	utAssert( fabs( val[20]+norm*0.08 ) < pow(10,-13), "Distribution of multiple species failed 2");
	utAssert( fabs( val[1] -norm*0.28 ) < pow(10,-13), "Distribution of multiple species failed 3");
	utAssert( fabs( val[21]-norm*0.28 ) < pow(10,-13), "Distribution of multiple species failed 4");
	utAssert( fabs( val[5] -norm*0.58 ) < pow(10,-13), "Distribution of multiple species failed 5");
	utAssert( fabs( val[25]-norm*0.58 ) < pow(10,-13), "Distribution of multiple species failed 6");
	utAssert( fabs( val[6] -norm*0.22 ) < pow(10,-13), "Distribution of multiple species failed 7");
	utAssert( fabs( val[26]-norm*0.22 ) < pow(10,-13), "Distribution of multiple species failed 8");

	return 0;

}

static int testPuBndIdMigrantsXD(){

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","100,100,100");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","-1,1,2");
	iniparser_set(ini,"population:m","10,1,10");
	iniparser_set(ini,"grid:trueSize","8,8,8");
	iniparser_set(ini,"grid:nGhostLayers","1,1,1,1,1,1");
	iniparser_set(ini,"grid:thresholds","1,1,1,-1,-1,-1");
	iniparser_set(ini,"grid:nEmigrantsAlloc","10");

	Population *pop = pAlloc(ini);
	Grid *grid = gAlloc(ini,1);
	MpiInfo *mpiInfo = gAllocMpi(ini);
	gCreateNeighborhood(ini,mpiInfo,grid);

	double vel[] = {0,0,0};
	double pos[] = {0,5,5};

	// Placing a line of particles along x in the middle of the yz-plane
	// Placing two species to test that as well
	for(;pos[0]<=10;pos[0]+=0.5){
		pNew(pop,0,pos,vel);
		pNew(pop,1,pos,vel);
	}

	// Placing a particle in each region
	for(int z=-1;z<=+1;z++){
		for(int y=-1;y<=+1;y++){
			for(int x=-1;x<=+1;x++){
				adSet(pos,3,5+x*4.5,5+y*4.5,5+z*4.5);
				pNew(pop,0,pos,vel);
				pNew(pop,1,pos,vel);
			}
		}
	}

	long int **migrants = mpiInfo->migrants;
	long int *nEmigrants = mpiInfo->nEmigrants;

	long int *result = malloc(10*sizeof(*result));

	// Specifying number of particles
	long int *nEmigrantsResult = malloc(81*sizeof(*nEmigrantsResult));
	alSet(nEmigrantsResult,81,	1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,
								1,1,0,1,1,0,1,1,0,3,3,0,0,0,0,4,4,0,1,1,0,1,1,0,1,1,0,
								1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0);

	// Specifying start of each new species instead. More in line with PINC practice.
	// long int *nMigrantsResult = malloc(108*sizeof(*nMigrantsResult));
	// alSet(nMigrantsResult,108,	0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2,
	// 							0,1,2,2, 0,1,2,2, 0,1,2,2, 0,3,6,6, 0,0,0,0, 0,4,8,8, 0,1,2,2, 0,1,2,2, 0,1,2,2,
	// 							0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2, 0,1,2,2,);

	// 3D

	puBndIdMigrants3D(pop,mpiInfo);

	utAssert(alEq(nEmigrants,nEmigrantsResult,81),"Wrong count of particles migrated to each domain (3D-method)");

	alSet(result,2,63,363);	// don't care about the remaining 8 elements
	for(int i=0;i<=11;i++){
		utAssert(alEq(migrants[i],result,2),"Wrong migrants[%i] (3D-method)",i);
		alShift(result,2,3);
	}
	alShift(result,2,9);
	for(int i=15;i<=26;i++){
		utAssert(alEq(migrants[i],result,2),"Wrong migrants[%i] (3D-method)",i);
		alShift(result,2,3);
	}
	alSet(result,6,0,3,99,300,303,399);
	utAssert(alEq(migrants[12],result,6),"Wrong migrants[12] (3D-method)");
	alSet(result,8,54,57,60,105,354,357,360,405);
	utAssert(alEq(migrants[14],result,8),"Wrong migrants[14] (3D-method)");

	// ND

	puBndIdMigrantsND(pop,mpiInfo);

	utAssert(alEq(nEmigrants,nEmigrantsResult,81),"Wrong count of particles migrated to each domain (ND-method)");

	alSet(result,2,63,363);	// don't care about the remaining 8 elements
	for(int i=0;i<=11;i++){
		utAssert(alEq(migrants[i],result,2),"Wrong migrants[%i] (ND-method)",i);
		alShift(result,2,3);
	}
	alShift(result,2,9);
	for(int i=15;i<=26;i++){
		utAssert(alEq(migrants[i],result,2),"Wrong migrants[%i] (ND-method)",i);
		alShift(result,2,3);
	}
	alSet(result,6,0,3,99,300,303,399);
	utAssert(alEq(migrants[12],result,6),"Wrong migrants[12] (ND-method)");
	alSet(result,8,54,57,60,105,354,357,360,405);
	utAssert(alEq(migrants[14],result,8),"Wrong migrants[14] (ND-method)");

	return 0;
}

static int testExtractEmigrantsXD(){

	// CREATING INI AND STRUCTS

	dictionary *ini = iniGetDummy();
	iniparser_set(ini,"population:nAlloc","100,100,100");
	iniparser_set(ini,"population:nParticles","1,1,1");
	iniparser_set(ini,"population:q","-1,1,2");
	iniparser_set(ini,"population:m","10,1,10");
	iniparser_set(ini,"grid:trueSize","8,8,8");
	iniparser_set(ini,"grid:nGhostLayers","1,1,1,1,1,1");
	iniparser_set(ini,"grid:thresholds","1,1,1,-1,-1,-1");
	iniparser_set(ini,"grid:nEmigrantsAlloc","10");

	Population *pop = pAlloc(ini);
	Grid *grid = gAlloc(ini,1);
	MpiInfo *mpiInfo = gAllocMpi(ini);
	gCreateNeighborhood(ini,mpiInfo,grid);

	// PLACE PARTICLES

	double vel[] = {1,2,3};
	double pos[] = {0,5,5};

	// Placing a line of particles along x in the middle of the yz-plane
	// Placing two species to test that as well
	for(;pos[0]<=10;pos[0]+=0.5){
		pNew(pop,0,pos,vel);
		pNew(pop,1,pos,vel);
	}

	// Placing a particle in each region
	for(int z=-1;z<=+1;z++){
		for(int y=-1;y<=+1;y++){
			for(int x=-1;x<=+1;x++){
				adSet(pos,3,5+x*4.5,5+y*4.5,5+z*4.5);
				pNew(pop,0,pos,vel);
				pNew(pop,1,pos,vel);
			}
		}
	}

	double **emigrants = mpiInfo->emigrants;
	long int *nEmigrants = mpiInfo->nEmigrants;

	// Specifying number of particles
	long int *nEmigrantsResult = malloc(81*sizeof(*nEmigrantsResult));
	alSet(nEmigrantsResult,81,	1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,
								1,1,0,1,1,0,1,1,0,3,3,0,0,0,0,4,4,0,1,1,0,1,1,0,1,1,0,
								1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0);

	// 3D


	Timer *t = tAlloc(0);
	puExtractEmigrants3D(pop,mpiInfo);
	tMsg(t,"3D");

	//
	// int ne=14;
	// alPrint(&mpiInfo->nEmigrants[3*ne],3);
	// adPrint(mpiInfo->emigrants[ne],6);
	// adPrint(mpiInfo->emigrants[ne]+6,6);
	// adPrint(mpiInfo->emigrants[ne]+12,6);
	// adPrint(mpiInfo->emigrants[ne]+18,6);
	// adPrint(mpiInfo->emigrants[ne]+24,6);
	// adPrint(mpiInfo->emigrants[ne]+30,6);
	// adPrint(mpiInfo->emigrants[ne]+36,6);
	// adPrint(mpiInfo->emigrants[ne]+42,6);

	utAssert(alEq(nEmigrants,nEmigrantsResult,81),"Wrong count of particles migrated to each domain (3D-method)");

	double tol = 0;	// It's just a copy. No computation. Should be truly equal.

	double *result=malloc(48*sizeof(*result));
	adSet(result,36,0.0,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.0,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.);
	utAssert(adEq(emigrants[12],result,36,tol),"Wrong migrants[12] (3D-method)");

	adSet(result,48,9.5,5.,5.,1.,2.,3.,	// Note: Back-fill will shuffle particle order
					10.,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					9.0,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					10.,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					9.0,5.,5.,1.,2.,3.);
	utAssert(adEq(emigrants[14],result,48,tol),"Wrong migrants[14] (3D-method)");

	for(int z=-1;z<=+1;z++){
		for(int y=-1;y<=+1;y++){
			for(int x=-1;x<=+1;x++){
				int ne=(x+1)+(y+1)*3+(z+1)*9;
				if(ne<12 || ne>14){
					adSet(result,12,5+x*4.5,5+y*4.5,5+z*4.5,1.,2.,3.,
									5+x*4.5,5+y*4.5,5+z*4.5,1.,2.,3.);
					utAssert(adEq(emigrants[ne],result,12,tol),"Wrong migrants[%i] (3D-method)",ne);
				}
			}
		}
	}

	long int iStopExpected[] = {17,117,200};
	utAssert(alEq(pop->iStop,iStopExpected,3),"Wrong number of particles left after extraction");

	adSet(result,3,1.,5.,5.);
	for(int p=0;p<17*3;p+=3){
		if(p==0) result[0] = 5;		// Results get shuffled a bit due to back-fill
		if(p==3) result[0] = 8.5;
		if(p>=6) result[0] = (p/3.0-2)*0.5+1;
		utAssert(adEq(&pop->pos[p],result,3,tol),"Wrong particles left after extraction");
		utAssert(adEq(&pop->vel[p],vel,3,tol),"Wrong particles left after extraction");
		utAssert(adEq(&pop->pos[p+300],result,3,tol),"Wrong particles left after extraction");
		utAssert(adEq(&pop->vel[p+300],vel,3,tol),"Wrong particles left after extraction");
		result[0] += 0.5;
	}


	// ND

	// Deleting old bullshit for particle struct
	pop->iStop[0] = 0;
	pop->iStop[1] = 100;
	pop->iStop[2] = 200;


	// Placing a line of particles along x in the middle of the yz-plane
	// Placing two species to test that as well
	adSet(pos,3,0.,5.,5.);
	for(;pos[0]<=10;pos[0]+=0.5){
		pNew(pop,0,pos,vel);
		pNew(pop,1,pos,vel);
	}
	// Placing a particle in each region
	for(int z=-1;z<=+1;z++){
		for(int y=-1;y<=+1;y++){
			for(int x=-1;x<=+1;x++){
				adSet(pos,3,5+x*4.5,5+y*4.5,5+z*4.5);
				pNew(pop,0,pos,vel);
				pNew(pop,1,pos,vel);
			}
		}
	}

	// ND

	tMsg(t,NULL);
	puExtractEmigrantsND(pop,mpiInfo);
	tMsg(t,"ND");
	//
	// int ne=14;
	// alPrint(&mpiInfo->nEmigrants[3*ne],3);
	// adPrint(mpiInfo->emigrants[ne],6);
	// adPrint(mpiInfo->emigrants[ne]+6,6);
	// adPrint(mpiInfo->emigrants[ne]+12,6);
	// adPrint(mpiInfo->emigrants[ne]+18,6);
	// adPrint(mpiInfo->emigrants[ne]+24,6);
	// adPrint(mpiInfo->emigrants[ne]+30,6);
	// adPrint(mpiInfo->emigrants[ne]+36,6);
	// adPrint(mpiInfo->emigrants[ne]+42,6);

	utAssert(alEq(nEmigrants,nEmigrantsResult,81),"Wrong count of particles migrated to each domain (ND-method)");

	adSet(result,36,0.0,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.0,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.,
					0.5,5.,5.,1.,2.,3.);
	utAssert(adEq(emigrants[12],result,36,tol),"Wrong migrants[12] (ND-method)");

	adSet(result,48,9.5,5.,5.,1.,2.,3.,	// Note: Back-fill will shuffle particle order
					10.,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					9.0,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					10.,5.,5.,1.,2.,3.,
					9.5,5.,5.,1.,2.,3.,
					9.0,5.,5.,1.,2.,3.);
	utAssert(adEq(emigrants[14],result,48,tol),"Wrong migrants[14] (ND-method)");

	for(int z=-1;z<=+1;z++){
		for(int y=-1;y<=+1;y++){
			for(int x=-1;x<=+1;x++){
				int ne=(x+1)+(y+1)*3+(z+1)*9;
				if(ne<12 || ne>14){
					adSet(result,12,5+x*4.5,5+y*4.5,5+z*4.5,1.,2.,3.,
									5+x*4.5,5+y*4.5,5+z*4.5,1.,2.,3.);
					utAssert(adEq(emigrants[ne],result,12,tol),"Wrong migrants[%i] (ND-method)",ne);
				}
			}
		}
	}

	utAssert(alEq(pop->iStop,iStopExpected,3),"Wrong number of particles left after extraction (ND)");

	adSet(result,3,1.,5.,5.);
	for(int p=0;p<17*3;p+=3){
		if(p==0) result[0] = 5;		// Results get shuffled a bit due to back-fill
		if(p==3) result[0] = 8.5;
		if(p>=6) result[0] = (p/3.0-2)*0.5+1;
		utAssert(adEq(&pop->pos[p],result,3,tol),"Wrong particles left after extraction (ND)");
		utAssert(adEq(&pop->vel[p],vel,3,tol),"Wrong particles left after extraction (ND)");
		utAssert(adEq(&pop->pos[p+300],result,3,tol),"Wrong particles left after extraction (ND)");
		utAssert(adEq(&pop->vel[p+300],vel,3,tol),"Wrong particles left after extraction (ND)");
		result[0] += 0.5;
	}

	return 0;
}


// Test conversion between rank and neighbor
static int testPuRankNeighbor(){

	dictionary *ini = iniGetDummy();
	Grid *grid = gAlloc(ini,1);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	aiSet(mpiInfo->nSubdomains,3,5,4,3);
	aiSet(mpiInfo->nSubdomainsProd,4,1,5,20,60);
	aiSet(mpiInfo->subdomain,3,4,0,1);
	mpiInfo->mpiRank = 24;

	gCreateNeighborhood(ini,mpiInfo,grid);

	utAssert(puNeighborToRank(mpiInfo,12)==23,"puNeighborToRank malfunctioning (without wrap-around)");
	utAssert(puNeighborToRank(mpiInfo,2)==15,"puNeighborToRank malfunctioning (with wrap-around)");

	utAssert(puRankToNeighbor(mpiInfo,23)==12,"puRankToNeighbor malfunctioning (without wrap-around)");
	utAssert(puRankToNeighbor(mpiInfo,15)==2,"puRankToNeighbor malfunctioning (with wrap-around)");

	return 0;

}

// All tests for pusher.c is contained in this function
void testPusher(){
	utRun(&testPuAcc3D1);
	utRun(&testPuDistr3D1);
	utRun(&testPuDistr3D1renorm);
	utRun(&testConstE);
	utRun(&testPuBndIdMigrantsXD);
	utRun(&testExtractEmigrantsXD);
	utRun(&testPuRankNeighbor);
}
