/**
 * @file		population.test.c
 * @brief		Unit tests for population.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */

#include "core.h"
#include "test.h"
#include <math.h>

static int testPCut(){

	dictionary *ini = iniGetDummy();

	iniparser_set(ini,"population:nAlloc","10,10");
	iniparser_set(ini,"population:nParticles","0,0");
	iniparser_set(ini,"population:q","-1,1");
	iniparser_set(ini,"population:m","1,100");
	Population *pop = pAlloc(ini);

	double posV[] = {0,1,2};
	double velV[] = {0,10,20};
	pNew(pop,1,posV,velV);

	adSet(posV,3,3.,4.,5.);
	adSet(velV,3,30.,40.,50.);
	pNew(pop,1,posV,velV);

	adSet(posV,3,6.,7.,8.);
	adSet(velV,3,60.,70.,80.);
	pNew(pop,1,posV,velV);

	adSet(posV,3,9.,10.,11.);
	adSet(velV,3,90.,100.,110.);
	pNew(pop,1,posV,velV);

	pCut(pop,1,33,posV,velV);

	double expected[] = {3,4,5};
	utAssert(adEq(posV,expected,3,pow(10,-14)),"Particle position extracted incorrectly");
	adSet(expected,3,30.,40.,50.);
	utAssert(adEq(velV,expected,3,pow(10,-14)),"Particle velocity extracted incorrectly");

	utAssert(pop->iStop[1]==13,"Particle counter not properly updated");

	pCut(pop,1,33,posV,velV);

	adSet(expected,3,9.,10.,11.);
	utAssert(adEq(posV,expected,3,pow(10,-14)),"Particle position fill-in malfunctioning");
	adSet(expected,3,90.,100.,110.);
	utAssert(adEq(velV,expected,3,pow(10,-14)),"Particle velocity fill-in malfunctioning");

	utAssert(pop->iStop[1]==12,"Particle counter not properly updated");

	return 0;

}

// All tests for io.c is contained in this function
void testPopulation(){
	utRun(&testPCut);
}
