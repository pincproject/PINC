/**
 * @file		units.c
 * @brief		Handles units and normalization.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 *
 * Core module for handling of units and normalization.
 */

#include "core.h"
#include "units.h"

const double elementaryCharge = 1.60217733e-19;    // [C]
const double electronMass = 9.10938188e-31;        // [kg]
const double vacuumPermittivity = 8.854187817e-12; // [F/m]

/**
 * @brief					Adds derived units
 * @param[in,out]	units	Pointer to Units
 *
 * Given that length, time, charge, mass and nDims is correctly set this
 * function will derive the rest of the units based on that.
 */
static void uAddDerivedUnits(Units *units);

/**
 * @brief					Parses indirect inputs
 * @param[in,out]	ini		Dictionary to input file
 *
 * This is the function that handles input parameters with suffixes such as "pc"
 * and "tot". It's old and should probably be replaced by something better, but
 * it works.
 */
static void parseIndirectInput(dictionary *ini);

/**
 * @brief					Creates SI Units based on the SemiSI scheme
 * @param[in,out]	ini		Dictionary to input file
 * @return					Units
 *
 * Warning: This function is not quite clean in that it does actually write to
 * the ini-file. It would be preferrable if uAlloc() didn't have to do that.
 */
static Units *uSemiSI(dictionary *ini);

/**
 * @brief					Creates SI Units based on the SI scheme
 * @param			ini		Dictionary to input file
 * @return					Units
 */
static Units *uSI(const dictionary *ini);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

void uFree(Units *units){
	free(units->weights);
	free(units);
}

Units *uAlloc(dictionary *ini){

	parseIndirectInput(ini);

	char *method = iniGetStr(ini, "methods:normalization");

	Units *units = NULL;
	if(!strcmp(method, "semiSI"))	units = uSemiSI(ini);
	else if(!strcmp(method, "SI"))	units = uSI(ini);
	else msg(ERROR, "methods:normalization not valid (must be SI or semiSI)");

	free(method);

	uAddDerivedUnits(units);
	return units;
}

void uNormalize(dictionary *ini, const Units *units){

	/*
	 * Normalization of charge, mass and density, and superparticle weighting.
	 */
	int nSpecies = units->nSpecies;
	double *weights = units->weights;

	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);

	// Simulation particle scaling
	for(int s=0; s<nSpecies; s++){
		charge[s]  *= weights[s];
		mass[s]    *= weights[s];
		density[s] /= weights[s];
	}

	// Normalization
	adScale(charge, nSpecies, 1.0/units->charge);
	adScale(mass,   nSpecies, 1.0/units->mass);
	adScale(density, nSpecies, 1.0/units->density);

	iniSetDoubleArr(ini, "population:charge", charge, nSpecies);
	iniSetDoubleArr(ini, "population:mass", mass, nSpecies);
	iniSetDoubleArr(ini, "population:density", density, nSpecies);

	free(charge);
	free(mass);
	free(density);

	/*
	 * Normalization of everything else.
	 * This part can be used as a template for other normalize functions.
	 */
	iniScaleDouble(ini, "population:thermalVelocity", 1.0/units->velocity);
	iniScaleDouble(ini, "population:drift", 1.0/units->velocity);
	iniScaleDouble(ini, "population:perturbAmplitude", 1.0/units->length);
	iniScaleDouble(ini, "fields:BExt", 1.0/units->bField);
	iniScaleDouble(ini, "fields:EExt", 1.0/units->eField);
	//msg(ERROR,"1.0/units->eField = %f",1.0/units->eField);

}
/*
void uPrintSummary(const Units *units){
	double X = units->length;
	double T = units->time;
	double Q = units->charge;
	double M = units->mass;
	msg(STATUS, "Characteristic length: %15g  (%a)",X,X);
	msg(STATUS, "Characteristic time  : %15g  (%a)",T,T);
	msg(STATUS, "Characteristic charge: %15g  (%a)",Q,Q);
	msg(STATUS, "Characteristic mass  : %15g  (%a)",M,M);
}
*/

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

static void parseIndirectInput(dictionary *ini){

	/*
	 * APPLIES MULTIPLIERS TO ELEMENTS WITH SUFFICES
	 */

	int nDims = iniGetInt(ini, "grid:nDims");

	double V = (double)gGetGlobalVolume(ini);

	int *L = gGetGlobalSize(ini);
	double *mul = malloc(nDims*sizeof(nDims));
	for(int i=0;i<nDims;i++) mul[i] = 1.0/L[i];

	iniApplySuffix(ini, "population:nParticles", "pc", &V, 1);
	iniApplySuffix(ini, "population:nAlloc", "pc", &V, 1);
	iniApplySuffix(ini, "grid:nEmigrantsAlloc", "pc", &V, 1);
	iniApplySuffix(ini, "grid:stepSize", "tot", mul, nDims);

	free(mul);
}
static Units *uSemiSI(dictionary *ini){

	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);
	double timeStep = iniGetDouble(ini, "time:timeStep");

	const double tol = 1e-10;
	if(abs(charge[0]+1)>tol)
		msg(ERROR, "Species 0 must have charge -1 with this normalization");
	if(abs(mass[0]-1)>tol)
		msg(ERROR, "Species 0 must have mass 1 with this normalization");

	adScale(charge, nSpecies, elementaryCharge);
	adScale(mass, nSpecies, electronMass);

	double wpe = sqrt(pow(elementaryCharge,2)*density[0]/
			(vacuumPermittivity*electronMass));
	timeStep /= wpe;

	iniSetDoubleArr(ini, "population:charge", charge, nSpecies);
	iniSetDoubleArr(ini, "population:mass", mass, nSpecies);
	iniSetDouble(ini, "time:timeStep", timeStep);

	free(charge);
	free(mass);
	free(density);

	return uSI(ini);
}

static Units *uSI(const dictionary *ini){


	int nDims = iniGetInt(ini, "grid:nDims");
	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double timeStep = iniGetDouble(ini, "time:timeStep");
	double *stepSize = iniGetDoubleArr(ini, "grid:stepSize", nDims);
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);
	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);

	double V  = gGetGlobalVolume(ini)*pow(stepSize[0],nDims);

	double *weights = (double*)malloc(nSpecies*sizeof(*weights));
	for(int s=0; s<nSpecies; s++){
		weights[s] = density[s]*V/nParticles[s];
	}

	double X  = stepSize[0];
	double T  = timeStep;
	double Q  = weights[0]*fabs(charge[0]);
	double M  = pow(T*Q,2)/(vacuumPermittivity*pow(X,nDims));

	free(charge);
	free(mass);
	free(density);
	free(stepSize);
	free(nParticles);
	

	Units *units = malloc(sizeof(*units));
	units->nDims = nDims;
	units->nSpecies = nSpecies;
	units->weights = weights;
	units->length = X;
	units->time = T;
	units->charge = Q;
	units->mass = M;

	return units;
}

static void uAddDerivedUnits(Units *units){

	double nDims = units->nDims;
	double charge = units->charge;
	double mass = units->mass;
	double length = units->length;
	double time = units->time;

	units->hyperArea     = pow(length,nDims-1);
	units->hyperVolume   = pow(length,nDims);
	units->frequency     = 1.0/time;
	units->velocity      = length/time;
	units->acceleration  = length/pow(time,2);
	units->density       = 1.0/pow(length,nDims);
	units->chargeDensity = charge/pow(length,nDims);
	units->potential     = pow(length/time,2)*mass/charge;
	units->eField        = length*mass/(pow(time,2)*charge);
	units->bField        = mass/(time*charge);
	units->energy        = mass*pow(length/time,2);
}
