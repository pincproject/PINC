/**
 * @file		units.c
 * @brief		Handles units and normalization.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 *
 * Core module for handling of units and normalization.
 */

#include "core.h"
#include "units.h"

static void uAddDerivedUnits(Units *units);

/******************************************************************************
 * DEFINING NORMALIZATION FUNCTIONS
 *****************************************************************************/

void parseIndirectInput(dictionary *ini){

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
Units *uSemiSI(dictionary *ini){

	const double elementaryCharge = 1.60217733e-19; // [C]
	const double electronMass = 9.10938188e-31; // [kg]
	const double vacuumPermittivity = 8.854187817e-12; // [F/m]

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

void uFree(Units *units){
	free(units);
}

static void uAddDerivedUnits(Units *units){

	double nDims = units->nDims;
	double charge = units->charge;
	double mass = units->mass;
	double length = units->length;
	double time = units->time;

	units->velocity = length/time;
	units->acceleration = length/pow(time,2);
	units->chargeDensity = charge/pow(length,nDims);
	units->potential = pow(length/time,2)*mass/charge;
	units->eField = length*mass/(pow(time,2)*charge);
	units->bField = mass/(time*charge);
	units->energy = mass*pow(length/time,2);
}

Units *uSI(dictionary *ini){
	
	const double vacuumPermittivity = 8.854187817e-12; // [F/m]

	int nDims = iniGetInt(ini, "grid:nDims");
	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double timeStep = iniGetDouble(ini, "time:timeStep");
	double *stepSize = iniGetDoubleArr(ini, "grid:stepSize", nDims);
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);
	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);

	double V  = gGetGlobalVolume(ini)*pow(stepSize[0],nDims);

	double *K = (double*)malloc(nSpecies*sizeof(*K));
	for(int s=0; s<nSpecies; s++){
		K[s] = density[s]*V/nParticles[s];
	}

	double X  = stepSize[0];
	double T  = timeStep;
	double Q  = K[0]*fabs(charge[0]);
	double M  = pow(T*Q,2)/(vacuumPermittivity*pow(X,nDims));
	msg(STATUS, "Characteristic length: %15g m  (%a)",X,X);
	msg(STATUS, "Characteristic time  : %15g s  (%a)",T,T);
	msg(STATUS, "Characteristic charge: %15g C  (%a)",Q,Q);
	msg(STATUS, "Characteristic mass  : %15g kg (%a)",M,M);

	adMul(charge, K, charge, nSpecies);
	adMul(mass,   K, mass,   nSpecies);
	adScale(charge, nSpecies, 1.0/Q);
	adScale(mass,   nSpecies, 1.0/M);
	iniSetDoubleArr(ini, "population:charge", charge, nSpecies);
	iniSetDoubleArr(ini, "population:mass", mass, nSpecies);


	free(K);
	free(charge);
	free(mass);
	free(density);
	free(stepSize);
	free(nParticles);

	Units *units = malloc(sizeof(*units));
	units->length = X;
	units->time = T;
	units->charge = Q;
	units->mass = M;
	uAddDerivedUnits(units);
	uNormalize(ini, units);
	return units;
}

void uNormalize(dictionary *ini, const Units *units){

	iniScaleDouble(ini, "population:thermalVelocity", 1.0/units->velocity);
	iniScaleDouble(ini, "population:drift", 1.0/units->velocity);
	iniScaleDouble(ini, "population:perturbAmplitude", 1.0/units->length);
	iniScaleDouble(ini, "fields:BExt", 1.0/units->bField);
	iniScaleDouble(ini, "fields:EExt", 1.0/units->eField);

}


