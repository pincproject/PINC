/**
 * @file		collsisions.h
 * @brief		MCC Colisional module.
 * @author		Steffen Mattias Brask <steffemb@fys.uio.no>
 *
 * Functions to apply Monte Carlo Collisions between particles
 */

#ifndef COLLISIONS_H
#define COLLISIONS_H

/**
 * @brief Contains varibles used for MCC (collisions)
 */
typedef struct{
	double pMaxElectron;		/// < Max possible collsion probability
	double maxFreqElectron;	/// < Max possible collsion freq
	double pMaxIon;			/// < Max possible collsion probability
	double maxFreqIon;			/// < Max possible collsion freq
	double artificialLoss;		/// < Artificial energyloss electrons
	double nt;
	double NvelThermal;
	//crossect values
	double mccSigmaCEX;
	double mccSigmaIonElastic;
	double mccSigmaElectronElastic;
	double energyConvFactor; // convert from PINC to eV
	double massUnit; // convert from PINC to real


} MccVars;

/**
 * @brief Contains Method functions for MCC (collisions)
 */
typedef struct{
	///< Function pointer to ConstantCrossect coll method function
	void (*mccCollideConstantCrossect)(const dictionary *ini, Population *pop,
		MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);

} MccMethod;



/**
 * @brief Run designed for testing MCC
 * @param 	ini
 *
 *  Run designed for testing MCC
 *
 */

// void mccTestMode2(dictionary *ini);
// funPtr mccTestMode2_set(dictionary *ini);

void mccTestMode(dictionary *ini);
funPtr mccTestMode_set(dictionary *ini);

funPtr constCrossect_set(dictionary *ini);

funPtr collissionsOff_set(dictionary *ini);

/**
 * @brief brief description
 * @param one 		one description
 * @param two		two description
 *
 * Looong description of mccTest
 */


int mccTest(int one, int two);

MccVars *mccAlloc(const dictionary *ini, const Units *units);

void mccGetPmaxElectronConstantFrq(const dictionary *ini,
	MccVars *mccVars,Population *pop,MpiInfo *mpiInfo);

void mccGetPmaxIonConstantFrq(const dictionary *ini,MccVars *mccVars,
	Population *pop,MpiInfo *mpiInfo);

void mccGetPmaxIonFunctional(const dictionary *ini,
	MccVars *mccVars,Population *pop,
	MpiInfo *mpiInfo);

void mccGetPmaxElectronFunctional(const dictionary *ini,
	MccVars *mccVars, Population *pop,
	MpiInfo *mpiInfo);

void mccGetPmaxIonStatic(const dictionary *ini,MccVars *mccVars,
	Population *pop, MpiInfo *mpiInfo);

void mccGetPmaxElectronStatic(const dictionary *ini,
	MccVars *mccVars, Population *pop,
	MpiInfo *mpiInfo);

void mccCollideElectronStatic(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo);

void mccCollideIonStatic(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);

void mccCollideElectronFunctional(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo);

void mccCollideIonFunctional(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,MpiInfo *mpiInfo);

void mccCollideElectronConstantFrq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo);

void mccCollideIonConstantFrq(const dictionary *ini, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);


#endif // COLLISIONS_H
