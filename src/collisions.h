/**
 * @file		collsisions.h
 * @brief		MCC Colisional module.
 * @author		Steffen Mattias Brask <steffen.brask@fys.uio.no>
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
	double collFrqCex;
	double collFrqIonElastic;
	double collFrqElectronElastic;
	double CEX_a;
	double CEX_b;
	double ion_elastic_a;
	double ion_elastic_b;
	double electron_a;
	double electron_b;
	double energyConvFactor; // convert from PINC to eV
	double *neutralDrift;


} MccVars;


//void neutTest(dictionary *ini);
funPtr neutTest_set();

//void oCollMode(dictionary *ini);
funPtr oCollMode_set();

//void mccMode(dictionary *ini);
funPtr mccMode_set(dictionary *ini);


funPtr mccConstCrossect_set(dictionary *ini);
funPtr mccFunctionalCrossect_set(dictionary *ini);
funPtr mccConstFreq_set(dictionary *ini);
funPtr mccCollissionsOff_set(dictionary *ini);


/**
 * @brief allocates necesarry variables for the collision module.
 * @param ini 		pointer to input file as dictionary
 * @param units		pointer to units
 *
 * ini, and units must exist and be initialized before this function.
 */

MccVars *mccAlloc(const dictionary *ini, const Units *units);


//
// /**
//  * @brief updates Pmax for electrons in the constant collision freq. model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
//
// void mccGetPmaxElectronConstantFrq(const dictionary *ini,
// 	MccVars *mccVars,Population *pop,MpiInfo *mpiInfo);
//
//
// /**
//  * @brief updates Pmax for Ions in the constant collision freq. model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
// void mccGetPmaxIonConstantFrq(const dictionary *ini,MccVars *mccVars,
// 	Population *pop,MpiInfo *mpiInfo);
//
//
// /**
//  * @brief updates Pmax for Ions in the functional cross-sect collision model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
// void mccGetPmaxIonFunctional(const dictionary *ini,
// 	MccVars *mccVars,Population *pop,
// 	MpiInfo *mpiInfo);
//
// /**
//  * @brief updates Pmax for Electrons in the functional cross-sect collision model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
// void mccGetPmaxElectronFunctional(const dictionary *ini,
// 	MccVars *mccVars, Population *pop,
// 	MpiInfo *mpiInfo);
//
//
//
// /**
//  * @brief updates Pmax for Ions in the constant cross-sect collision model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  * @param *rng		pointer to random number generator
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
//  void mccGetPmaxIonStatic(const dictionary *ini,MccVars *mccVars,Grid *rhoNeutral,
//  	Population *pop, MpiInfo *mpiInfo,const gsl_rng *rng);
//
//
// /**
//  * @brief updates Pmax for Electrons in the constant cross-sect collision model
//  * @param *ini		pointer to input file as dictionary
//  * @param *mccVars		pointer to mcc specific variables
//  * @param *pop		pointer to particle population
//  * @param *mpiInfo		pointer to mpi specific variables
//  *
//  * mccAlloc() must be called prior to this function. Function is called in the
//  * main collission function.
//  */
// void mccGetPmaxElectronStatic(const dictionary *ini,
// 	MccVars *mccVars,Grid *rhoNeutral, Population *pop,
// 	MpiInfo *mpiInfo);
//

/**
 * @brief main collision function for electron constant cross-sect collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the constant cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideElectronStatic(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);


/**
 * @brief main collision function for Ion constant cross-sect collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the constant cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideIonStatic(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);


/**
 * @brief main collision function for Electrons functional cross-sect collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the functional cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideElectronFunctional(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo);


/**
 * @brief main collision function for Ion functional cross-sect collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the functional cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideIonFunctional(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,MpiInfo *mpiInfo);



/**
 * @brief main collision function for electron constant frequency collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the constant cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideElectronConstantFrq(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng,
	MpiInfo *mpiInfo);



/**
 * @brief main collision function for ion constant frequency collision model
 * @param *ini		pointer to input file as dictionary
 * @param *pop		pointer to particle population
 * @param *mccVars		pointer to mcc specific variables
 * @param *rng		pointer to random number generator
 * @param *mpiInfo		pointer to mpi specific variables
 *
 * This function is called in the method handler for the constant cross-sect
 * collision model. Initialization, i.e. Allocation must be performed prior to
 * this call. In practice only a call to collide() should be necessary in the
 * main time loop, and the rest is defined in the input file.
 */
void mccCollideIonConstantFrq(const dictionary *ini,Grid *rhoNeutral, Population *pop,
	MccVars *mccVars, const gsl_rng *rng, MpiInfo *mpiInfo);


#endif // COLLISIONS_H
