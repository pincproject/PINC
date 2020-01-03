/**
 * @file		population.h
 * @brief		Input/output functions.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */

#ifndef NEUTRALS_H
#define NEUTRALS_H



/**
 * @brief Contains a population of neutral charged particles.
 *
 * The position and velocity of particles is stored in a flat manner, such that
 * (x,y,z) of particle 0 comes first, then (x,y,z) of particle 1, and so on
 * (assuming 3D). As an example, printing the (x,y,z) position of the first 3
 * particles:
 *
 * @code
 *	NeutralPopulation neutralPop;
 *  ...
 *	for(int i=0;i<3;i++){
 *		double *iPos = &neutralPop.pos[i*neutralPop.nDims];
 *		printf("Particle %i is located at (%f,%f,%f).\n",i,iPos[0],iPos[1],iPos[2]);
 *	}
 * @endcode
 *
 *
 */
typedef struct{
	double *pos;		///< Position
	double *vel;		///< Velocity
	long int *iStart;	///< First index of specie s (nSpecies+1 elements)
	long int *iStop;	///< First index not of specie s (nSpecies elements)
	double *mass;		///< Mass (nSpecies elements)
	//double *kinEnergy;	///< Kinetic energy (nSpecies+1 elements)
	//double *potEnergy;	///< Potential energy (nSpecies+1 elements)
	int nSpeciesNeutral;		///< Number of species
	int nDims;			///< Number of dimensions (usually 3)
	bndType *bnd;		/// type of boundaries for particles, 2*nDims
	//hid_t h5;			///< HDF5 file handler
	double stiffnessConstant; //stiffnessConstant decides magnitude of pressure.
	double rho0;
} NeutralPopulation;


/**
 * @brief	Allocates memory for Neutral Population according to ini-file
 * @param	ini		Dictionary to input file
 * @return			NeutralPopulation
 *
 * Allocates memory for as many particles and species as specified in
 * collisions:nSpeciesNeutral and population:nAlloc in ini-file. This function only
 * allocates the memory for the particles, it does not generate them.
 *
 * Remember to call pNeutralFree() to free memory.
 */
NeutralPopulation *pNeutralAlloc(const dictionary *ini,const MpiInfo *mpiInfo);

/**
 * @brief					Frees memory for Neutral Population
 * @param[in,out]	pop		Pointer to population to be freed
 */
void pNeutralFree(NeutralPopulation *pop);



 //#########################################
 // Distributer
 // ########################################

funPtr NeutralDistr3D1_set(dictionary *ini);
void NeutralDistr3D1(const NeutralPopulation *pop, Grid *rho);

funPtr NeutralDistr3D1Vector_set(dictionary *ini);
void NeutralDistr3D1Vector(const NeutralPopulation *pop, Grid *rho);


//######################################


void neFillGhost(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo);

void nePurgeGhost(NeutralPopulation *pop, const Grid *grid);

void neVelMaxwell(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng);

void nePosLattice(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo);

void nePosUniform(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo, const gsl_rng *rng);

void neVelAssertMax(const NeutralPopulation *pop, double max);

void neInjectParticles(int slicePos,int dim,int multiplyDens,const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo);
//#########################################
// Mover/Accelerator
// ########################################

void neMove(NeutralPopulation *pop);
void neAcc3D1(NeutralPopulation *pop, Grid *Pgrad, Grid *divBulkV);
funPtr neAcc3D1_set(dictionary *ini);

//#########################################
// Finite diff functions
// ########################################

void divFinDiff1st(const Grid *scalar, Grid *field);


//#############################
// Migration
//#############################
void neExtractEmigrants3DOpen(NeutralPopulation *pop, MpiInfo *mpiInfo);
funPtr neExtractEmigrants3DOpen_set(const dictionary *ini);
void neMigrate(NeutralPopulation *pop, MpiInfo *mpiInfo, Grid *grid);
void nePosAssertInLocalFrame(const NeutralPopulation *pop, const Grid *grid);

//#############################
// Pressure solver
//#############################

void nePressureSolve3D(Grid *rhoNeutral,Grid *P,NeutralPopulation *pop, const MpiInfo *mpiInfo);

void neSetBndSlices(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo);

void neSetBndSlicesRho(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo);


//#############################
// Object functions
//#############################


void nuObjectCollide(NeutralPopulation *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo);
void nuObjectpurge(NeutralPopulation *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo);
#endif // POPULATION_H
