#ifndef PUSHER_H
#define PUSHER_H

/**
 * @file		pusher.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle Pusher
 * @date		16.12.15
 *
 * All parts of a standard particle pusher is implemented herein, including
 * particle mover, particle accelerators (leapfrog, boris) and interpolation
 * schemes of various orders.
 */

/******************************************************************************
 * DECLARING DATATYPES
 *****************************************************************************/
/**
 * @brief Contains a population of particles.
 *
 *
 */
//  typedef struct{
// 	long int *migrants;		///< Index p of migrants (size specified in input file)
// 	int nNeigh;				///< =pow(3,nDims)-1
// 	int here;				///< Number of this domain itself
// 	long int *neighStart;	///< At which index in 'migrants' the neighbours start (nNeigh elements)
// 	long int *nMigrants;	///< Number of migrants of each specie for each neighbour (nSpecies*nNeigh elements)
// 	double *bufferPos;		///< Buffer for receiving position
// 	double *bufferVel;		///< Buffer for receiving velocity
// 	double limits;			///< limits of migration [xl, yl, zl, xu, yu, zu]
// } Migrants;


/**
 * @brief Moves particles one time-step forward
 * @param[in,out]	pop		Population
 * @return					void
 *
 * No boundary conditions are enforced and particles may therefore travel out of
 * bounds. Other functions must be called subsequently to enforce boundary
 * conditions or transfer them to other sub-domains as appropriate. Otherwise
 * PINC may fail ungracefully.
 */
void puMove(Population *pop);

/**
 * @brief Enforce periodic particle boundary conditions
 * @param[in,out]	pop		Population
 * @param			grid	Grid
 *
 * Run after puMove() to enforce conditions. Only the grid size is used from the
 * grid input so it may be any quantity with the correct dimensions. Typically
 * the E-field.
 */
void puBndPeriodic(Population *pop, const Grid *grid);

/**
 * @brief Distributes charge density on grid using 1st order interpolation. Fixed to 3D.
 * @param			pop		Population
 * @param[in,out]	rho		Charge density
 * @return					void
 *
 * Assuming particles are correctly placed prior to calling this function.
 * puMove() and some boundary enforcing function should therefore be called
 * first.
 */
void puDistr3D1(const Population *pop, Grid *rho);

/**
 * @brief Accelerates particles using 1st order interpolation. Fixed to 3D.
 * @param[in,out]	pop		Population
 * @param			E		E-field to determine acceleration from
 * @return					void
 *
 * Accelerates particles by interpolating the E-field and updating the velocity
 * of the particles by one time-step. A half time-step may be achieved by
 * multiplying the E-field by 0.5 prior to acceleration using gMulDouble().
 * This is useful to initialize a leapfrog algorithm. Remember to multiply
 * E-field by 2 again afterwards.
 *
 * E is not constified since it will be re-normalized (scaled) for each specie
 * in order to speed up calculations. However, when this function is done it is
 * scaled back to original again.
 *
 * NB: Only works on 3D
 */
void puAcc3D1(Population *pop, Grid *E);
void puAcc3D1KE(Population *pop, Grid *E);

void puIdMigrants3D(Population *pop, MpiInfo *mpiInfo);
void puIdMigrantsND(Population *pop, MpiInfo *mpiInfo);

void puExtractEmigrantsND(Population *pop, MpiInfo *mpiInfo);
void puExtractEmigrants3D(Population *pop, MpiInfo *mpiInfo);

void puMigrate(Population *pop, MpiInfo *mpiInfo, Grid *grid);

int puRankToNeighbor(MpiInfo *mpiInfo, int rank);
int puNeighborToRank(MpiInfo *mpiInfo, int neighbor);
int puNeighborToReciprocal(int neighbor, int nDims);
void puBndIdMigrants3D(Population *pop, MpiInfo *mpiInfo);
void puBndIdMigrantsND(Population *pop, MpiInfo *mpiInfo);

#endif // PUSHER_H
