/**
 * @file		pusher.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle Pusher
 * @date		16.12.15
 *
 * All parts of a standard explicit particle pusher is implemented herein,
 * including particle mover, particle accelerators (leapfrog, boris) and
 * interpolation schemes of various orders.
 */

#ifndef PUSHER_H
#define PUSHER_H

/******************************************************************************
 * DECLARING DATATYPES
 *****************************************************************************/

/**
 * @brief Moves particles one timestep forward
 * @param[in,out]	pop		Population
 * @return					void
 *
 * No boundary conditions are enforced and particles may therefore travel out of
 * bounds. Other functions must be called subsequently to enforce boundary
 * conditions or transfer them to other sub-domains as appropriate. Otherwise
 * PINC may fail ungracefully.
 */
void puMove(Population *pop);

/** @name Accelerators (with interpolation)
 * This group of functions accelerates the particles in a population by
 * advancing their velocities one time step. The interpolation of field
 * quantities from the grid to the particles is embedded into the accelerators
 * for better performence, despite a loss of modularity\footnote{Having
 * separate accelerators and interpolators require either (a) that the
 * interpolated field values of all particles be stored until subsequently
 * calling a function for accelerating the particles which is very memory
 * demanding, or (b) having a loop iterating through the particles with one
 * interpolator and one accelerator function, yielding two function calls per
 * particle. While function calls are normally not considered that costly doing
 * it per particle is likely disastrous. The operations per particle should
 * largely be limited to arithmetic operations requiring only minimal memory
 * look-up. Internally, though, the interpolators *are* in fact separate inline
 * functions acting per particle giving some feeling of modularity. Since they
 * are inlined no function calls should actually be made.}. The following
 * naming convention exist, according to which algorithms they use:
 *
 *	name				| description
 *	--------------------|-----------------------------------------------------------------------------------
 *	puAccXDY()			| Simple increment in velocity (e.g. leapfrog) (no external B-field)
 *	puAccXDYKE()		| Same as above but computes kinetic energy for each specie at the mid-step
 *	puBorisXDY()		| Boris algorithm (for external homogeneous B-field, S and T are vectors)
 *	puBorisXDYKE()		| Same as above but computes kinetic energy for each specie at the mid-step
 *	puBorisInhXDY()		| Boris algorithm (for external inhomogeneous B-field, S and T are Grid quantities)
 *	puBorisInhXDYKE()	| Same as above but computes kinetic energy for each specie at the mid-step
 *
 * Where X indicates the dimensionality and Y the order of interpolation used
 * in wheighting the field(s) from the grid nodes to the particles, e.g. Y=0
 * for the NGP method and Y=1 for the PIC/CIC method. The
 * interpolation is carried out by the underlying functions named
 * puInterpXDY() according to the same convention. Functions with X=N works on
 * configurations of arbitrary dimensionality, which is commonplace in PINC.
 * However, since N-dimensional interpolation is significantly more
 * time-consuming than algorithms with fixed dimensionality (at least for order
 * higher than 0) some fixed dimensionality algorithms are included. For
 * instance, puInterp3D1() is much faster than puInterpND1().
 *
 * Remember that Boris and leapfrog methods require the velocities to be
 * located at half-integer steps. This initialization of the velocities can be
 * performed by multiplying E (and S and T in case of Boris) by 0.5,
 * accelerating once, and restoring E (and S and T) by multiplying them by 2.
 * For instance to get a leapfrog iteration:
 *
 * @code
 *	// Assume position and velocity initialized at timestep 0 here
 *
 *	gMul(E,0.5);
 *	puAcc3D1KE(pop,E); // Increment velocity to timestep 0.5
 *	gMul(E,2);
 *
 *	for(int n=1; n<=nTimeSteps; n++){ // Mind the range of n
 *
 *		puMove(pop);  // Advance position to timestep n
 *
 *		// Enforce boundaries here (see migrate.h)
 *
 *		puDistr3D1(pop, rho); // Distribute charges onto grid
 *
 *		// Solve field here (see e.g. multigrid.h)
 *
 *		puAcc3D1KE(pop, E); // Advance velocity to timestep n+0.5
 *	}
 * @endcode
 *
 * The kinetic energy-computing accelerators will here compute the energy at
 * integer steps (n) leaving the total kinetic energy per specie, per
 * subdomain, in the variable pop. Summing across the subdomains and storing to
 * file can be done by pWriteEnergy().
 *
 * @param[in,out]	pop		Population
 * @param			E		Electric field
 * @param			S		Rotation parameter (Boris only)
 * @param			T		Rotation parameter (Boris only)
 * @return					void
 *
 * The E input is usually not constified since it is rescaled several times
 * inside the function (specie-specific renormalization). By the end of the
 * function call, however, it should be restored to its initial value (to within
 * machine precision).
 *
 * The rotation parameters S and T for the homogeneous Boris methods are
 * generated from the external B-field before the loop by
 * puGet3DRotationParameters(). For inhomogeneous fields S and T are Grid
 * quantities (no function to create them yet). For slowly time-varying
 * magnetic fields S and T can be regenerated each iteration. However, using a
 * Poisson solver does not properly deal with electromagnetic effects, so if
 * considering this the magnetic field must be kept quasi-static (and
 * quasi-homogeneous?).
 */
///@{
void puAcc3D1(Population *pop, Grid *E);
void puAcc3D1KE(Population *pop, Grid *E);
void puAccND1(Population *pop, Grid *E);
void puAccND1KE(Population *pop, Grid *E);
void puAccND0(Population *pop, Grid *E);
void puAccND0KE(Population *pop, Grid *E);
void puBoris3D1(Population *pop, Grid *E, const double *T, const double *S);
void puBoris3D1KE(Population *pop, Grid *E, const double *T, const double *S);

funPtr puAcc3D1_set(dictionary *ini);
funPtr puAcc3D1KE_set(dictionary *ini);
funPtr puAccND1_set(dictionary *ini);
funPtr puAccND1KE_set(dictionary *ini);
funPtr puAccND0_set(dictionary *ini);
funPtr puAccND0KE_set(dictionary *ini);
///@}

/**
 * @brief Generates rotation parameters for puBoris-funcitons.
 * @param			ini		Input file
 * @param[out]		T		Rotation parameter named t in B&L
 * @param[out]		S		Rotation parameter named s in B&L
 *
 * S and T must be pre-allocated to hold 3*nSpecies doubles each.
 * This functions needs some cleanup.
 */
void puGet3DRotationParameters(dictionary *ini, double *T, double *S);


/** @name Accelerators (with interpolation)
 * These functions distributes or deposits charges onto the charge densty grid
 * by interpolating the charges onto the nearest grid points. They are named
 * puDistrXDY() where X signifies the dimensionality and Y the order of
 * interpolation (similart to puInterpXDY() and the accelerator functions).
 *
 * These functions will crash ungracefully if particles are placed
 * out-of-bounds or out-of-threshold area. Make sure to migrate particles to
 * other subdomains before calling.
 *
 * @param			pop		Population
 * @param[in,out]	rho		Charge density
 * @return					void
 */
///@{
void puDistr3D1(const Population *pop, Grid *rho);
void puDistrND1(const Population *pop, Grid *rho);
void puDistrND0(const Population *pop, Grid *rho);

funPtr puDistr3D1_set(dictionary *ini);
funPtr puDistrND1_set(dictionary *ini);
funPtr puDistrND0_set(dictionary *ini);
///@}

// EVERYTHING BELOW THIS SHOULD MOVE TO SEPARATE MIGRATION.H MODULE.

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
