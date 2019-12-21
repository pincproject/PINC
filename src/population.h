/**
 * @file		population.h
 * @brief		Input/output functions.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */

#ifndef POPULATION_H
#define POPULATION_H

/**
 * @brief	Allocates memory for Population according to ini-file
 * @param	ini		Dictionary to input file
 * @return			Population
 *
 * Allocates memory for as many particles and species as specified in
 * populations:nSpecies and population:nAlloc in ini-file. This function only
 * allocates the memory for the particles, it does not generate them.
 *
 * Remember to call pFree() to free memory.
 */
Population *pAlloc(const dictionary *ini,const MpiInfo *mpiInfo);

/**
 * @brief					Frees memory for Population
 * @param[in,out]	pop		Pointer to population to be freed
 */
void pFree(Population *pop);

/**
 * @brief	Assign particles uniformly distributed positions
 * @param			ini		Dictionary to input file
 * @param[in,out]	pop		Population of particles
 * @param			grid	Grid in which the particles are to be distributed
 * @param			rngSync	Synchronized random number generator
 * @return			void
 *
 * The amount of particles specified by population:nParticles in ini will be
 * generated with uniformly distributed random positions within the simulation
 * domain (global reference frame). In case of multiple subdomains only
 * particles residing in this MPI node's subdomain will be stored, and will be
 * transformed to its local reference frame. The rng should have the same seed
 * (be synchronized) on all MPI nodes when calling this function as that will
 * ensure that all nodes generates the same particles and discards particles
 * not belonging to their subdomain. Failure to do so may lead to the number of
 * particles generated being different than specified in ini.
 *
 * Beware that this function do not assign any velocity to the particles.
 * @see pVelMaxwell()
 */
void pPosUniform(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const gsl_rng *rngSync);

void pPosLattice(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo);

/**
 * @brief	Assign particles artificial positions suitable for debugging
 * @param			ini		Dictionary to input file
 * @param[in,out]	pop		Population
 *
 * The amount of particles specified by population:nParticles in ini will be
 * generated with values given by the following code:
 *
 * @code
 *	pos[i*nDims+d] = 1000*mpiRank + i + (double)d/10 + (double)s/100;
 * @endcode
 */
void pPosDebug(const dictionary *ini, Population *pop);


void pPosPerturb(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo);

/**
 * @brief Set the same velocity to all particles
 * @param[in,out]	pop		Population
 * @param			vel		Velocity to set (expected to be pop->nDims long)
 */
void pVelSet(Population *pop, const double *vel);

void pVelZero(Population *pop);
void pVelConstant(const dictionary *ini, Population *pop, double constant1, \
   double constant2);

void pPosAssertInLocalFrame(const Population *pop, const Grid *grid);
void pVelAssertMax(const Population *pop, double max);
void pSumPotEnergy(Population *pop);
void pSumKinEnergy(Population *pop);

/**
 * @brief	Assign particles Maxwellian distributed velocities
 * @param			ini		Dictionary to input file
 * @param[in,out]	pop		Population of particles
 * @param			rng		Random number generator
 * @return			void
 *
 * Iterates through all particles belonging to pop and assignes Maxwellian
 * distributed velocities to them, according to the temperature specified in
 * ini. Contrary to in posUniform() rng should at this point _not_ have the same
 * seed as that will lead to particles in different sub-domains having identical
 * velocities.
 */
void pVelMaxwell(const dictionary *ini, Population *pop, const gsl_rng *rng);


/**
 * @brief	Add new particles on boundary edge according to drift.
 * @param			ini		Dictionary to input file
 * @param[in,out]	pop		Population of particles
 * @param			rng		Random number generator
 * @param			mpiInfo MpiInfo
 * @return			void
 *
 * TODO:
 */
void pInfluxDrift(const dictionary *ini, Population *pop, const gsl_rng *rng, const MpiInfo *mpiInfo);
void pPurgeGhost(Population *pop, const Grid *grid);
void pFillGhost(const dictionary *ini, Population *pop, const gsl_rng *rng, const MpiInfo *mpiInfo);


/**
 * @brief	Add new particle to population
 * @param[in,out]	pop		Population
 * @param			s		Specie of new particle
 * @param			pos		Position of new particle (nDims elements)
 * @param			vel		Velocity of new particle (nDims elements)
 * @return			void
 */
void pNew(Population *pop, int s, const double *pos, const double *vel);

/**
 * @brief	Cut a particle from a population
 * @param[in,out]	pop		Population
 * @param			s		Which specie the particle belongs to
 * @param			p		Which index p the particle starts at
 * @param[out]		pos		The particle's position
 * @param[out]		vel		The particle's velocity
 * @return					void
 *
 * Note that the particle to fetch is adressed by the array index p. Thus,
 * if the third particle in the population is to be fetched p=2*nDims. In
 * addition, the specie it belongs to, s, must be specified. This rather odd
 * way of specifying which particle to fetch is because p and s of the particle
 * in quest is often already known and calculating p and s from for instance the
 * particle number i is then unnecessary amount of operations. Failure to
 * provide valid values of p and s results in unpredictable behaviour, with
 * the likely consequence of corrupting the whole population.
 */
void pCut(Population *pop, int s, long int p, double *pos, double *vel);

/**
 * @brief	Creates .pop.h5-file to store population in
 * @param	ini				Dictionary to input file
 * @param	pop[in,out]		Population
 * @param	fName			Filename
 * @return	void
 * @see pWriteH5(), pCloseH5()
 *
 * An output file is created whose filename is as explained in openH5File().
 * Remember to call pCloseH5().
 *
 * The file will have one group "/pos" for position data and one group "/vel"
 * for velocity data. Each of these will have groups "specie <s>" for each
 * specie. For each time-step, the population data will be stored in a dataset
 * named "n=<timestep>" where <timestep> is signified with one decimal allowing
 * interleaved quantities.
 *
 * The population data is stored in normalized values. To obtain physical
 * values, the output contains two attributes: The "Position denormalization
 * factor" and "Velocity denormalization factor". Upon multiplication by these,
 * the position/velocity will become physical values.
 */
void pOpenH5(	const dictionary *ini, Population *pop, const Units *units,
	   			const char *fName);

/**
 * @brief	Stores particles in Population in .pop.h5-file
 * @param	pop		Population
 * @param	mpiInfo	MpiInfo
 * @param	posN	Timestep of position data to be stored
 * @param	velN	Timestep of velocity data to be stored
 * @return			void
 *
 * The position and velocity of all particles are stored, referred to global
 * reference frame. The function takes care of merging the particles from all
 * MPI nodes to one file.
 *
 * NB: pop is not constified because all particles are transformed to global
 * reference frame before writing to .h5-file. However, they are transferred
 * back to local reference frame after writing so pop should remain unchanged to
 * within machine precision.
 */
void pWriteH5(Population *pop, const MpiInfo *mpiInfo, double posN, double velN);

/**
 * @brief	Closes .pop.h5-file
 * @param	pop		Population
 * @return	void
 */
void pCloseH5(Population *pop);

/**
 * @brief Transforms particles to local reference frame
 * @param	pop			Population of particles
 * @param	mpiInfo		MPI information about the reference frames
 * @return	void
 * @see toGlobalFrame()
 */
void pToLocalFrame(Population *pop, const MpiInfo *mpiInfo);

/**
 * @brief Transforms particles to global reference frame
 * @param	pop			Population of particles
 * @param	mpiInfo		MPI information about the reference frames
 * @return	void
 * @see toLocalFrame()
 *
 * For parallelization by means of configuration space decomposition, the
 * the particles' positions are usually specified with respect to a local
 * reference frame to that subdomain in order to ease computation. Some
 * operations may require the position in global reference frame (e.g. when
 * storing to file) for which purpose it can be transformed using this function.
 */
void pToGlobalFrame(Population *pop, const MpiInfo *mpiInfo);

/**
 * @brief Creates datasets in .xy.h5-file for storing energy
 * @param	xy		.xy.h5-identifier
 * @param	pop		Population
 * @return			void
 *
 * Uses xyCreateDataset() to create datasets corresponding to the potential and
 * kinetic energies stored in Population. Datasets for both kinetic and
 * potential energy is created for each specie separately, as well as total for
 * for all species, e.g. for two species:
 *	- /energy/kinetic/specie 0
 *  - /energy/kinetic/specie 1
 *	- /energy/kinetic/total
 *	- /energy/potential/specie 0
 *	- /energy/potential/specie 1
 *	- /energy/potential/total
 *
 * pWriteEnergy() can be used to populate these datasets.
 */
void pCreateEnergyDatasets(hid_t xy, Population *pop);

void pCreateTemperatureDatasets(hid_t xy, Population *pop);

/**
 * @brief Writes energies to .xy.h5-file
 * @param	xy		.xy.h5-identifier
 * @param	pop		Population
 * @return			void
 *
 * Uses xyWrite() to write potential and kinetic energies stored in Population
 * to .xy.h5 datasets. These datasets must first be created using
 * pCreateEnergyDatasets(). Note that this function does not populate the energy
 * variables in Population with meaningful values, energy-computing functions
 * such as gPotEnergy() and puAcc3D1KE() must be used for that. Beware though
 * that energy computing functions typically only populate the energy per
 * specie, or if that is unobtainable by the algorithm, the summed (total)
 * energy for all species. In the former case, the total energy can be obtained
 * simply by addition during post-processing.
 */
void pWriteEnergy(hid_t xy, Population *pop, double x,Units *units);

/**
 * Same as pWriteEnergy except temperature instead of kinetic energy.
 */
void pWriteTemperature(hid_t xy, Population *pop, double x,Units *units,dictionary *ini);

/*
Finding particles in population close to object, discards
particles that will not intersect object next timestep
pop->vicinity contains index of particles that are close
to the object
*/
void pVicinityParticles(Population *pop, const Object *obj);

/*
Integral of radiance for specific wavenumber to infinity, finds
total energy in Joule in the band per timetep
*/


/*
Integral of radiance for specific wavenumber to infinity, finds
number of photons that cause electrons to be emitted per timestep
*/
double *pPlanckPhotonIntegral(dictionary *ini, const Units *units, Object *obj);

/*
emits photoELectrons from object surface every timestep
*/
void pPhotoElectrons(dictionary *ini, Population *pop, const Object *obj, 
                     const Units *units, const gsl_rng *rng, const MpiInfo *MpiInfo);

/*
Find Collision type is called during puMove; returns
pSecondaryElectron/pReflect/pAdhere
*/
funPtr pFindCollisionType(Population *pop, Object *obj, long int n);
void pSecondaryElectron(Population *pop, const Object *obj, long int n);
void pReflect(Population *pop, const Object *obj, long int n, const MpiInfo *mpiInfo);
void pAdhere(Population *pop, const Object *obj, long int n);

#endif // POPULATION_H
