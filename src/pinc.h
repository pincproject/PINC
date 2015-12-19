/**
 * @file		pinc.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 				Gullik Vetvik Killie <gullikvk@gys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC main header.
 * @date		11.10.15
 *
 * The PINC main header file constitutes a framework of function and struct
 * declarations central to all PINC modules.
 */

#ifndef PINC_H
#define PINC_H

#include "iniparser.h"
#include <time.h>
#include <hdf5.h>
#include <gsl/gsl_rng.h>

/******************************************************************************
 * DEFINING PHYSICAL CONSTANTS (following SI standard)
 *****************************************************************************/

#define BOLTZMANN 1.3806488e-23
#define ELECTRON_MASS 9.10938215e-31

/******************************************************************************
 * DECLARING DATATYPES
 *****************************************************************************/
/**
 * @brief Contains a population of particles.
 *
 * The position and velocity of particles is stored in a flat manner, such that
 * (x,y,z) of particle 0 comes first, then (x,y,z) of particle 1, and so on
 * (assuming 3D). As an example, printing the (x,y,z) position of the first 3
 * particles:
 *
 * @code
 *	Population pop;
 *  ...
 *	for(int i=0;i<3;i++){
 *		double *iPos = &pop.pos[i*pop.nDims];
 *		printf("Particle %i is located at (%f,%f,%f).\n",i,iPos[0],iPos[1],iPos[2]);
 *	}
 * @endcode
 *
 * The subset of particles belonging to specie s start at particle i=iStart[s]
 * and stop at i=iStop[s]-1 (the -1 simplifies many calculations and is more
 * consistent with the C way of counting). Due to particle increase/decrease
 * there is allocated space for more particles than are present. Thus the
 * _allocated_ space for specie s start at i=iStart[s] and stop at
 * i=iStart[s+1]-1. For convenience, iStart has nSpecies+1 elements such that
 * this is true also for the last specie. The last element is then simply the
 * number of particles allocated in total.
 *
 * The position of the particles is normalized with respect to 'stepSize' in
 * Grid, such that a particle with local position (1,2,3) is located _on_ node
 * (1,2,3) in the grid. Particles are usually specified in local frame but may
 * temporarily be expressed in global frame. See MpiInfo.
 *
 * energy[i] simply holds the kinetic energy of particle i, or if energy is not
 * computed, energy points to NULL.
 *
 * If a population h5 output file is created, the handler to this file is
 * stored in h5.
 */
typedef struct{
	double *pos;		///< Position
	double *vel;		///< Velocity
	double *energy;		///< Kinetic energy
	long int *iStart;	///< First index of specie s (nSpecies+1 elements)
	long int *iStop;	///< First index not of specie s (nSpecies elements)
	double *renormRho;	///< Re-normalization factors for rho (nSpecies elements)
	double *renormE;	///< Re-normalization factors for E (nSpecies elements)
	int nSpecies;		///< Number of species
	int nDims;			///< Number of dimensions (usually 3)
	hid_t h5;			///< HDF5 file handler
} Population;

/**
 * @brief Contains information regarding how the PIC code is parallelized.
 *
 * The total simulation domain can be split across several subdomains, one for
 * each MPI node, where position (J,K,L) (in case of 3D) of each subdomain is
 * stored in 'subomain'. 'nSubdomains' represents the number of subdomains along
 * each dimension, and nSubdomainsProd is the cumulative product of nSubdomains
 * similarly as sizeProd in Grid.
 *
 * The local reference frame is defined such that a particle with integer
 * position (j,k,l) would be located _on_ grid point (j,k,l). This makes it fast
 * to determine the index of the nodes surrounding a particle. See Population.
 *
 * offset is the offset of this MPI node's subdomain with respect to the global
 * reference frame. Adding/subtracting this to a position converts to/from the
 * global reference frame. See toLocalFrame(), toGlobalFrame().
 *
 * posToSubdomain is a factor which can be used to determine which subdomain a
 * globally specified position belongs to, e.g. for 1D:
 *
 * @code
 *	int J = (int)(posToNode[0]*pos[0]);
 * @endcode
 */
typedef struct{
	int mpiRank;				///< MPI rank
	int mpiSize;				///< MPI size
	int nDims;					///< Number of dimensions
	int *subdomain;				///< MPI node (nDims elements)
	int *nSubdomains;			///< Number of MPI nodes (nDims elements)
	int *nSubdomainsProd;		///< Cumulative product of nSubdomains (nDims+1 elements)
	int *offset;				///< Offset from global reference frame (nDims elements)
	double *posToSubdomain;		///< Factor for converting position to subdomain (nDims elements)
} MpiInfo;

/**
 * @brief A grid-valued quantity, for instance charge density or E-field.
 *
 * This datatype can represent both scalar fields and vector fields on an N-
 * dimensional grid. Basically, it can be thought of as an arbitrarily
 * dimensioned array object which is stored flat/linearly in memory but which
 * contains supporting variables in order to work with it efficiently.
 *
 * The values in the array is stored in natural/lexicographical ordering in
 * 'val', such that for a 3D array, the elements has the following order:
 *
 *	(0,0,0), (1,0,0), (2,0,0), ..., (0,1,0), (1,1,0), (2,1,0), ...
 *
 * 'rank' is the number of dimensions of the array, while 'size' is the size
 * along each dimension. A 128x128x128 array therefore has rank=3 and
 * size={128,128,128}. The following is an example of how to access element
 * (a,b,c) in the array:
 *
 * @code
 *	Grid *rho;
 *	...
 *	int *size = rho->size;
 *	long int p = a + b*size[0] + c*size[0]*size[1];
 *	double element = rho->val[p];
 * @endcode
 *
 * For convenience and to speed up computations, 'sizeProd' is the cumulative
 * product of 'size' starting at 1. For our example,
 * sizeProd={1,128,128*128,128*128*128}. Using this, the linear index p in the
 * above code can equivalently be computed as:
 *
 * @code
 *	long int p = a*sizeProd[0] + b*sizeProd[1] + c*sizeProd[2];
 * @endcode
 *
 * Also note that adding sizeProd[d] to a linear index p corresponds to to
 * incrementing one step along dimension d in the array. This can be utilized
 * for speeding up certain calculations. sizeProd[rank] is the total number of
 * elements in the array.
 *
 * While this struct is generic and may store any kind of array, the standard
 * way of storing quantities on a grid in PINC is that the first dimension in
 * the array represents the field component (e.g. x, y or z-component of a 3D
 * electric field) while the consecutive dimensions represents the physical
 * dimensions in the grid. For instance a vector field on a 3D grid of
 * 128x128x128 grid points would be represented by a 3x128x128x128 array of
 * rank 4. This means that all field components are lumped together, stored one
 * grid point at a time. A scalar field on the same grid would be represented by
 * a 1x128x128x128 array also of rank 4. The rank is therefore one more than the
 * number of dimensions of the grid. sizeProd[2], for instance, then represents
 * an increment in y-direction in the grid. The following more involved example
 * prints the 3D vector value E on grid points (2,k,3) for all k avoiding
 * complete re-computation of the linear index p each time:
 *
 * @code
 *	Grid *E;
 *
 *	...
 *
 *	int *size = E->size;
 *	long int *sizeProd = E->sizeProd;
 *
 *	int j = 2;
 *	int l = 3;
 *	long int p = j*sizeProd[1] + l*sizeProd[3];
 *
 *	for(int k=0;k<size[2];k++){
 *
 *		double *value = E->val[p];
 *		printf("E(%i,%i,%i)=(%f,%f,%f)\n",j,k,l,value[0],value[1],value[2]);
 *		p += sizeProd[2];
 *	}
 * @endcode
 *
 * In the case of paralellization through domain decomposition this datatype
 * represents the local sub-domain and an MpiInfo object keeps track of how the
 * sub-domaines are related to one another. In the case of domain decomposition,
 * however, each MPI node will need data on grid points outside of its own
 * sub-domain. Thus the outermost layers of grid points will often be ghost
 * points of neighbouring sub-domains. Also, one might need ghost points for
 * some boundary condition implementations regardless of domain decomposition.
 * In any case, 'nGhostLayers' represents the number of ghost layers at each
 * boundary. For 3D (ND), the first 3 (N) elements of 'nGhostLayers' indicates
 * the number of ghost points along the lower boundaries of dimensions x, y, and
 * z, respectively. The last 3 (N) elements indicates the number of ghost points
 * along the upper boundaries. 'trueSize' indicates how many grid points truly
 * belongs to this sub-domain.
 *
 * Consider for instance a vector field on a sub-domain of 128x128x128 grid
 * padded with 1 layer of ghost nodes along the whole sub-domain. This will lead
 * to the struct having the following member values:
 *	- rank = 4
 *	- size = {3,130,130,130}
 *	- trueSize = {3,128,128,128}
 *	- sizeProd = {1,3,390,50700,6591000}
 *	- nGhostLayers = {0,1,1,1,0,1,1,1}
 *
 * Naturally, domain decomposition doesn't affect the first (non-physical)
 * dimension in the array.
 *
 * 'stepSize' is the step-size of each dimension in terms of Debye lengths (or
 * possibly some other quantity in the future). This doesn't make sense for the
 * first non-physical dimension and is therefore arbitrarily set to 1. Thence
 * the product of all elements in 'stepSize' is the volume of a cell.
 *
 * 'h5' is a HDF5 file identifier used to store the grid quantity to an .h5-file
 * and are used by gWriteH5(). The other two h5-variables are also used by
 * gWriteH5() since they only needs to be computed once.
 *
 * 'slice' is a buffer which is large enough to store any slice cut through the
 * array using getSlice().
 */

typedef struct{
	double *val;				///< Array of values on the grid
	int rank;					///< Number of dimensions of array (not grid)
	int *size;					///< Size of array (including ghosts) (rank elements)
	int *trueSize;				///< Size of array (excluding ghosts) (rank elements)
	long int *sizeProd;			///< Cumulative product of size (rank+1 elements)
	int *nGhostLayers;			///< Number of ghost layers in grid (2*rank elements)
	double *stepSize;			///< Step-sizes in Debye lengths (rank elements)

	double *slice;				///< Slice buffer of the grid sent to other subdomains
	hid_t h5;					///< HDF5 file handler
	hid_t h5MemSpace;			///< HDF5 memory space description
	hid_t h5FileSpace;			///< HDF5 file space description
} Grid;

typedef struct timespec TimeSpec;

/**
 * @brief	Timer struct for simple profiling
 * @see allocTimer(), freeTimer(), tMsg()
 */
typedef struct{
	TimeSpec previous;			///< Time of previous call
	int rank;					///< Rank of node or negative to turn off timer
} Timer;

/******************************************************************************
 * DEFINED IN POPULATION.C
 *****************************************************************************/

/**
 * @brief	Allocates memory for Population according to ini-file
 * @param	ini		Dictionary to input file
 * @see freePopulation(), posUniform(), velMaxwell()
 *
 * Allocates memory for as many particles and species as specified in
 * populations:nSpecies and population:nAlloc in ini-file. This function only
 * allocates the memory for the particles, it does not generate them.
 *
 * Remember to call freePopulation() to free memory.
 */
Population *pAlloc(const dictionary *ini);

/**
 * @brief					Frees memory for Population
 * @param[in,out]	pop		Pointer to population to be freed
 * @see allocPopulation()
 */
void pFree(Population *pop);

/**
 * @brief	Sets specie-specific normalization parameters in Population
 * @param	ini				Dictionary to input file
 * @param	pop[in,out]		Population
 * @param	timeStepMul		What multiple time step is to be used
 * @param	factor			Which multiplicative factor is used for rho
 *
 * timeStepMul=1 for whole time steps and 0.5 for half and so on.
 * This function makes use of computeRenormE() and computeRenormPhi(). See those
 * for more in-depth explanation.
 *
 * timeStepMul and factor both defaults to 1 during allocation of Population.
 * Calling this function is only necessary for other values.
 */
void pSetSpecieNorm(const dictionary *ini, Population *pop, double timeStepMul, double factor);

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

/**
 * @brief Set the same velocity to all particles
 * @param[in,out]	pop		Population
 * @param			vel		Velocity to set (expected to be pop->nDims long)
 */
void pVelSet(Population *pop, const double *vel);

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
 * @brief	Creates .pop.h5-file to store population in
 * @param	ini				Dictionary to input file
 * @param	pop[in,out]		Population
 * @param	mpiInfo			MpiInfo
 * @param	fName			Filename
 * @return	void
 * @see pWriteH5(), pCloseH5()
 *
 * An output file is created whose filename is as explained in createH5File().
 * Remember to call pCloseH5().
 *
 * The file will have one group "/pos" for position data and one group "/vel"
 * for velocity data. Each of these will have groups "specie <s>" for each
 * specie. For each time-step, the population data will be stored in a dataset
 * named "n=<timestep>" where <timestep> is signified with one decimal allowing
 * interleaved quantities.
 *
 * In PINC it is made an distinction between _non-dimensionalizing_ and
 * _normalizing_. Input quantities are non-dimensionalized by specifying them
 * in terms of Debye lengths, plasma frequency, elementary charges and so on
 * rather than using SI-units. Further on, the program normalizes them with
 * respect to for instance cell size in order to make computations as fast as
 * possible. The data stored in .pop.h5 is non-dimensionalized _and_ normalized
 * as it is often much cheaper to just rescale the axes in the visualization
 * tool rather than re-scaling all quantities in PINC.
 *
 * The file will have four attributes of size nDims attached to the root group
 * ("/") which is useful for interpreting the data. These are:
 *	- Position denormalization factor
 *	- Position dimensionalizing factor
 *	- Velocity denormalization factor
 *	- Velocity dimensionalizing factor
 *
 * The position denormalization factor can be multiplied to the integer axis to
 * convert it to be in terms of Debye lengths. Another multiplication by axis
 * dimensionalizing factor converts the axes to meters. Likewise for the
 * velocity factors.
 */
void pCreateH5(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const char *fName);

/**
 * @brief	Stores particles in Population in .pop.h5-file
 * @param	pop		Population
 * @param	mpiInfo	MpiInfo
 * @param	posN	Timestep of position data to be stored
 * @param	velN	Timestep of velocity data to be stored
 * @return			void
 * @see pCreateH5()
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

/******************************************************************************
 * DEFINED IN GRID.C
 *****************************************************************************/

/**
 * @brief Allocates a Grid object as specified in the input file
 * @param	ini			Input file
 * @param	nValues		Number of values per grid point
 * @return				Pointer to Grid
 *
 * Use nValues=1 for scalar field, nValues=3 for 3D vector field and so on.
 *
 * Remember to free using gFree().
 *
 * NB! Assumes 1 ghost point on all edges for now.
 */

Grid *gAlloc(const dictionary *ini, int nValues);

/**
 * @brief Frees allocated grid
 * @param	grid	Grid
 * @return	void
 */
void gFree(Grid *grid);

/**
 * @brief Allocates the memory for an MpiInfo struct according to input file
 * @param	ini		Input file dictionary
 * @return	Pointer to MpiInfo
 */
MpiInfo *gAllocMpi(const dictionary *ini);

/**
 * @brief Frees the memory of an MpiInfo struct
 * @param	mpiInfo		MpiInfo
 * @return	void
 */
void gFreeMpi(MpiInfo *mpiInfo);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param ini			dictionary
 * @param *Grid	Grid struct
 * @param *mpiInfo		MpiInfo struct
 * @param d				Along which dimension it should exhange ghost cells
 *
 * This function gets the ghost layers for the subdomains, by extracting the
 * outer layer of the true grid (total grid - ghost layer), and then sending it
 * to the neighboring subdomains. Then the subdomain places the ghost layers it
 * recieves in it's outer layer. It exhanges the ghost layer perpendicular to
 * the dimension it recieves as an input parameter.
 *
 * NB! Only works with 1 ghost layer.
 * @see getSendRecvSetSlice
 */
//void swapHalo(dictionary *ini, Grid *Grid, MpiInfo *mpiInfo, int d);

/**
 * @brief Set all values in grid to zero
 * @param	grid	Grid
 * @return			void
 */
void gZero(Grid *grid);

/**
 * @brief Set grid quantity to vector (or scalar) for all grid points
 * @param	grid	Grid
 * @param	value	Array (vector) of values to set
 *
 * Each grid point is set to have the vector value specified by 'value'. Hence
 * value is expected to have length grid->size[0]
 */
void gSet(Grid *grid, const double *value);

/**
 * @brief Multiply all values in grid by a double
 * @param	grid	Grid
 * @param	num		Double to multiply by
 * @return			void
 */
void gMulDouble(Grid *grid, double num);

/**
 * @brief Performs a central space finite difference on a grid
 * @param 	scalar 	Value to do the finite differencing on
 * @return	field	Field returned after derivating
 */

void gFinDiff1st(const Grid *scalar, Grid *field);

/**
 * @brief Performs a 2nd order central space finite difference on a grid
 * @param 	rho 	Value to do the finite differencing on
 * @return	phi		Field returned after derivating
 *
 *
 */
void gFinDiff2nd3D(Grid *phi,const Grid *rho);

void gFinDiff2nd2D(Grid *phi,const Grid *rho);

 /**
 * @brief Normalize E-field
 * @param	ini		Input file dictionary
 * @param	E		E-field
 *
 * Normalizes an non-normalized (but non-dimensional) E-field store in a Grid
 * object according to step-size, time step and mass and charge of specie 0.
 */
void gNormalizeE(const dictionary *ini, Grid *E);

/******************************************************************************
 * DEFINED IN IO.C
 *****************************************************************************/

/**
 * @brief Defines different types of messages
 * @see msg()
 */
typedef enum{
	STATUS = 0x00,		///< Normal status output about the progress of execution.
	WARNING = 0x01,		///< Warning. Something might not be like the user intended.
	ERROR = 0x02,		///< Error which makes the program unable to proceed. Program will stop.
	ONCE = 0x10			///< Output message from all MPI-nodes. To be bitwise ORed.
} msgKind;

/**
 * @brief	The PINC equivalent of printf().
 * @param	kind	STATUS, WARNING or ERROR depending on what to output.
 * @param	format	printf-like format specifier
 * @param	...		printf-like arguments
 * @return	void
 * @see msgKind, fMsg(), printf()
 *
 * This replaces printf() in the PINC context. Similar syntax to printf().
 * In the case of an ERROR, the program is terminated. Appends end-of-line
 * automatically at the end.
 *
 * The message will by default be printed by all nodes calling msg(), however
 * kind can be bitwise ORed with ONCE to only allow the master to display this
 * message, e.g. STATUS|ONCE.
 */
void msg(msgKind kind, const char* restrict format,...);

/**
 * @brief	Prints message to file given by msgfiles:<fNameKey> in ini-file.
 * @param	ini			Input file dictionary
 * @param	fNameKey	Name of key holding the filename
 * @param	format		printf-like format specifier
 * @param	...			printf-like arguments
 * @return	void
 * @see		msg(), printf(), fprintf()
 *
 * This replaces fprintf() in the PINC context being the preferred way to write
 * messages (not data) to files. This guarantees that no file names are hard
 * coded but configurable from the ini-file.
 *
 * Beware that this is not a high performance writing function, and very
 * frequent invocations (for instance per particle) should be avoided.
 */
void fMsg(dictionary *ini, const char* restrict fNameKey, const char* restrict format, ...);

/**
 * @brief	Opens PINC input ini-file specified in PINC's arguments.
 * @brief	argc	Argument count (as passed to PINC)
 * @param	argv	Arugment vector (as passed to PINC)
 * @return	Allocated iniparser dictionary holding data from ini-file
 *
 * Performs sanity check on argc and argv and opens the specified input file.
 * It also empties all files specified in msgfiles.
 */
dictionary* iniOpen(int argc, char *argv[]);

/**
 * @brief Get the number of elements in an array
 * @param 	ini		Dictionary to search
 * @param	key		Key string to look for
 * @return	Number of elements in entry. 0 if entry does not exist.
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
int iniGetNElements(const dictionary* ini, const char* key);

/**
 * @brief Get the array of strings associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements in returned array
 * @return			NULL-terminated array of NULL-terminated strings
 * @see				iniparser_getstring(), freeStrArr(), listToStrArr()
 *
 * Output is similar to listToStrArr(). Remember to free resulting string array
 * using freeStrArr().
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
char** iniGetStrArr(const dictionary *ini, const char *key, int *nElements);

/**
 * @brief Get the array of integers associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements in returned array
 * @return			array of integers
 * @see				iniGetLongIntArr(), iniparser_getint()
 *
 * Remember to free result array using free().
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
int* iniGetIntArr(const dictionary *ini, const char *key, int *nElements);

/**
 * @brief Get the array of long integers associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements in returned array
 * @return			array of integers
 * @see				iniGetIntArr()
 *
 * Remember to free result array using free().
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
long int* iniGetLongIntArr(const dictionary *ini, const char *key, int *nElements);

/**
 * @brief Get the array of doubles associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements in returned array
 * @return			array of doubles
 * @see				iniparser_getdouble()
 *
 * Remember to free result array using free().
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
double* iniGetDoubleArr(const dictionary *ini, const char *key, int *nElements);

/**
 * @brief Assert that a number of entries are arrays of equal length.
 * @param			ini			Dictionary to search
 * @param			nKey		Number of keys to search for
 * @param			...			Keys to search for
 * @return			Number of elements in arrays
 * @see				iniGetNElements()
 *
 * Example:
 * @code
 *	iniAssertEqualNElements(ini,3,"mySec:a","mySec:b","mySec:c");
 * @endcode
 *
 * This code does nothing if the specified entries have equal length. Unequal
 * lengths will lead to an error.
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
int iniAssertEqualNElements(const dictionary *ini, int nKeys, ...);

/**
 * @brief Frees dynamically allocated NULL-terminated array of strings
 * @param	strArr	Pointer to array of strings
 * @return	void
 * @see listToStrArr(), iniGetStrArr()
 */
void freeStrArr(char** strArr);

/**
 * @brief Creates .h5-file
 * @param	ini		Input file dictionary
 * @param	fName	File name
 * @param	fSubExt	File sub-extension
 * @return	HDF5 file identifier
 * @see		H5Fcreate(), H5Fclose()
 *
 * This is the PINC equivalent of H5Fcreate() to create HDF5 files. The file
 * will be named <fName>.<fSubExt>.h5 or <prefix>_<fName>.<fSubExt>.h5.
 * Similarly as .h5 indicates the file type being h5, <fSubExt> indicates
 * _what kind_ of .h5-file it is, i.e. what kind of data it contains, typically
 * grid data (grid) or population data (pop).
 *
 * fName is a name which is hard-coded to further specify not only the _kind_ of
 * data but exactly which quantity. E.g. a file named "rho.grid.h5" specifies
 * charge density whereas "E.grid.h5" is for electric field.
 *
 * The file will be stored in the folder specified by "files:output" in the
 * input file. Examples of valid values of "files:output":
 *	- output = data/
 *  - output = ./data/
 *	- output = ../data/
 *	- output = ~/data/
 *	- output = /home/me/data/
 *
 * However, this variable also allows the user to specify a prefix to each file.
 * This is indicated by _not_ ending the variable with '/'. For instance:
 *	- output = prefix
 *	- output = data/prefix
 *
 * This will output files such a "prefix_rho.grid.h5".
 *
 * Contrary to H5Fcreate() this function creates parent directories unless they
 * already exists. If a file already exists it will fail but contrary to
 * H5Fcreate() it will fail gracefully with an ERROR.
 *
 * Close return value using H5Fclose().
 */
hid_t createH5File(const dictionary* ini, const char *fName, const char *fSubExt);

/**
 * @brief	Assign particles artificial positions suitable for debugging
 * @param			ini				Input file dictionary
 * @param[in,out]	grid			Grid
 *
 * A quantity will be artificially assigned values depending on the position of
 * the grid points. For scalar valued quantities, node j will have value j,
 * node (j,k) will have value j+k*10 and so on for higher dimensions. For
 * instance at node (j,k,l)=(4,5,6) the value will be 456.
 *
 * For vector valued quantities the integer part of all values at a given point
 * are the same (unless there is 10 values or more), but the decimal part
 * increments by 0.1 for each value.
 *
 * In addition, 1000*mpiRank is added to all values. For instance the third
 * value (z-component) at grid point (j,k,l)=(4,5,6) of the subdomain with MPI
 * rank 3 will be 3654.3.
 */
void gValDebug(Grid *grid, const MpiInfo *mpiInfo);

/**
 * @brief	Creates .grid.h5-file to store population in
 * @param	ini				Dictionary to input file
 * @param	grid			Grid
 * @param	mpiInfo			MpiInfo
 * @param	denorm			Quantity denormalization factors
 * @param	dimen			Quantity dimensionalizing factors
 * @param	fName			Filename
 * @return	void
 * @see gWriteH5(), gCloseH5()
 *
 * An output file is created whose filename is as explained in createH5File().
 * Remember to call gCloseH5().
 *
 * The file will have one dataset in the root group for each time-step a grid
 * quantity is stored, named "n=<timestep>" where <timestep> is signified with
 * one decimal allowing for interleaved quantities.
 *
 * In PINC it is made an distinction between _non-dimensionalizing_ and
 * _normalizing_. Input quantities are non-dimensionalized by specifying them
 * in terms of Debye lengths, plasma frequency, elementary charges and so on
 * rather than using SI-units. Further on, the program normalizes them with
 * respect to for instance cell size in order to make computations as fast as
 * possible. The data stored in .grid.h5 is non-dimensionalized _and_ normalized
 * as it is often much cheaper to just rescale the axis in the visualization
 * tool rather than re-scaling all quantities in PINC.
 *
 * The file will have four attributes of size nDims attached to the root group
 * ("/") which is useful for interpreting the data. These are:
 *	- Axis denormalization factor
 *	- Axis dimensionalizing factor
 *	- Quantity denormalization factor
 *	- Quantity dimensionalizing factor
 *
 * The axis denormalization factor can be multiplied to the integer axis to
 * convert it to be in terms of Debye lengths. Another multiplication by axis
 * dimensionalizing factor converts the axes to meters. Likewise for the
 * quantity factors. However, since the quantity factors depend on which
 * quantity it is (e.g. charge density or electric field), it must be specified
 * in the inputs denorm and dimen in this function. They are expected to be of
 * length nDims.
 */
void gCreateH5(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo, const double *denorm, const double *dimen, const char *fName);

/**
 * @brief	Stores quantity in Grid in .grid.h5-file
 * @param	grid			Grid
 * @param	mpiInfo			MpiInfo
 * @param	n				Timestep of quantity to be stored
 * @return	void
 * @see gCreateH5()
 *
 * The position and velocity of all particles are stored, referred to global
 * reference frame. The function takes care of merging the particles from all
 * MPI nodes to one file.
 *
 */
void gWriteH5(const Grid *grid, const MpiInfo *mpiInfo, double n);

/**
 * @brief	Closes .grid.h5-file
 * @param	grid		Grid
 * @return	void
 */
void gCloseH5(Grid *grid);

/******************************************************************************
 * DEFINED IN AUX.C
 *****************************************************************************/

/**
 * @brief	Allocates a Timer struct
 * @param	rank	Rank of the node where the timer is active (-1 for all)
 * @return	Pointer to Timer struct
 * @see		Timer, freeTimer, tMsg()
 *
 * Remember to free using freeTimer().
 */
Timer *allocTimer(int rank);

/**
 * @brief	Frees a Timer struct allocated with allocTimer()
 * @param	timer 	Pointer to Timer struct
 * @see		Timer, allocTimer()
 */
void freeTimer(Timer *timer);

/**
 * @brief	Prints a message along with timing information
 * @param	timer		Pointer to timer
 * @param	format		printf-like format specifier
 * @param	...			printf-like arguments
 * @return	void
 * @see 	Timer, allocTimer(), printf()
 *
 * Useful for testing execution speed of chunks of code. Each call to tMsg()
 * prints the total time the program has been running, along with the time since
 * last call to tMsg() before it resets the timer.
 *
 * To reset the timer without printing set format=NULL.
 *
 * Only MPI nodes for which the timer is activated by allocTimer() will print
 * the messages.
 */
void tMsg(Timer *timer, const char *restrict format, ...);

/**
 * @brief Returns the product of all elements in an integer array
 * @param	a			Pointer to array
 * @param	nElements	Number of elements in array
 * @return	Product
 */
int intArrProd(const int *a, int nElements);

/**
 * @brief Returns the product of all elements in a double array
 * @param	a			Pointer to array
 * @param	nElements	Number of elements in array
 * @return	Product
 */
double doubleArrProd(const double *a, int nElements);

/**
 * @brief Returns the cumulative product of the elements in an integer array
 * @param	a			pointer to array
 * @param	nElements	Number of elements in array
 * @return	Pointer to allocated array of size nElements+1
 *
 * The result will be given by:
 * @code
 *	result[0] = 1;
 *	result[1] = a[0];
 *	result[2] = a[0]*a[1];
 *	...
 * @endcode
 */
int *intArrCumProd(const int *a, int nElements);

/**
 * @brief Returns the cumulative product of the elements in an long integer array
 * @param	a			pointer to array
 * @param	nElements	Number of elements in array
 * @return	Pointer to allocated array of size nElements+1
 *
 * The result will be given by:
 * @code
 *	result[0] = 1;
 *	result[1] = a[0];
 *	result[2] = a[0]*a[1];
 *	...
 * @endcode
 */
long int *longIntArrCumProd(const int *a, int nElements);

/**
 * @brief Returns the product of all elements in an integer array
 * @param	a			pointer to array
 * @param	b			pointer to array
 * @param	nElements	Number of elements in arrays
 * @return	Pointer to allocated array of size nElements
 */
int *intArrMul(const int *a, const int *b, int nElements);

/**
 * @brief Concatenates strings
 * @param	n	Number of strings to concatenate
 * @param	...	List of pointer to strings to concatenate
 * @return		Pointer to concatenated string
 * @see	strcat()
 *
 * Contrary to strcat() this function allocates a new string of suitable length.
 * Make sure to call free().
 */
char *strAllocCat(int n,...);

#endif // PINC_H
