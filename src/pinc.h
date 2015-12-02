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
 * The position and velocity is stored in a flat manner, such that (x,y,z) of
 * particle 0 comes first, then (x,y,z) of particle 1, and so on (assuming 3D).
 * As an example, printing the (x,y,z) position of the first 3 particles:
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
 * and stop at i=iStop[s]. Due to particle increase/decrease there is allocated
 * space for more particles than are present. Thus the _allocated_ space for
 * specie s start at i=iStart[s] and stop at i=iStart[s+1]-1. For convenience,
 * iStart has nSpecies+1 elements such that this is true also for the last
 * specie. The last element is then simply the number of particles allocated
 * in total.
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
	long int *iStop;	///< Last index of specie s (nSpecies elements)
	double *q;			///< Specie charge [elementary charges] (nSpecies elements)
	double *m;			///< Specie mass [electron masses] (nSpecies elements)
	int nSpecies;		///< Number of species
	int nDims;			///< Number of dimensions (usually 3)
	hid_t h5;			///< HDF5 file handler
} Population;

/**
 * @brief Contains information regarding how the PIC code is parallelized.
 *
 * The total simulation domain can be split across several subdomains, one for
 * each MPI node, where position J,K,L (in case of 3D) of each subdomain is
 * stored in subdomain. nSubdomains represents the number of subdomains along
 * each dimension, and nSubdomainsProd is the cumulative product of nSubdomains
 * similarly as nGPointsProd in Grid.
 * @see Grid
 *
 * Position of particles and objects is normalized with respect to stepsize
 * and can be specified in a global reference frame, or local to one subdomain.
 * The local reference frame is defined such that a particle with integer
 * position (j,k,l) would be located _on_ grid point (j,k,l). This
 * makes it fast to determine the index of the nodes surrounding a particle.
 *
 * offset is the offset of this MPI node's subdomain with respect to the global
 * reference frame. Adding/subtracting this to a position converts to/from the
 * global reference frame.
 * @see toLocalFrame(), toGlobalFrame()
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
 * @brief Contains specification of the grid.
 *
 * nDims and nGPoints specifies the number of dimensions and the number of grid
 * points along each of them. In the case of decomposing the domain into sub-
 * domains using MPI, each MPI node will have its local grid. The number of
 * grid points include true grid points plus ghost grid points that are copied
 * from the neighbouring MPI node. The number of ghost points copied from the
 * neighbours is given in nGhosts. Consider the following 1D example grid:
 *
 * 	g	x	x	x	x	g	g
 *
 * It consists of nGPoints[0]=7 grid points but the leftmost and the two right-
 * most are simply copied from the neighbours. Hence, nGhosts={1,2}. For
 * several dimensions nGPoints first lists all ghosts on the lower edge along
 * all dimensions then all ghosts on the upper edge.
 *
 * Finally, nGPointsProd contains the cumulative products of nGPoints:
 *
 * @code
 *	nGPointsProd[0]=1;
 *	nGPointsProd[1]=1*nGPoints[0];
 *	nGPointsProd[2]=1*nGPoints[0]*nGPoints[1];
 *  ...
 * @endcode
 *
 * nGPointsProd[nDims] will be the total number of grid points.
 *
 * @see allocGrid
 * @see GridQuantity
 */
typedef struct{
	int nDims;					///< Number of dimensions (usually 1-3)
	int *nGPoints;				///< The number of grid points per dimension (nDims elements)
	int *nTGPoints;				///< The number of true grid points (nDims elements)
	long int *nGPointsProd;		///< Cumulative product of nGPoints (nDims+1 elements)
	int *nGhosts;				///< Number of ghost layers in grid (2*nDims elements)
	double *dr;					///< Step-size (nDims elements)
} Grid;

/**
 * @brief A quantity defined on a grid, for instance charge density.
 *
 * Can represent both scalar and vector quantities. For scalar quantities
 * nValues=1, whereas vector fields typically have nValues=nDims=3. The nodes
 * are lexicographically ordered, i.e. node (0,0,0), (1,0,0), (2,0,0), ...,
 * (0,1,0), .... Thus for a 3-dimensional scalar quantity rho, the field value
 * at node (j,k,l) is accessed in the following manner:
 *
 * @code
 *	Grid *rho;
 *	...
 *	int *nGPointsProd = rho->grid->nGPointsProd;
 *	int p = j + k*nGPointsProd[1] + l*nGPointsProd[2];
 *	printf("rho(%i,%i,%i)=%f\n",j,k,l,rho->val[p]);
 * @endcode
 *
 * Note that nGPointsProd[d] represents an increment of 1 in direction d.
 *
 * For vector fields all values are stored one node at a time. For instance, to
 * print the vector value of E at node (j,k,l):
 *
 * @code
 *	Grid E;
 *	...
 *  int *nGPointsProd = E->grid->nGPointsProd;
 *  int p = j + k*nGPointsProd[1] + l*nGPointsProd[2];
 *	double *pVal = &E.val[p*E.nValues];
 *	printf("The field at node (%i,%i,%i) is (%f,%f,%f).\n",j,k,l,pVal[0],pVal[1],pVal[2]);
 * @endcode
 *
 * slice is a buffer which can be used to temporarily store an (N-1)-dimensional
 * slice of the grid quantity.
 *
 * The variables starting with h5 is used to store information regarding writing
 * to h5 output file if such a file is created.
 * @see createGridQuantityH5()
 */
typedef struct{
	double *val;		///< The values on the grid
	double *slice;		///< Slice buffer of the grid sent to other subdomains
	int nValues;		///< Number of values per grid point (usually 1 or 3)
	Grid *grid;			///< Specifications of the grid
	hid_t h5;			///< HDF5 file handler
	hid_t h5MemSpace;	///< HDF5 memory space description
	hid_t h5FileSpace;	///< HDF5 file space description
} GridQuantity;

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
Population *allocPopulation(const dictionary *ini);

/**
 * @brief	Frees memory for Population
 * @param	pop		Pointer to population to be freed
 * @see allocPopulation()
 */
void freePopulation(Population *pop);

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
 * @see velMaxwell()
 */
void posUniform(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const gsl_rng *rng);

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
void posDebug(const dictionary *ini, Population *pop);

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
void velMaxwell(const dictionary *ini, Population *pop, const gsl_rng *rng);

/**
 * @brief	Creates .pop.h5-file to store population in
 * @param	ini		Dictionary to input file
 * @param	pop		Population
 * @param	mpiInfo	MpiInfo
 * @param	fName	Filename
 * @return	void
 * @see writePopulationH5(), closePopulationH5()
 *
 * An output file is created whose filename is as explained in createH5File().
 * Remember to call closePopulationH5().
 *
 * The file will have one group "/pos" for position data and one group "/vel"
 * for velocity data. Each of these will have groups "specie <s>" for each
 * specie. For each time-stepd, the population data will be stored in a dataset
 * named "n=<timestep>" where <timestip> is signified with one decimal allowing
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
 *	Position denormalization factor
 *	Position dimensionalizing factor
 *	Velocity denormalization factor
 *	Velocity dimensionalizing factor
 *
 * The position denormalization factor can be multiplied to the integer axis to
 * convert it to be in terms of Debye lengths. Another multiplication by axis
 * dimensionalizing factor converts the axes to meters. Likewise for the
 * velocity factors.
 */
void createPopulationH5(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const char *fName);

/**
 * @brief	Stores particles in Population in .pop.h5-file
 * @param	pop		Population
 * @param	mpiInfo	MpiInfo
 * @param	posN	Timestep of position data to be stored
 * @param	velN	Timestep of velocity data to be stored
 * @return			void
 * @see createPopulationH5()
 *
 * The position and velocity of all particles are stored, referred to global
 * reference frame. The function takes care of merging the particles from all
 * MPI nodes to one file.
 *
 * NB: pop is not constified because all particles are transformed to global
 * reference frame before writing to .h5-file. However, they are transferred
 * back to local reference frame after writing so pop should remain unchanged to
 * within numerical precision.
 */
void writePopulationH5(Population *pop, const MpiInfo *mpiInfo, double posN, double velN);

/**
 * @brief	Closes .pop.h5-file
 * @param	pop		Population
 * @return	void
 */
void closePopulationH5(Population *pop);

/******************************************************************************
 * DEFINED IN GRID.C
 *****************************************************************************/

/**
 * @brief Allocates the memory for the a grid struct as specified in the input file
 * @param	ini		input file
 * @return	Pointer to Grid
 *
 * Remember to free using freeGrid().
 *
 * @see Grid
 */

Grid *allocGrid(const dictionary *ini);

/**
 * @brief Frees allocated grid
 * @param	grid	Grid
 * @return	void
 */
void freeGrid(Grid *grid);

/**
 * @brief Allocates the memory for a GridQuantity as specified in the input file
 * @param	ini 	dictionary of the input file
 * @param	grid 	Grid
 * @param	nValues	Number of values per grid point.
 * @return	Pointer to GridQuantity
 *
 * Remember to free using freeGridQuantity().
 *
 * @see Grid
 * @see gridQuantity
 */
GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues);

/**
 * @brief Frees the memory of a GridQuantity struct
 * @param	gridQuantity	GridQuantity
 * @return	void
 */
void freeGridQuantity(GridQuantity *gridQuantity);

/**
 * @brief Allocates the memory for an MpiInfo struct according to input file
 * @param	ini		Input file dictionary
 * @return	Pointer to MpiInfo
 */
MpiInfo *allocMpiInfo(const dictionary *ini);

/**
 * @brief Frees the memory of an MpiInfo struct
 * @param	mpiInfo		MpiInfo
 * @return	void
 */
void freeMpiInfo(MpiInfo *mpiInfo);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param	ini				Input file dictionary
 * @param	gridQuantity 	GridQuantity
 *
 * TBD
 * Frees the memory of the GridQuantity struct, since the Grid member of the struct
 * can be shared by several gridQuantity'ies it needs to be freed seperately
 * @see freeGrid
 */
void swapHalo(dictionary *ini, GridQuantity *gridQuantity);

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
 *	output = data/
 *  output = ./data/
 *	output = ../data/
 *	output = ~/data/
 *	output = /home/me/data/
 *
 * However, this variable also allows the user to specify a prefix to each file.
 * This is indicated by _not_ ending the variable with '/'. For instance:
 *	output = prefix
 *	output = data/prefix
 * will output files such a "prefix_rho.grid.h5".
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
 * @param[in,out]	gridQuantity	GridQuantity
 *
 * A quantity will be artificially assigned values depending on the position of
 * the grid points. For scalar valued quantities, node j will have value j,
 * node (j,k) will have value j*10+k and so on for higher dimensions. For
 * instance at node (j,k,l)=(4,5,6) the value will be 456.
 *
 * For vector valued quantities the integer part of all values at a given point
 * are the same (unless there is 10 values or more), but the decimal part
 * increments by 0.1 for each value.
 *
 * In addition, 1000*mpiRank is added to all values. For instance the third
 * value (z-component) at grid point (j,k,l)=(4,5,6) of the subdomain with MPI
 * rank 3 will be 3456.3.
 */
void gridValDebug(GridQuantity *gridQuantity, const MpiInfo *mpiInfo);

/**
 * @brief	Creates .grid.h5-file to store population in
 * @param	ini				Dictionary to input file
 * @param	gridQuantity	GridQuantity
 * @param	mpiInfo			MpiInfo
 * @param	denorm			Denormalization factors
 * @param	dimen			Dimensionalizing factors
 * @param	fName			Filename
 * @return	void
 * @see writeGridQuantityH5(), closeGridQuantityH5()
 *
 * An output file is created whose filename is as explained in createH5File().
 * Remember to call closeGridQuantityH5().
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
 *	Axis denormalization factor
 *	Axis dimensionalizing factor
 *	Quantity denormalization factor
 *	Quantity dimensionalizing factor
 *
 * The axis denormalization factor can be multiplied to the integer axis to
 * convert it to be in terms of Debye lengths. Another multiplication by axis
 * dimensionalizing factor converts the axes to meters. Likewise for the
 * quantity factors. However, since the quantity factors depend on which
 * quantity it is (e.g. charge density or electric field), it must be specified
 * in the inputs denorm and dimen in this function. They are expected to be of
 * length nDims.
 */
void createGridQuantityH5(const dictionary *ini, GridQuantity *gridQuantity, const MpiInfo *mpiInfo, const double *denorm, const double *dimen, const char *fName);

/**
 * @brief	Stores quantity in gridQuantity in .grid.h5-file
 * @param	gridQuantity	GridQuantity
 * @param	mpiInfo			MpiInfo
 * @param	n				Timestep of quantity to be stored
 * @return	void
 * @see createGridQuantityH5()
 *
 * The position and velocity of all particles are stored, referred to global
 * reference frame. The function takes care of merging the particles from all
 * MPI nodes to one file.
 *
 */
void writeGridQuantityH5(const GridQuantity *gridQuantity, const MpiInfo *mpiInfo, double n);

/**
 * @brief	Closes .grid.h5-file
 * @param	gridQuantity	GridQuantity
 * @return	void
 */
void closeGridQuantityH5(GridQuantity *gridQuantity);

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
 * @brief Makes all parent directories of URL path
 * @param	path	Path
 * @return	0 for success, 1 for failure
 *
 * Examples:
 *	path="dir/dir/file"	generates the folder "./dir/dir/"
 *  path="dir/dir/dir/" generates the folder "./dir/dir/dir/"
 *  path="../dir/file" generates the folder "../dir/"
 *	path="/dir/file" generates the folder "/dir/"
 *
 * Already existing folders are left as-is. This function can be used to ensure
 * that the parent directories of its path exists.
 */
int makeParentPath(const char *path);

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
