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
	int mpiRank;		///< MPI Rank
	int mpiSize;		///< MPI Size
	hid_t h5;			///< h5-file
 } Population;

 /**
  * @brief Contains specification on how the grid is structured and decomposed.
  *
  * The total simulation domain can be split across several MPI nodes where the
  * position J,K,L (in case of 3D) of each sub-domain is stored in node. nSubdomains
  * represents the number of sub-domains along each dimension.
  *
  * nDims and nGPoints specifies the number of dimensions and the number of grid
  * points along each of them in this node. The number of grid points include
  * true grid points plus ghost grid points that are copied from the
  * neighbouring MPI node. The number of ghost points copied from the neighbours
  * is given in nGhosts. Consider the following 1D example grid:
  *
  * 	g	x	x	x	x	g	g
  *
  * It consists of nGPoints[0]=7 grid points but the leftmost and the two right-
  * most are simply copied from the neighbours. Hence, nGhosts={1,2}. For
  * several dimensions nGPoints first lists all ghosts on the lower edge along
  * all dimensions then all ghosts on the upper edge.
  *
  * Position of particles and objects is normalized with respect to stepsize
  * and can be specified in a global reference frame, or one local to the node.
  * The local reference frame is defined such that a particle with integer
  * position (j,k,l) would be located _on_ grid point (j,k,l). This makes it
  * fast to determine the index of the nodes surrounding a particle. The offset
  * between the global and the local reference frames is stored in offset for
  * easy conversion, and the variable posToNode is handy for determining the
  * node to which a particle belong (if it exits the sub-domain of this node),
  * e.g. for 1D:
  *
  * @code
  *	J = (int)(posToNode[0]*pos[0]);
  * @endcode
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
  * nGPointsProd[nDims-1] will be the total number of grid points.
  *
  * @see allocGrid
  * @see GridQuantity
  */
typedef struct{
	int nDims;					///< Number of dimensions (usually 1-3)
	int *nGPoints;				///< The number of grid points per dimension (nDims elements)
	long int *nGPointsProd;		///< Cumulative product of nGPoints (nDims+1 elements)
	int *nGhosts;				///< Number of ghost layers in grid
	double *dr;					///< Step-size (nDims elements)
} Grid;

typedef struct{
	int mpiRank;				///< MPI rank
	int mpiSize;				///< MPI size
	int *subdomain;				///< MPI node (nDims elements)
	int *nSubdomains;			///< Number of MPI nodes (nDims elements)
	int *nSubdomainsProd;		///< Cumulative product of nSubdomains
	int *offset;				///< Offset from global reference frame (nDims elements)
	double *posToSubdomain;		///< Factor for converting position to subdomain (nDims elements)
} MpiInfo;

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
  *  int *nGPointsProd = rho->grid->nGPointsProd;
  *  int p = j + k*nGPointsProd[1] + l*nGPointsProd[2];
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
  */
typedef struct{
	double *val;				///< The values on the grid
	double *slice;			///< Slices of the grid sent to other subdomains
	int nValues;				///< Number of values per grid point (usually 1 or 3)
	Grid *grid;					///< Specifications of the grid
	hid_t h5;					///< h5-file
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
 * @param			rng		Random number generator
 * @return			void
 * @see velMaxwell()
 *
 * The amount of particles specified by population:nParticles in ini will be
 * generated with uniformly distributed random positions within the simulation
 * domain (global reference frame). In case of multiple MPI nodes only particles
 * residing in the given MPI node's sub-domain will be stored, and will be
 * transformed to its local reference frame. The rng should have the same seed
 * on all MPI nodes when calling this function as that will ensure that all
 * nodes generates the same particles and discards particles not belonging to
 * their sub-domain consecutively. Failure to do so may lead to the number of
 * particles generated being slightly different than specified in ini.
 *
 * Beware that this function do not assign any velocity to the particles.
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
 * @see velMaxwell()
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
 * @param	ini		Dictionary to input file specifying file name
 * @param	pop		Population
 * @return			HDF5 file identifier
 */
void createPopulationH5(const dictionary *ini, Population *pop);

/**
 * @brief	Stores particles in Population in file
 * @param	pop		Population
 * @param	file	.pop.h5-file created with h5openPopulation()
 * @param	n		Timestep to store as
 * @return			void
 * @see	h5openPopulation()
 */
void writePopulationH5(Population *pop, double posN, double velN);

/******************************************************************************
 * DEFINED IN GRID.C
 *****************************************************************************/

/**
 * @brief Allocates the memory for the a grid struct as specified in the input file
 * @param dictionary *ini
 * @return Grid *grid
 *
 * Allocates the memory for a grid struct and returns a pointer to it
 *
 * @see Grid
 */

Grid *allocGrid(const dictionary *ini);

/**
 * @brief Frees a grid struct (TBD)
 */

void freeGrid(Grid *grid);

/**
 * @brief Allocates the memory for the a grid struct as specified in the input file
 * @param ini 				dictionary of the input file
 * @param grid 				grid struct
 * @param nValues			number of values per grid point.
 * @return gridQuantity 	GridQuantity struct
 *
 *
 *
 * Allocates the memory for a GridQuantity struct and returns a pointer to it
 *
 *
 * NB! Assumes 1 ghost point on all edges for now.
 * @see Grid
 * @see gridQuantity
 */

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param ini			dictionary
 * @param *gridQuantity	GridQuantity struct
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

void swapHalo(dictionary *ini, GridQuantity *gridQuantity, MpiInfo *mpiInfo, int d);

/**
 * @brief Frees the memory of a GridQuantity struct
 * @param gridQuantity 		Gridquantity struct
 * @return void
 *
 * Frees the memory of the GridQuantity struct, since the Grid member of the struct
 * can be shared by several gridQuantity'ies it needs to be freed seperately
 * @see freeGrid
 */

void freeGridQuantity(GridQuantity *gridQuantity);

void freeMpiInfo(MpiInfo *mpiInfo);
MpiInfo *allocMpiInfo(const dictionary *ini);

/******************************************************************************
 *		DEFINED IN MULTIGRID.C
 *****************************************************************************/


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
 * @param	format	Error message with specification of how to interpret data.
 * @param	...		Data to be interpreted in message.
 * @return	void
 * @see msgKind, fMsg(), printf()
 *
 * This replaces printf() in the PINC context. Similar syntax to printf().
 * In the case of an ERROR, the program is terminated. Appends end-of-line
 * automatically at the end.
 */
void msg(msgKind kind, const char* restrict format,...);

/**
 * @brief	Prints message to file given by msgfiles:<fNameKey> in ini-file.
 * @param	ini			ini-file dictionary
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
 * Beware that this is not high performance writing function, and very frequent
 * invocations (for instance per particle) should be avoided.
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
 * @brief	Completes the time-section in the ini-file.
 * @param	ini		ini-file dictionary
 * @return	void
 * @see		ini_complete_grid()
 *
 * Computes the missing parameter Nt, T or dt from the other two and adds
 * it to the dictionary. Sanity checks included.
 */
//void ini_complete_time(dictionary *ini);

/**
 * @brief	Completes the grid-section in the ini-file.
 * @param	ini		ini-file dictionary
 * @return	void
 * @see		ini_complete_time()
 *
 * Computes the missing parameter Ng, L or dx from the other two and adds
 * it to the dictionary. Sanity checks included.
 */
//void ini_complete_grid(dictionary *ini);

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
 * @see				iniGetIntArr(), iniGetDoubleArr(), iniparser_getstring()
 * @see				freeStrArr(), listToStrArr()
 *
 * Output is similar to listToStrArr(). Remember to free result string array
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
 * @see				iniGetLongIntArr(), iniGetStrArr(), iniGetDoubleArr(), iniparser_getint()
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
 * @return			array of idoubles
 * @see				iniGetStrArr(), iniGetIntArr(), iniparser_getdouble()
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
 * @see listToStrArr()
 */
void freeStrArr(char** strArr);

/**
 * @brief Creates .h5-file
 * @param	fName	File name
 * @param	fcpList	File Creation Property List
 * @param	acpList	File Access Property List
 * @return	HDF5 file identifier
 * @see		H5Fcreate(), H5Fclose()
 *
 * Contrary to H5Fcreate() this function creates parent directories unless they
 * already exists. If a file already exists it will fail but contrary to
 * H5Fcreate() it will fail gracefully with an ERROR.
 *
 * Close file with H5Fclose() as usual.
 */
hid_t createH5File(const char *fName, hid_t fcpList, hid_t fapList);

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
 * @brief Concatenates two strings
 * @param	a	First string
 * @param	b	Second string
 * @return		Pointer to concatenated string
 * @see	strcat()
 *
 * Contrary to strcat() this function allocates a new string of suitable length.
 * Make sure to run free().
 */
char *strAllocCat(const char *a, const char *b);


#endif // PINC_H
