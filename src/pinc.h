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
} Population;

/****************************************************************************
***			Defined in grid.c
****************************************************************************/
/**
 * @brief Contains specification on how the grid is structured and decomposed.
 *
 * The total simulation domain can be split across several MPI nodes where the
 * position J,K,L (in case of 3D) of each sub-domain is stored in node. nNodes
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
	int nDims;					///< Number of dimensions (usually 3)
	int *nGPoints;				///< The number of grid points per dimension (nDims elements)
	int *nGPointsProd;			///< Cumulative product of nGPoints (nDims+1 elements)
	int *nGhosts;				///< Number of ghost grid points (2*nDims elements)
	int *node;					///< MPI node (nDims elements)
	int *nNodes;				///< Number of MPI nodes (nDims elements)
	int *offset;				///< Offset from global reference frame (nDims elements)
	double *posToNode;			///< Factor for converting position to node (nDims elements)
} Grid;

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
	double *halo;				///< Ghost cells surrounding true cells
	int nValues;				///< Number of values per grid point (usually 1 or 3)
	Grid *grid;					///< Specifications of the grid
} GridQuantity;

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
 * @see Grid
 * @see gridQuantity
 */

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues);

/**
 * @brief Frees the memory of a GridQuantity struct
 * @param gridQuantity 		Gridquantity struct
 *
 * Frees the memory of the GridQuantity struct, since the Grid member of the struct
 * can be shared by several gridQuantity'ies it needs to be freed seperately
 * @see freeGrid
 */

void freeGridQuantity(GridQuantity *gridQuantity);

/**
 * @brief Gather the edges of a grid and returns a vector with them
 * @param 	ini 				dictionary of the input file
 * @param 	gridQuantity		GridQuantity struct
 * @return 	ghostEdge 			vector containing ghost values (*double)
 *
 * From a grid this function gathers the ghost layer, stores it in a 1D array
 * and returns it.
 *
 * In the case where the grid has several dimensions first the lower edge is stored
 * for all dimensions, and then the upper edge is stored for each dimension
 * So for a 2D case the vector looks like:
 *		\f[
 * 		Edge = [\partial \vec{x}_{min}\;\;|\;\;\partial \vec{y}_{min}\;\;|\;\;\partial \vec{x}_{max}
 *				\;\;|\;\;\partial \vec{y}_{max}]
 *		\f]
 * In 3D it will be:
 * 		\f[
 *			Edges = [\;\; \partial \vec{x}_{min} \;\;|\;\; \partial \vec{y}_min \;\;|\;\; \partial \vec{z}_{min}
 *					\;\;|\;\; \partial \vec{x}_{max} \;\;|\;\;  \partial \vec{y}_{max} \;\;|\;\;
 *					\partial \vec{z}_{max} ]
 *		\f]
 *
 * Note: 	ghostEdge vector should be considered for a member of grid structs
 * 			to avoid allocating it each time
 *
 */
double *getHalo(dictionary *ini, GridQuantity *gridQuantity);

/**
 * @brief 	Places the ghost vector on the grid again after swapping
 *
 * More documentation TBD
 */
void distributeHalo(dictionary *ini, GridQuantity *gridQuantity);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param dictionary 	*ini
 * @param GridQuantity 	*gridQuantity
 *
 * TBD
 */

 void swapHalo(dictionary *ini, GridQuantity *gridQuantity);

/**
 * @brief Writes information about the grid structs to a parsefile
 * @param ini 		dictionary of the input file
 * @param grid 		grid struct
 * @param nValues	number of values per grid point.
 * @return void
 *
 * TBD
 * 2 subdomains 2D grid example:
 * @code
	111111	222222
	111111	222222
	111111	222222
	111111	222222

	swapGhosts(ini, gridQuantity);

	111112	122222
	111112	122222
	111112	122222
	111112	122222
  @endcode
 */
void gridParseDump(dictionary *ini, Grid *grid, GridQuantity *gridQuantity);



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
 * @brief	kind	STATUS, WARNING or ERROR depending on what to output.
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

/******************************************************************************
 * DEFINED IN AUX.C
 *****************************************************************************/

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
 * @brief Returns the product of all elements in an integer array
 * @param	a			pointer to array
 * @param	b			pointer to array
 * @param	nElements	Number of elements in arrays
 * @return	Pointer to allocated array of size nElements
 */
int *intArrMul(const int *a, const int *b, int nElements);

/******************************************************************************
 * DEFINED IN POPULATION.C
 *****************************************************************************/
Population *allocPopulation(const dictionary *ini);
void freePopulation(Population *pop);
/******************************************************************************
 * DEFINED IN GRID.C
 *****************************************************************************/



#endif // PINC_H
