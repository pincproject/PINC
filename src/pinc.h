/**
 * @file		pinc.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
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
 * particle 0 comes first, then (x,y,z) of particle 1, and so on (in the case
 * of nDim=3 dimensions). As an example, printing the (x,y,z) position of the
 * first 3 particles:
 *
 * @code
 *	Population pop;
 *  ...
 *	for(int i=0;i<3;i++){
 *		double *iPos = &pop.pos[i*pop.nDim];
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
 * energy[i] simply holds the energy of particle i, or if energy is not
 * computed, energy points to NULL.
 *
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

/**
 * @brief Contains grid functions, for instance charge density or E-field.
 * 
 * Can represent both scalar fields and vector fields. For scalar fields
 * nValues=1, whereas vector fields typically have nValues=nDim=3. The nodes
 * are lexicographically ordered, i.e.\ node (0,0,0), (1,0,0), (2,0,0), ...,
 * (0,1,0), .... Thus for a 3-dimensional scalar field rho, the field value at
 * node (j,k,l) is accessed in the following manner:
 *
 * @code
 *	Grid rho;
 *	...
 *  int p = j + k*rho.nNodes[0] + l*rho.nNodes[0]*rho.nNodes[1];
 *	printf("rho(%i,%i,%i)=%f\n",j,k,l,rho.val[p]);
 * @endcode
 *
 * For convenience, nNodesProd contains the cumulative product of nNodes:
 *
 * @code
 *	nNodesProd[0]=1;						// Corresponding to a j-increment
 *	nNodesProd[1]=1*nNodes[0];				// Corresponding to a k-increment
 *	nNodesProd[2]=1*nNodes[0]*nNodes[1];	// Corresponding to an l-increment
 *  ...
 * @endcode
 *
 * with the last element being the total number of nodes. Note that
 * nNodesProd[d] corresponds to an increment in direction d, and that this
 * can be useful for saving operations. For instance if rho(j,k,l) and
 * rho(j,k+1,l) is to be printed:
 *
 * @code
 *  int p = j + k*rho.nNodesProd[1] + l*rho.nNodesProd[2];
 *	printf("rho(%i,%i,%i)=%f\n",j,k,l,rho.val[p]);
 *  p += rho.nNodesProd[1];
 *	printf("rho(%i,%i,%i)=%f\n",j,k+1,l,rho.val[p]);
 * @endcode
 * 
 * For vector fields all values are stored one node at a time. For instance, to
 * print the vector value of E at node (j,k,l):
 *
 * @code
 *	Grid E;
 *	...
 *  int p = j + k*E.nNodesProd[1] + l*E.nNodesProd[2];
 *	double *pVal = &E.value[p*E.nValues];
 *	printf("The field at node (%i,%i,%i) is (%f,%f,%f).\n",j,k,l,pVal[0],pVal[1],pVal[2]);
 * @endcode
 */
/*typedef struct{
	double *val;		///< The values on the grid
	int *nNodes;		///< The number of nodes per direction (nDim elements)
	int *nNodesProd;	///< Cumulative product of nNodes (nDim+1 elements)
	int *compNode;		///< Computational node (nDim elements)
	int *nCompNodes;	///< Number of computational nodes (nDim elements)
	int nDims;			///< Number of dimensions (usually 3)
	int nValues;		///< Number of values per node (usually 1 or 3)
} Grid;*/

/******************************************************************************
 * DEFINED IN IO.C
 *****************************************************************************/

/**
 * Enumeration for message kinds
 * @see msg()
 */
typedef enum{
	STATUS,		///< Normal status output about the progress of execution.
	WARNING,	///< Warning. Something might not be like the user intended.
	ERROR		///< Error which makes the program unable to proceed. Program will stop.
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
 * @see				iniGetStrArr(), iniGetDoubleArr(), iniparser_getint()
 *
 * Remember to free result array using free(). 
 *
 * This function can be seen as an extension to iniparser in order to parse
 * comma-separated entries as arrays.
 */
long int* iniGetIntArr(const dictionary *ini, const char *key, int *nElements);

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

void freePopulation(Population *pop);
Population *allocPopulation(dictionary *ini);

#endif // PINC_H
