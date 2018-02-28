/**
 * @file		core.h
 * @brief		PINC main header.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 * 				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 *
 * The PINC main header file constitutes a framework of function and struct
 * declarations central to all PINC modules.
 */

#ifndef CORE_H
#define CORE_H

#include "iniparser.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <hdf5.h>
#include <gsl/gsl_rng.h>
#include "version.h"

/******************************************************************************
 * DEFINING PHYSICAL CONSTANTS (following SI standard)
 *****************************************************************************/

/* #define BOLTZMANN 1.3806488e-23 */
/* #define ELECTRON_MASS 9.10938215e-31 */

/******************************************************************************
 * DEFINING CORE DATATYPES (used by several modules)
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
 * The position of the particles is normalized with respect to the step size in
 * Grid, such that a particle with local position (1,2,3) is located _on_ node
 * (1,2,3) in the grid. Particles are usually specified in local frame but may
 * temporarily be expressed in global frame. See MpiInfo.
 *
 * kinEnergy and potEnergy stores the kinetic and potential energy of the
 * particles if an energy-computing function is utilized (see e.g. puAcc3D1KE
 * and gPotEnergy). Some energy-computing functions may be able to compute the
 * energy per specie, in which case kinEnergy[s] and potEnergy[s] is the energy
 * for specie s. Others may only be able to compute the net energy for all
 * species, in which case this is stored in the last element (
 * kinEnergy[nSpecies] or potEnergy[nSpecies]). Note that this is only the
 * energy for the current subdomain, and a separate function must be employed
 * to sum the energy across the subdomains and store it to an .h5-file.
 *
 * If a population h5 output file is created, the handler to this file is
 * stored in h5.
 */
typedef struct{
	double *pos;		///< Position
	double *vel;		///< Velocity
	long int *iStart;	///< First index of specie s (nSpecies+1 elements)
	long int *iStop;	///< First index not of specie s (nSpecies elements)
	double *charge;		///< Normalized charge (q-bar)
	double *mass;		///< Normalized mass (m-bar)
	double *kinEnergy;	///< Kinetic energy (nSpecies+1 elements)
	double *potEnergy;	///< Potential energy (nSpecies+1 elements)
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

	int nSpecies;				///< Number of species
	int nNeighbors;				///< Number of neighbors (3^nDims-1) TBD: Omit if it's faster to recompute each time
	int neighborhoodCenter;		///< Index of center/self in neighborhood
	long int **migrants;		///< nMigrants (DEPRECATED)
	long int **migrantsDummy;	///< Useful in computations (DEPRECATED)
	long int *nEmigrants;		///< Number of migrants of each specie to each neighbor (nSpecies*nNeighbor elements)
	long int *nEmigrantsAlloc;	///< Number of migrants allocated for to each neighbor (nNeighbor elements)
	long int *nImmigrants;		///< Number of immigrants of each specie from each neighbour (nSpecies*nNeighbor elements)
	long int nImmigrantsAlloc;
	double **emigrants;			///< Buffer to house emigrants
	double **emigrantsDummy;	///< YAY
	double *immigrants;			///< Buffer to house immigrants
	double *thresholds;			///< Threshold for migration (2*nDims elements)

	MPI_Request *send;
	MPI_Request *recv;
} MpiInfo;

/**
 * @brief Defines different types of boundary conditons
 * @see gAlloc
 * @see gBnd
 */
typedef enum{
	PERIODIC = 0x01,		///< Periodic boundary conditions.
	DIRICHLET = 0x02,		///< Dirichlet boundary condtions.
	NEUMANN = 0x03,			///< Neumann boundary conditons.
	NONE = 0x10				///< For nValues
} bndType;


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

	double *sendSlice;			///< Slice buffer of the grid sent to other
	double *recvSlice;			///< Slice buffer of the grid sent to other
	double *bndSlice;			///< Slices used by Dirichlet and Neumann boundaries
	hid_t h5;					///< HDF5 file handler
	hid_t h5MemSpace;			///< HDF5 memory space description
	hid_t h5FileSpace;			///< HDF5 file space description

	bndType *bnd;				///< Array storing boundary conditions
} Grid;

/**
 * @brief Contains characteristic scales to be used for normalization and
 * denormalization in PINC.
 *
 * This datatype contains characteristic scales that can be used to normalize or
 * denormalize quantities in PINC. These characteristic scales have physical
 * units (normally SI units) so if for instance a variable has the dimension of
 * "velocity", it is normalize by dividing by units->velocity. Likewise, a
 * simulation variable E for the E-field, can be obtained in physical units by
 * multiplying by units->eField. 
 *
 * The (new) convention in PINC is that all normalization takes place during
 * initialization
 *
 * @code
 * 	Units *units = uAlloc(ini);
 * 	uNormalize(ini, units);
 *
 * 	Grid *rho = gAlloc(ini, SCALAR);
 * 	Population *pop = pAlloc(ini);
 *
 * 	pOpenH5(ini, pop, units, "pop");
 * 	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
 *
 * 	uFree(units);
 * @endcode
 */
typedef struct{

	double nDims;		///< Number of spatial dimensions
	double weight;		///< Number of physical particle per simulation particle

	// Characteristic SI base units (with charge instead of current)
	double charge;			///< Charge
	double mass;			///< Mass
	double length;			///< Length
	double time;			///< Time
	
	// Derived units
	double velocity;		///< Velocity
	double acceleration;	///< Acceleration
	double chargeDensity;	///< Electric charge density
	double potential;		///< Electric potential
	double eField;			///< Electric field
	double bField;			///< Magnetic flux density
	double energy;			///< Energy

} Units;
/**
 * @brief	Timer struct for simple profiling
 * @see allocTimer(), freeTimer(), tMsg()
 *
 *	Simple timer struct to keep track of time.
 *	A simple example of where the time to add 10 to an integer 10 times is
 *	computed:
 *	\code
 Timer *t = tAlloc();

 int k = 0;
 for (int i = 0; i < 10; i++){
	 tStart(t);
	 k += 10;
	 tStop(t);
	 tMsg(t->total, "Hello: ");
 }

 tFree(t);
 *	\endcode
 */
typedef struct{
	unsigned long long int total;		/// Total time
	unsigned long long int start;		/// Previous start time
} Timer;
//
// unsigned long long int getNanoSec();
// void tMsg(int rank, Timer *timer, format....);
// void tStart(...);
// void tStop(...);
// void tic();
// void toc();


/**
 * @brief Defines different types of messages
 * @see msg()
 */
typedef enum{
	STATUS = 0x00,	///< Normal status output about the progress of execution.
	WARNING = 0x01,	///< Warning. Something might not be like the user intended.
	ERROR = 0x02,	///< Error which makes the program unable to proceed. Program will stop.
	TIMER = 0x03,	///< Printing out formatted timing result
	ALL = 0x10		///< Output message from all MPI-nodes. To be bitwise ORed.
} msgKind;

/**
 * @brief Pointer to function returning void. (All parameter lists works).
 */
typedef void (*funPtr)();


/******************************************************************************
 * INCLUDING CORE MODULES
 *****************************************************************************/

#include "population.h"
#include "grid.h"
#include "io.h"
#include "aux.h"
#include "units.h"

#endif // CORE_H
