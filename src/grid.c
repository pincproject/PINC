/**
 * @file		grid.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Grid-struct handling.
 * @date		30.10.15
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#include "pinc.h"
#include <mpi.h>

/******************************************************************************
 * DECLARATIONS
 *****************************************************************************/
/**
 * @brief Returns the ND-index of this MPI node in the global reference frame
 * @param	ini		input settings
 * @return	The N-dimensional index of this MPI node
 */
static int *getNode(const dictionary *ini);

/******************************************************************************
 * DEFINITIONS
 *****************************************************************************/

static int *getNode(const dictionary *ini){

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Get ini info
	int nDims;
	int *nNodes = (int*)iniGetIntArr(ini,"grid:nNodes",&nDims);

	// Sanity check
	int totalNNodes = intArrProd(nNodes,nDims);
	if(totalNNodes!=size)
		msg(ERROR,"The product of grid:nNodes does not match the number of MPI processes");

	// Determine node
	int *node = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++){

		node[d] = rank % nNodes[d];
		rank /= nNodes[d];

	}

	free(nNodes);
	return node;

}
