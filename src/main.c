/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 * Replaces old DiP3D main.c file by Wojciech Jacek Miloch.
 */

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include <gsl/gsl_rng.h>

int main(int argc, char *argv[]){

	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	msg(STATUS,"PINC started.");

	// Random Number Generator (RNG)
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	dictionary *ini = iniOpen(argc,argv);
	Population *pop = allocPopulation(ini);

	/*
	 * TEST ZONE
	 */

	Grid *grid = malloc(sizeof(Grid));
	grid->nDims = 3;
	grid->node = malloc(grid->nDims*sizeof(int));
	grid->node[0] = 1;
	grid->node[1] = 1;
	grid->node[2] = 1;
	grid->nNodes = malloc(grid->nDims*sizeof(int));
	grid->nNodes[0] = 1;
	grid->nNodes[1] = 1;
	grid->nNodes[2] = 1;
	grid->nGPoints = malloc(grid->nDims*sizeof(int));
	grid->nGPoints[0] = 16;
	grid->nGPoints[1] = 16;
	grid->nGPoints[2] = 16;
	grid->nGhosts = malloc(grid->nDims*sizeof(int));
	grid->nGhosts[0] = 0;
	grid->nGhosts[1] = 0;
	grid->nGhosts[2] = 0;
	grid->posToNode = malloc(grid->nDims*sizeof(double));
	grid->posToNode[0] = 1.0/16;
	grid->posToNode[1] = 1.0/16;
	grid->posToNode[2] = 1.0/16;

	populateUniformly(ini,pop,grid,rng);


	/*
	 * FINALIZE PINC VARIABLES
	 */
	freePopulation(pop);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	gsl_rng_free(rng);
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}


