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

#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"

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

	int nDims = 3;

	Grid *grid = malloc(sizeof(Grid));
	grid->nDims = nDims;
	grid->node = malloc(nDims*sizeof(int));
	grid->node[0] = 1;
	grid->node[1] = 1;
	grid->node[2] = 1;
	grid->nNodes = malloc(nDims*sizeof(int));
	grid->nNodes[0] = 4;
	grid->nNodes[1] = 4;
	grid->nNodes[2] = 4;
	grid->posToNode = malloc(nDims*sizeof(double));
	grid->posToNode[0] = 1.0/128;
	grid->posToNode[1] = 1.0/128;
	grid->posToNode[2] = 1.0/128;

	populateUniformly(ini,pop,grid,rng);

	for(int s=0;s<pop->nSpecies;s++){
		msg(STATUS,"specie %i",s);
		msg(STATUS,"iStart=%i, iStop=%i, particles=%i",pop->iStart[s],pop->iStop[s],pop->iStop[s]-pop->iStart[s]+1);
		for(long int i=pop->iStart[s];i<=pop->iStop[s];i++){
			double *pos=&pop->pos[i*nDims];
//			msg(STATUS,"particle at %f,%f,%f",pos[0],pos[1],pos[2]);
		}
	}


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


