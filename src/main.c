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
	msg(STATUS|ONCE,"PINC started.");

	int mpiRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Random Number Generator (RNG)
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	dictionary *ini = iniOpen(argc,argv);
	Population *pop = allocPopulation(ini);
	char *dataPath = "data/";

	/*
	 * TEST ZONE
	 */

	Grid *grid = allocGrid(ini);
	int *nGPoints = grid->nGPoints;
	int *nGPointsProd = grid->nGPointsProd;
	int *nNodes = grid->nNodes;
	int *node = grid->node;
	double *posToNode = grid->posToNode;
	int nDims = grid->nDims;
	msg(STATUS|ONCE,"nDims=%i",nDims);
	msg(STATUS|ONCE,"nGPoints={%i,%i,%i}",nGPoints[0],nGPoints[1],nGPoints[2]);
	msg(STATUS|ONCE,"nGPointsProd={%i,%i,%i,%i}",nGPointsProd[0],nGPointsProd[1],nGPointsProd[2],nGPointsProd[3]);
	msg(STATUS,"node={%i,%i,%i}",node[0],node[1],node[2]);
	msg(STATUS|ONCE,"nNodes={%i,%i,%i}",nNodes[0],nNodes[1],nNodes[2]);
	msg(STATUS|ONCE,"posToNode={%f,%f,%f}",posToNode[0],posToNode[1],posToNode[2]);
	msg(STATUS|ONCE,"1/128=%f",1.0/128.0);

	posUniform(ini,pop,grid,rng);
	gsl_rng_set(rng,mpiRank);
	velMaxwell(ini,pop,rng);

	for(int s=0;s<pop->nSpecies;s++){
		msg(STATUS,"specie %i",s);
		msg(STATUS,"iStart=%i, iStop=%i, particles=%i",pop->iStart[s],pop->iStop[s],pop->iStop[s]-pop->iStart[s]+1);
		for(long int i=pop->iStart[s];i<=pop->iStop[s];i++){
			double *pos=&pop->pos[i*nDims];
			double *vel=&pop->vel[i*nDims];
			msg(STATUS,"particle %i: r={%3.2f,%3.2f,%3.2f}, v={%2.2f,%2.2f,%2.2f}",i,pos[0],pos[1],pos[2],vel[0],vel[1],vel[2]);
		}
	}

	writePopulation(dataPath,pop);

	/*
	 * FINALIZE PINC VARIABLES
	 */
	freePopulation(pop);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	gsl_rng_free(rng);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
