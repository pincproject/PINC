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
	Grid *grid = allocGrid(ini);

	/*
	 * TEST ZONE
	 */
//	Timer *timer = allocTimer(0);
	posUniform(ini,pop,grid,rng);
//	tMsg(timer,"position");
	gsl_rng_set(rng,mpiRank);
	velMaxwell(ini,pop,rng);
//	tMsg(timer,"velocity");
//	free(timer);

	Timer *timer = allocTimer(0);
	for(int i=0;i<1000;i++){
		writePopulation("test",pop);
	}
	tMsg(timer,"writePopulation times 1000");

	freeTimer(timer);


	/*
	 * FINALIZE PINC VARIABLES
	 */
//	freePopulation(pop);
//	freeGrid(grid);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	gsl_rng_free(rng);

	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
