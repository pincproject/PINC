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
#include <hdf5.h>

int main(int argc, char *argv[]){


	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	msg(STATUS|ONCE,"PINC started.");
	MPI_Barrier(MPI_COMM_WORLD);

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
	gsl_rng_set(rng,grid->mpiRank);
	velMaxwell(ini,pop,rng);
//	tMsg(timer,"velocity");
//	free(timer);

//	int nDims = pop->nDims;
//	for(long int i=0;i<1000;i++){
//		for(int d=0;d<nDims;d++){
//			pop->pos[i*nDims+d]= i + (double)d/10;
//		}
//	}
	posDebug(ini,pop);

/*
	Timer *timer = allocTimer(0);
	createPopulationH5(ini,pop);
	tMsg(timer,"allocate");
	for(int n=0;n<3;n++){
		writePopulationH5(pop,n,n+0.5);
	}
	tMsg(timer,"method 1");
	for(int n=0;n<3;n++){
		writePopulationH52(pop,n,n+0.5,mpiInfo);
	}
	tMsg(timer,"method 2");
	H5Fclose(pop->h5);
	freeTimer(timer);
*/

	/*
	 * FINALIZE PINC VARIABLES
	 */
	freePopulation(pop);
	freeGrid(grid);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	gsl_rng_free(rng);

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
