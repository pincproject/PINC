/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@fys.uio.no>
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
#include <hdf5.h>

#include "pinc.h"
#include "iniparser.h"
#include "pusher.h"
// #include "multigrid.h"
// #include "test.h"

int main(int argc, char *argv[]){


	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	msg(STATUS|ONCE,"PINC started.");
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * INITIALIZE PINC VARIABLES
	 */

	dictionary *ini = iniOpen(argc,argv);
	Population *pop = pAlloc(ini);
	Grid *grid = gAlloc(ini,2);
	MpiInfo *mpiInfo = gAllocMpi(ini);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);

	/*
	 *	TEST AREA
	 */


	pPosUniform(ini,pop,mpiInfo,rng);

	Grid *rho = gAlloc(ini,1);

	Grid *E = gAlloc(ini,3);
	double val[] = {1,0,0};
	gSet(E,val);
	gNormalizeE(ini,E);

	// double *vel = pop->vel;
	double velc[] = {0,0,0};
	pVelSet(pop,velc);

	double *pos = pop->pos;
	pos[0] = 100;
	pos[1] = 0;
	pos[2] = 0;
	pos[30] = 100;
	pos[31] = 0;
	pos[32] = 0;
	pos[60] = 100;
	pos[61] = 0;
	pos[62] = 0;

	gMul(E,0.5);
	puAcc3D1(pop,E);
	gMul(E,2);

	int N = 10;
	for(int n=0;n<N;n++){
		// printf("pos(%i)=(%6.2f,%6.2f,%6.2f) vel(%.1f)=(%6.2f,%6.2f,%6.2f)\n",n,pos[0],pos[1],pos[2],n+0.5,vel[0],vel[1],vel[2]);
		printf("1: (%6.2f,%6.2f,%6.2f), 2: (%6.2f,%6.2f,%6.2f), 3: (%6.2f,%6.2f,%6.2f)\n",
				pos[0],pos[1],pos[2],pos[30],pos[31],pos[32],pos[60],pos[61],pos[62]);
		puMove(pop);
		puBndPeriodic(pop,E);
		puDistr3D1(pop,rho);
		puAcc3D1(pop,E);
	}


	/*
	 * FINALIZE PINC VARIABLES
	 */

	gFreeMpi(mpiInfo);
	pFree(pop);
	gFree(grid);
	iniparser_freedict(ini);
	gsl_rng_free(rng);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	msg(STATUS|ONCE,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
