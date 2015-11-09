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
#include "multigrid.h" 

int main(int argc, char *argv[]){



	/*
	 * INITIALIZE THIRD PARTY LIBRARIES
	 */
	MPI_Init(&argc,&argv);
	msg(STATUS,"PINC started.");

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	 
	dictionary *ini = iniOpen(argc,argv);

	Grid *grid = allocGrid(ini);
	GridQuantity *gridQuantity = allocGridQuantity(ini, grid, 1);
	Multigrid *multigrid = allocMultigrid(ini, gridQuantity);


	/*
	 * 	Test Area
	 */
	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}

	//Writing information to parsedump (debugging)
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);

	//Changing them from the multigrid struct
	for(int i = 0; i < 5; i++){
		multigrid->gridQuantities[0]->val[i] = 5.;
	}

	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES FROM MULTIGRID \n");

	//Check that both are changed
	//Writing information to parsedump (debugging)
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);


	fMsg(ini, "parsedump", "\n CHANGING GRID VALUES BACK FROM GRIDQUANTITY \n");

	for(int i = 0; i < 5; i++){
		gridQuantity->val[i] = 1.;
	}
	gridParseDump(ini, grid, gridQuantity);
	multigridParseDump(ini, multigrid);


	/*
	 *		TESTING DONE
	 *		Results: accessible from both structs and cahnging them works for both places
	 */



	/*
	 * FINALIZE PINC VARIABLES
	 */
	freeGrid(grid);
	freeGridQuantity(gridQuantity);
	iniparser_freedict(ini);

	/*
	 * FINALIZE THIRD PARTY LIBRARIES
	 */
	msg(STATUS,"PINC completed successfully!"); // Needs MPI
	MPI_Finalize();

	return 0;
}
