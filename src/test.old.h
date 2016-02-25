/**
 * @file	    test.h
 * @author	    Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC temp tests
 * @date        09.10.15
 *
 * Test Area for PINC model, a place where some small test are made to see that small
 * parts of the code works as suspected. Mainly here to keep main.c clean and uncluttered
 *
 * 		NOT TO BE INCLUDED IN FINAL PRODUCT
 */


#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "pinc.h"
#include "iniparser.h"
#include "multigrid.h"

/**
 * @brief Tests that grids are accessible and editable from
 * Grid and Multigrid structs
 *
 * @param 	ini				dictionary
 * @param 	Grid	Grid
 * @param 	multigrid 		Multigrid
 *
 * Populates a few grid points from Grid, dumps them in the parsefile,
 * changes them from multigrid, dumps them again. Then they are changed from Grid
 * and dumped again
 *
 */

// void testGridAndMGStructs(dictionary *ini, Grid *Grid, Multigrid *multigrid);

/**
 * @brief Temp manual test of boundary send/recieve function
 */

// void testBoundarySendRecieve(dictionary *ini, Grid *Grid, Multigrid *multigrid);

/*
 * temporary for testing slice function
 */

void testGetSlice(dictionary *ini, Grid *Grid);

void testSwapHalo(dictionary *ini, Grid *Grid, MpiInfo *mpiInfo);

/**
 * Short testing routine of Gauss-Seidel RB algorithm
 */
void testmgGS(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,MpiInfo *mpiInfo);

void testRestriction(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,
 					MpiInfo *mpiInfo);

void testMultigrid(dictionary *ini, Multigrid *multiRho, Multigrid *multiPhi,
					MpiInfo *mpiInfo);

void testDerivatives(dictionary *ini, Grid *rho, Grid *phi);

void dumpWholeGridIndexes(dictionary *ini, Grid *grid);


/**
 * @brief Writes information about the grid structs to a parsefile
 * @param ini 		dictionary of the input file
 * @param grid 		grid struct
 * @param nValues	number of values per grid point.
 * @return void
 *
 * TBD
 * 2 subdomains 2D grid example:
 * @code
	111111	222222
	111111	222222
	111111	222222
	111111	222222
 * @endcode
 */
void dumpWholeGrid(dictionary *ini, Grid *grid);
