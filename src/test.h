/**
 * @file	    test.h
 * @author	    Gullik Vetvik Killie <gullikvk@fys.uio.no>
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
 * GridQuantity and Multigrid structs
 * 
 * @param 	ini				dictionary
 * @param 	gridQuantity	GridQuantity
 * @param 	multigrid 		Multigrid
 *
 * Populates a few grid points from gridQuantity, dumps them in the parsefile, 
 * changes them from multigrid, dumps them again. Then they are changed from gridQuantity
 * and dumped again
 *
 */

void testGridAndMGStructs(dictionary *ini, GridQuantity *gridQuantity, Multigrid *multigrid);

/**
 * @brief Temp manual test of boundary send/recieve function
 */

void testBoundarySendRecieve(dictionary *ini, GridQuantity *gridQuantity, Multigrid *multigrid);

/**
 * @brief Dumps a 2D grid to parsefile
 */
void dump2DGrid(dictionary *ini, GridQuantity *gridQuantity);
