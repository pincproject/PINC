
#ifndef MULTIGRID_H
#define MULTIGRID_H

/**
 * @file    multigrid.h
 * @author    Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright University of Oslo, Norway
 * @brief   Poisson Solver, Multigrid.
 * @date    29.10.15
 *
 * Functions dealing with the initialisation and destruction of multigrid structures and
 * a multigrid solver containing restriction, prolongation operatorors and smoothers 
 * 
 */

#include "pinc.h"
#include "iniparser.h"



/**
 * @brief Contains the grids needed in a multigrid solver, as well as other specifications
 * 
 * 
 *
 *
 */
typedef struct {
   Grid *grids;         		///< Array of Grid structs of decreasing coarseness
   int nLevels;         		///< #Grid levels
   int nCycles;         		///< Multigrid cycles we want to run

   void (*coarseSolv)(void);	///< Function pointer to a Coarse Grid Solver function
   void (*postSmooth)(void);	///< Function pointer to a Post Smooth function
   void (*preSmooth)(void);		///< Function pointer to a Pre Smooth function
} Multigrid;


//Initialisers for the grid and multigrid structs
Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity);

//Dumps the 

void gaussSeidel(void);
void jacobian(void);

#endif // MULTIGRID_H
