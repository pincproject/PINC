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
 * @param ini 			Input file, contains specification for the run
 * @param gridQuantity  Grid with quantities and grid specifications as a memmber
 *
 * The finest grid, grid 0, links to the grid used in the rest of the program, the other grids
 * in the grid array is coarser grids used in the
 * The preSmooth, postSmooth and coarseSolv is set in the input.ini file, for now the only
 * options are Gauss-Seidel Red-Black (gaussSeidel). More TBD.
 */
typedef struct {
   GridQuantity **gridQuantities;   ///< Array of Grid structs of decreasing coarseness
   int nLevels;         			///< #Grid levels
   int nCycles;         			///< Multigrid cycles we want to run

   void (*coarseSolv)(void);	///< Function pointer to a Coarse Grid Solver function
   void (*postSmooth)(void);	///< Function pointer to a Post Smooth function
   void (*preSmooth)(void);		///< Function pointer to a Pre Smooth function
} Multigrid;

//Initialisers for the grid and multigrid structs
Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity);

/**
 * @brief Free multigrid struct, top gridQuantity needs to be freed seperately
 * @param Multigrid *multigrid
 *
 * Since the finest grid is allocated seperately and used on it's own without
 * the multigrid struct, it is not freed in this destructor.
 * Variables freed: gridQuantity [1->end]
 *
 */
void freeMultigrid(Multigrid *multigrid);


void gaussSeidel(void);
void jacobian(void);

#endif // MULTIGRID_H
