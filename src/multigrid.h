
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


typedef struct {
   Grid *grids;         ///< Grids levels in the multigrid. Each grid level has half nNodes of the previous
   int nLevels;         ///< #Grid levels
   int nCycles;         ///< Multigrid cycles we want to run
   //Possible to store functions in an array, so we can define which is used 
   //at initialization?
   //
   //func Presmoother;
   //func PostSmoother;
   //func CoarseGridSolver;
   void (*coarseSolv)(void);
   void (*postSmooth)(void);
   void (*preSmooth)(void);
} Multigrid;


//Initialisers for the grid and multigrid structs
Grid *allocGrid(const dictionary *ini);

GridQuantity *allocGridQuantity(const dictionary *ini, Grid *grid, int nValues);

Multigrid *allocMultigrid(const dictionary *ini, GridQuantity *gridQuantity);


//Dumps the 
void gridParseDump(dictionary *ini, Grid *grid);

void gaussSeidel(void);
void jacobian(void);

void freeGrid(Grid *grid);

#endif // MULTIGRID_H
