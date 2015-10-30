
#ifndef MULTIGRID_H
#define MULTIGRID_H

/**
 * @file    multigrid.c
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

/*
 * Not needed anymore, grid is in pinc.h
 * Storing it here for temporary convenience
*/

typedef struct{
  double *val;    ///< The values on the grid
  int *nGPoints;    ///< The number of nodes per direction (nDim elements)
  int *nGPointProd;  ///< Cumulative product of nNodes (nDim+1 elements)
  int *node;    ///< Computational node (nDim elements)
  int *nNodes;  ///< Number of computational nodes (nDim elements)
  int nDims;      ///< Number of dimensions (usually 3)
  int nValues;    ///< Number of values per node (usually 1 or 3)
} Grid;


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
/*   void (*preSmooth)(void);
*/} Multigrid;

//Useful functions
int *getCompNode(const dictionary *ini);

//Initialisers for the grid and multigrid structs
Grid *allocGrid(const dictionary *ini,const int nValues);

void freeGrid(Grid *grid);

void allocMultigrid(void);
/*
void GaussSeidel(void);*/


#endif // MULTIGRID_H
