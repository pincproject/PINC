#ifndef MGSTRUCTS_H
#define MGSTRUCTS_H


typedef struct Grid{
  int *values;        ///< Values stored in grid nNodes^nDim*nValues
  int *nNodes;        ///< The number of nodes per direction
  int *nNodesProd;    ///< Cumulative product of nNodes
  int nDim;           ///< Number of dimensions
  int nValues;        ///< Number of values per dimension
} Grid;


typedef struct Multigrid{
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

//Initialisers for the grid and multigrid structs
void gridInit(Grid *grid, int *nNodes, int nDim, int nValues);

void multigridInit(Multigrid *multigrid, Grid *fineGrid, int nLevels);
/*
void GaussSeidel(void);*/


#endif // MGSTRUCTS_H
