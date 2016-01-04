#ifndef MULTIGRID_H
#define MULTIGRID_H

/**
 * @file		multigrid.h
 * @author		Gullik Vetvik Killie <gullikvk@student.matnat.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Poisson Solver, multigrid.
 * @date		26.10.15
 *
 *
 * Functions dealing with the initialisation and destruction of multigrid structures and
 * a multigrid solver containing restriction, prolongation operatorors and smoothers
 *
 */


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
    Grid **grids;   ///< Array of Grid structs of decreasing coarseness
    int nLevels;         			///< #Grid levels
    int nMGCycles;         			///< Multigrid cycles we want to run
	int nPreSmooth;					///<
	int nPostSmooth;
	int nCoarseSolve;

    void (*coarseSolv)(Grid *phi, const Grid *rho, const int nCycles);	///< Function pointer to a Coarse Grid Solver function
    void (*postSmooth)(Grid *phi, const Grid *rho, const int nCycles);	///< Function pointer to a Post Smooth function
    void (*preSmooth)(Grid *phi, const Grid *rho, const int nCycles);	///< Function pointer to a Pre Smooth function
	void (*restrictor)(const Grid *fine, Grid *coarse);	///< Function pointer to restrictor
	void (*prolongator)(Grid *fine, const Grid *coarse, const MpiInfo *mpiInfo);	///< Function pointer to prolongator
} Multigrid;

Multigrid *mgAlloc(const dictionary *ini, Grid *grid);

 /**
  * @brief Free multigrid struct, top gridQuantity needs to be freed seperately
  * @param Multigrid *multigrid
  *
  * Since the finest grid is allocated seperately and used on it's own without
  * the multigrid struct, it is not freed in this destructor.
  * Variables freed: gridQuantity [1->end]
  *
  */
void mgFree(Multigrid *multigrid);

void linearMGSolv(Multigrid *multiRho, Multigrid *multiPhi);

void gaussSeidel2D(Grid *phi, const Grid *rho, const int nCycles);
void jacobian(Grid *phi, const Grid *rho, const int nCycles);

//Restriction and prolongators
void halfWeightRestrict2D(const Grid *fine, Grid *coarse);
void halfWeightRestrict3D(const Grid *fine, Grid *coarse);
void bilinearProlong2D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);
void bilinearProlong3D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);


 #endif // MULTIGRID_H
