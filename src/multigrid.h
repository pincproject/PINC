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
  * options are Gauss-Seidel Red-Black (mgGS). More TBD.
  */
 typedef struct {
    Grid **grids;   ///< Array of Grid structs of decreasing coarseness
    int nLevels;         			///< #Grid levels
    int nMGCycles;         			///< Multigrid cycles we want to run
	int nPreSmooth;					///<
	int nPostSmooth;
	int nCoarseSolve;

    void (*coarseSolv)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Coarse Grid Solver function
    void (*postSmooth)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Post Smooth function
    void (*preSmooth)(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);	///< Function pointer to a Pre Smooth function
	void (*restrictor)(const Grid *fine, Grid *coarse);	///< Function pointer to restrictor
	void (*prolongator)(Grid *fine, const Grid *coarse, const MpiInfo *mpiInfo);	///< Function pointer to prolongator

} Multigrid;

/**
 * @brief	Function pointers for the different slice operations
 * @see gHaloOp
 */
typedef void (*MgAlgo)(int level,int bottom, int top, Multigrid *mgRho, Multigrid *mgPhi,
									Multigrid *mgRes, const MpiInfo *mpiInfo);



/**
 * @brief Allocates multigrid struct
 * @param grid 		Finest grid
 * @param ini		Dictionary
 *
 *	Allocates memory for a multgrid struct. A multigrid struct consist of an array of
 *	grids, where each subsequent grid has half the true grid points of the previous.
 *
 *	In addition it stores a few additional parameters that are useful during the
 *	handling of the multigrid algorithm.
 * 	nLevels:	Number of grids
 *	nMGCycles:	Number of MG cycles to run
 *	nPreSmooth:	Number of cycles the presmoother to run
 *	nPostSmooth:Number of cycles the postsmoother to run
 *	nCoarseSolve:Number of cycles for the coarse solver to run
 *
 *	The algorithms for the solver, restrictors and prolongators are set in the allocation according
 *	to a input file, then it is handled by a function pointer stored in the MG struct.
 *'	So if we want to use the prolongator, the function pointer can just be used as a regular function
 * 	and the chosen prolongator algorithm will be used.
 * 	\code
		multigrid->prolongator(fine, coarse, mpiInfo);
 *	\endcode
 *
 *
 *	NB!The number of true grid points used in the finest grid needs to be a multiple
 *	nLevels*2, to make it possible to half the grid points down to the coarsest grid.
 */

Multigrid *mgAlloc(const dictionary *ini, Grid *grid);

 /**
  * @brief Free multigrid struct, top gridQuantity needs to be freed seperately
  * @param 	multigrid
  *
  * Since the finest grid is allocated seperately and used on it's own without
  * the multigrid struct, it is not freed in this destructor.
  * Variables freed: gridQuantity [1->end]
  *
  */
void mgFree(Multigrid *multigrid);

MgAlgo getMgAlgo(const dictionary *ini);


/**
 * @brief Performs a Multigrid run on a test case, used to optimize
 * @param 	ini
 *
 *  This run performs Multigrid runs and times how long it takes to
 *  reach a certain error of the solution. The time spent is saved as a history
 *  file. It is mainly used in the multigrid parameter optimizer, mgOptimize.py,
 *  found in the script folder.
 *
 */
void mgRun(dictionary *ini);

/**
 * @brief Performs a Multigrid run on a test case, used to optimize
 * @param 	ini
 *
 * This run solves sinusoidal test case, then it compares the solution
 * to a analytical solution. It is used by a framework routine to
 * investigate the scaling of the error compared to the stepsize.
 *
 */

void mgErrorScaling(dictionary *ini);


/**
 * @brief Performs a multigrid V cycle
 * @param   level       Grid level the V cycle starts on
 * @param   bottom      Grid level at the bottom of the cycle
 * @param   top         Grid level at the top of the cycle
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mpiInfo     MpiInfo struct containing subdomain information
 */
void mgVRegular(int level,int bottom, int top, Multigrid *mgRho, Multigrid *mgPhi,
 									Multigrid *mgRes, const MpiInfo *mpiInfo);
/**
 * @brief Performs a recursive multigrid V cycle
 * @param   level       Grid level the V cycle starts on
 * @param   bottom      Grid level at the bottom of the cycle
 * @param   top         Grid level at the top of the cycle
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mpiInfo     MpiInfo struct containing subdomain information
 */
void mgVRecursive(int level, int bottom, int top, Multigrid *mgRho, Multigrid *mgPhi,
 					Multigrid *mgRes, const MpiInfo *mpiInfo);

/**
 * @brief Performs a Full multigrid cycle
 * @param   level       Grid level the V cycle starts on
 * @param   bottom      Grid level at the bottom of the cycle
 * @param   top         Grid level at the top of the cycle
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mpiInfo     MpiInfo struct containing subdomain information
 */
void mgFMG(int level, int bottom, int top, Multigrid *mgRho, Multigrid *mgPhi,
 					Multigrid *mgRes, const MpiInfo *mpiInfo);
/**
 * @brief Performs a multigrid W cycle
 * @param   level       Grid level the V cycle starts on
 * @param   bottom      Grid level at the bottom of the cycle
 * @param   top         Grid level at the top of the cycle
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mgRho       MgGrid struct containing rho
 * @param   mpiInfo     MpiInfo struct containing subdomain information
 */
void mgW(int level, int bottom, int top, Multigrid *mgRho, Multigrid *mgPhi,
 					Multigrid *mgRes, const MpiInfo *mpiInfo);



/**
 * @brief Solves Poissons equation for electric potential, with multigrid V cycles
 * @param	mrRho	Source term
 * @param	mgPhi	Solution term
 * @param	mgRes	Residual
 * @param	mpiInfo	Subdomain information
 * @return	mgPhi
 *
 *	This is an implementation of a Multigrid V Cycle solver. See "DOC" for more information.
 */

void mgSolve(MgAlgo mgAlgo, Multigrid *mgRho, Multigrid *mgPhi, Multigrid *mgRes, const MpiInfo *mpiInfo);

funPtr mgSolve_set(dictionary *ini);

/**
 * @brief Gauss-Seidel Red and Black 3D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	3D dimensional implementation of Gauss-Seidel RB, which does several sweeps through the
 *	grid trying to simplify the iteration through the grid.
 *
 */
void mgGS3DNew(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief Gauss-Seidel Red and Black 3D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	3D dimensional implementation of Gauss-Seidel RB, which does one sweep through the grid for
 *	each color, but has slightly more complicated behaviour on the edges, due to needing to skip
 *	the ghostlayers.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void mgGS3D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief Gauss-Seidel Red and Black 2D
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 *	2D dimensional implementation of Gauss-Seidel RB, which does several sweeps through the
 *	grid trying to simplify the iteration through the grid.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void mgGS2D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);


void mgGSND(Grid *phi, const Grid *rho, int nCycles, const MpiInfo *mpiInfo);

/**
 * @brief mgJacob method
 * @param	rho		Source term
 * @param	phi		Solution term
 * @param	mpiInfo	Subdomain information
 * @return	phi
 *
 * Non-optimized implementation of a mgJacob2D algorithm to solve poissons equation.
 *
 *	NB! Assumes 1 ghost layer, and even number of grid points.
 */
void mgJacobND(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);
void mgJacob1D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);
void mgJacob3D(Grid *phi, const Grid *rho, const int nCycles, const MpiInfo *mpiInfo);


/**
 * @brief Half weight restriction, 2D
 * @param	fine	Source term
 * @param	coarse	Solution term
 * @param	mpiInfo	Subdomain information
 * @return	coarse
 *
 *	Implementation of a half weight restriction scheme, copying the fine
 *	down to the coarser grid.
 *
 */

void mgHalfRestrict2D(const Grid *fine, Grid *coarse);


/**
 * @brief Half weight restriction, 3D
 * @param	fine	Source term
 * @param	coarse	Solution term
 * @return	coarse
 *
 *	Implementation of a half weight restriction scheme, copying the fine
 *	down to the coarser grid.
 *
 */
void mgHalfRestrict3D(const Grid *fine, Grid *coarse);

/**
 * @brief Half weight restriction, ND
 * @param	fine	Source term
 * @param	coarse	Solution term
 * @return	coarse
 *
 *	Implementation of a half weight restriction scheme, copying the fine
 *	down to the coarser grid.
 *
 */

void mgHalfRestrictND(const Grid *fine, Grid *coarse);



/**
 * @brief Bilinear interpolation, 2D
 * @param	fine	Fine grid
 * @param	coarse	Coarse grid
 * @param	mpiInfo	Subdomain information
 * @return	fine
 *
 *	Implementation of a bilnear interpolation scheme, interpolating the coarse
 *	grid onto the fine grid.
 *
 */

void mgBilinProl2D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);

/**
 * @brief Bilinear interpolation, 3D
 * @param	fine	Fine grid
 * @param	coarse	Coarse grid
 * @param	mpiInfo	Subdomain information
 * @return	fine
 *
 *	Implementation of a bilnear interpolation scheme, instructterpolating the coarse
 *	grid onto the fine grid.
 *
 */
void mgBilinProl3D(Grid *fine,const Grid *coarse, const MpiInfo *mpiInfo);
/**
 * @brief Bilinear interpolation, ND
 * @param	fine	Fine grid
 * @param	coarse	Coarse grid
 * @param	mpiInfo	Subdomain information
 * @return	fine
 *
 *	Implementation of a bilnear interpolation scheme, instructterpolating the coarse
 *	grid onto the fine grid.
 *
 */

void mgBilinProlND(Grid *fine, const Grid *coarse,const  MpiInfo *mpiInfo);

/*
 * @brief Restrict boundary conditions down to coarser grids
 * @param mgGrid   Multigrid struct
 *
 *  OBS, WARNING!
 *  INJECTION PROBABLY ONLY WORKS FOR CONSTANT DIRICHLET AND NEUMANN CONDITIONS
 */
void mgRestrictBnd(Multigrid *mgGrid);


/**
 * @brief Computes residual
 * @param	res		Residual grid
 * @param	phi		Phi	grid
 * @param	rho		Rho grid
 * @param	mpiInfo	Subdomain information
 * @return	res
 *
 *	Computes the residual on a grid level.
 *	\f[
 *		d_l = \nabla^2_l\phi_l - \rho_l
 *	\f]
 */
void mgResidual(Grid *res, const Grid *rho, const Grid *phi,const MpiInfo *mpiInfo);

/**
 * @brief Returns mass of a grid
 * @param	grid		Grid struct
 * @param	ini			Dictionary
 *
 *	Computes a mass of a grid:
 *	\f[
 *		M = \frac{1}{N} \sum_N d^2_N
 *	\f]
 *
 * NB! Not written for efficiency
 */
 double mgResMass3D(Grid *grid, MpiInfo *mpiInfo);


 /**
  * @brief Compares numerical solution to an analytical solution
  * @param  numerical           Numerical solution
  * @param  analytical          Analytical solution
  * @param  error               Difference between solutions
  * @return error
  */

/**
 * @brief Computes avg error, returns in a percentage
 * @param  numerical           Numerical solution
 * @param  analytical          Analytical solution
 * @param  error               Difference between solutions
 * @return error
 */
void mgCompError(const Grid *numerical,const Grid *analytical, Grid *error);


double	mgAvgError(Grid *phi,Grid *sol,Grid *error,MpiInfo *mpiInfo);

/**
 * @brief Returns the square of the error on the true grid
 * @param  error               Difference between solutions
 * @param  mpiInfo             MpiInfo
 * @return error
 *
 *  WARNING!    Stores the squared values on the original grid
 *              Recompute error if needed
 */
double mgSumTrueSquared(Grid *error,const MpiInfo *mpiInfo);

 /**
  * @brief Writes out information about the MG cycles, used when optimizing the number of cycles
  * @param ini 					dictionary of the input file
  * @param multigrid 		multigrid struct
  *
  * Debug help
  */
 void parseMGOptim(dictionary *ini, Multigrid *multigrid);



 #endif // MULTIGRID_H
