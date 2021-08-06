/**
 * @file		hyprepoisson.h
 * @brief		hypre multigrid Solver for poisson problem.
 * @author		Steffen Brask <s.m.brask@fys.uio.no>
 *
 * Functions dealing with the ..
 */


#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
//#include <mpi.h>
#include "core.h"
#include "hyprepoisson.h"

/******************************************************************************
 * 				Local functions
 *****************************************************************************/



/*************************************************
 *		Inline functions
 ************************************************/





/*************************************************
 * 		ALLOCATIONS
 * 		DESTRUCTORS
 ************************************************/
 HypreSolver* hAlloc(const dictionary *ini, const Grid *rho, Grid *phi){

 	HypreSolver *solver = (HypreSolver *)malloc(sizeof(*solver));

 	// Alloc rest of struct

 	return solver;
 }

 void hFree(HypreSolver *solver){
  //free hypre stuff
 }

 /**
  * @brief Returns routines necessary for using solver
  * @param[out] solve			&hSolve()
  * @param[out] solverAlloc	&hAlloc()
  * @param[out] solverFree	&hFree()
  */
 void hSolver(	void (**solve)(),
 				HypreSolver *(**solverAlloc)(),
 				void (**solverFree)()){

 	*solve=hSolve;
 	*solverAlloc=hAlloc;
 	*solverFree=hFree;
 }

 funPtr hSolver_set(dictionary *ini){

  // Sanity checks:
 	//int nDims = iniGetInt(ini,"grid:nDims");
 	//if(nDims!=1) msg(ERROR,"sMode only works with grid:nDims=1");

 	return hSolver;
 }

 void hSolve(const HypreSolver *solver,
 	Grid *rho, Grid *phi){

 	int rank = rho->rank;
 	int *nGhostLayers = (int *)malloc(2*rank*sizeof(*nGhostLayers));
 	memcpy(nGhostLayers, rho->nGhostLayers, 2*rank*sizeof(*nGhostLayers));

 	gRemoveHalo(rho);
 	gRemoveHalo(phi);

  // rest of Hypre Solv here ...

 	gInsertHalo(rho, nGhostLayers);
 	gInsertHalo(phi, nGhostLayers);
 }
