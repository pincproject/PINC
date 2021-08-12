/**
 * @file		hyprepoisson.h
 * @brief		hypre multigrid Solver for poisson problem.
 * @author		Steffen Brask <s.m.brask@fys.uio.no>
 *
 * Functions dealing with the ..
 */


#ifndef HYPREPOISSON_H
#define HYPREPOISSON_H

#include <HYPRE_struct_ls.h>
//#include <HYPRE_krylov.h>

/**
 * @brief ..
 * @param listParams ..

 *
 * ..
 */

//DEFINE STRUCTS

/**
 * @brief Hypre Multigrid solver
 */
typedef struct {
	//HYPRE stuff here
  HYPRE_StructGrid     grid;
  HYPRE_StructStencil  stencil;
  HYPRE_StructMatrix   A;
  HYPRE_StructVector   b;
  HYPRE_StructVector   x;

  HYPRE_StructSolver   solver;
  HYPRE_StructSolver   precond;

  int *ilower;
  int *iupper;

} HypreSolver;

//DEFINE allocators destructors

/**
 * @brief Allocates and initializes HypreSolver
 * @param  ini Input file
 * @param  rho Charge density (source)
 * @param  phi Electric potential (unknown)
 * @return     HyprelSolver
 */
HypreSolver* hAlloc(const dictionary *ini, const Grid *rho, Grid *phi);

/**
 * @brief Frees HypreSolver
 * @param solver HypreSolver
 */
void hFree(HypreSolver *solver);


/**
 * @brief Solves phi given rho
 * @param solver  SpectralSolver
 * @param rho     Charge density (source)
 * @param phi     Electric potential (unknown)
 */
void hSolve(const HypreSolver *solver, Grid *rho, Grid *phi);

funPtr hSolver_set(dictionary *ini);

//DEFINE run modes, if needed?

//DEFINE additional global functions


 #endif // HYPREPOISSON_H
