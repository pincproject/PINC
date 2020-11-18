/**
 * @file		hypre.h
 * @brief		Hypre multigrid solver
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */
#ifndef HYPRE_H
#define HYPRE_H

/**
 * @brief Hypre multigrid solver
 */
typedef struct {
} HypreSolver;

/**
 * @brief Allocates and initializes HypreSolver
 * @param  ini Input file
 * @param  rho Charge density (source)
 * @param  phi Electric potential (unknown)
 * @return     HypreSolver
 */
HypreSolver* hyAlloc(const dictionary *ini, const Grid *rho, Grid *phi);

/**
 * @brief Frees HypreSolver
 * @param solver HypreSolver
 */
void hyFree(HypreSolver *solver);

/**
 * @brief Solves phi given rho
 * @param solver  HypreSolver
 * @param rho     Charge density (source)
 * @param phi     Electric potential (unknown)
 * @param mpiInfo MpiInfo
 */
void hySolve(const SpectralSolver *solver, Grid *rho, Grid *phi, const MpiInfo *mpiInfo);

funPtr hySolver_set(dictionary *ini);

#endif // HYPRE_H
