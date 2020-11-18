/**
 * @file		hypre.c
 * @brief		Hypre multigrid solver
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 */
#define _XOPEN_SOURCE 700

#include "core.h"
#include "spectral.h"

HypreSolver* hyAlloc(const dictionary *ini, const Grid *rho, Grid *phi){

	HypreSolver *solver = (HypreSolver *)malloc(sizeof(*solver));

	return solver;
}

void hyFree(HypreSolver *solver){

	free(solver);
}

/**
 * @brief Returns routines necessary for using solver
 * @param[out] solve		&hySolve()
 * @param[out] solverAlloc	&hyAlloc()
 * @param[out] solverFree	&hyFree()
 */
void hySolver(	void (**solve)(),
				HypreSolver *(**solverAlloc)(),
				void (**solverFree)()){

	*solve=hySolve;
	*solverAlloc=hyAlloc;
	*solverFree=hyFree;
}

funPtr sSolver_set(dictionary *ini){
	return hySolver;
}

void hySolve(const HypreSolver *solver,
	Grid *rho, Grid *phi, const MpiInfo *mpiInfo){

	gRemoveHalo(rho);
	gRemoveHalo(phi);

	gInsertHalo(rho, nGhostLayers);
	gInsertHalo(phi, nGhostLayers);
}
