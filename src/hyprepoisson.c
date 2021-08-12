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
#include <HYPRE_struct_ls.h>
//#include <HYPRE_krylov.h>
#include "core.h"
#include "grid.h"
#include "hyprepoisson.h"


/******************************************************************************
 * 				Local functions
 *****************************************************************************/
 static int *getSubdomain(const dictionary *ini){
  // copied from grid.c, this should be global ...
 	// Get MPI info
 	int mpiSize, mpiRank;
 	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);
 	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

 	// Get ini info
 	int nDims = iniGetInt(ini,"grid:nDims");
 	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);


 	// Sanity check
 	int totalNSubdomains = aiProd(nSubdomains,nDims);
 	if(totalNSubdomains!=mpiSize)
 		msg(ERROR,"The product of grid:nSubdomains does not match the number of MPI processes");

 	// Determine subdomain of this MPI node
 	int *subdomain = malloc(nDims*sizeof(*subdomain));
 	for(int d=0;d<nDims;d++){
 		subdomain[d] = mpiRank % nSubdomains[d];
 		mpiRank /= nSubdomains[d];
 	}

 	free(nSubdomains);
 	return subdomain;

 }


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
   HYPRE_Init();

   int nDims = iniGetInt(ini,"grid:nDims");
   //int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
   int *trueSize = iniGetIntArr(ini, "grid:trueSize", nDims);
   int *subdomain = getSubdomain(ini);

   // SET UP GRID

   int *ilower = malloc(nDims*sizeof(*ilower));// Needs indices of global lower corner for current subdim
   int *iupper = malloc(nDims*sizeof(*iupper));// Needs indices of global upper corner for current subdim

   //exit(0);
   for(int i = 0; i < nDims; i++) ilower[i]=subdomain[i]*trueSize[i];
   for(int i = 0; i < nDims; i++) iupper[i]=(subdomain[i]+1)*trueSize[i]-1;

   //DEBUGING
   //for(int i = 0; i < nDims; i++) printf("%i \n",ilower[i]);
   //for(int i = 0; i < nDims; i++) printf("%i \n",iupper[i]);
   //exit(0);

   HYPRE_StructGridCreate(MPI_COMM_WORLD, nDims, &solver->grid);
   HYPRE_StructGridSetExtents(solver->grid, ilower, iupper);
   HYPRE_StructGridAssemble(solver->grid);

   // DEFINE STENCIL
   int stencilSize = 1+nDims*2;
   HYPRE_StructStencilCreate(nDims, stencilSize, &solver->stencil);

   // Maybe a 3-D stencil works for 2, and 1-D also. The entries are the same ..
   int offsets[7][3] = {{0,0,0}, {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1},{0,0,1}};

   for (int i = 0; i < stencilSize; i++)
     HYPRE_StructStencilSetElement(solver->stencil, i, offsets[i]);

   // SET UP MATRIX (ignoring boundary conditions)

   HYPRE_StructMatrixCreate(MPI_COMM_WORLD, solver->grid, solver->stencil, &solver->A);
   HYPRE_StructMatrixInitialize(solver->A);

   int stencil_indices[7] = {0,1,2,3,4,5,6}; /* labels for the stencil entries -
   these correspond to the offsets
   defined above */
   int nvalues  = phi->sizeProd[nDims+1]*stencilSize;
   double *values = malloc(nvalues*sizeof(*values));
   //double *exact_values = malloc(nvalues*sizeof(*values));

   double stencilCoeff = nDims*2;// 2,4,6, Is this correct?
   for (long int i = 0; i < nvalues; i += stencilSize){
     values[i] = stencilCoeff;
     for (int j = 1; j < stencilSize; j++) values[i+j] = -1.0;
   }

   HYPRE_StructMatrixSetBoxValues(solver->A, ilower, iupper, stencilSize,
     stencil_indices, values);

  /* // SET UP BOUNDARY CONDITIONS */
  //
  // For now we only set periodic. This should be moved to separate function
  //
  int* nPeriod = malloc(nDims*sizeof(*nPeriod)); //periodicity for each dimension
  int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
  for (int i=0;i<nDims;i++) nPeriod[i]=nSubdomains[i]*trueSize[i];
  //HYPRE_StructGridSetPeriodic(solver->grid, nPeriod);
  free(nPeriod);

  //

  HYPRE_StructMatrixAssemble(solver->A);

  // INITIALIZE VECTORS

  HYPRE_StructVectorCreate(MPI_COMM_WORLD, solver->grid, &solver->b);
  HYPRE_StructVectorCreate(MPI_COMM_WORLD, solver->grid, &solver->x);
  HYPRE_StructVectorInitialize(solver->b);
  HYPRE_StructVectorInitialize(solver->x);


  // Create solver
  HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver->solver);

  HYPRE_GMRESSetMaxIter((HYPRE_Solver) solver->solver, 500 );
  HYPRE_GMRESSetKDim((HYPRE_Solver) solver->solver,30);
  HYPRE_GMRESSetTol((HYPRE_Solver) solver->solver, 1.0e-06 );
  HYPRE_GMRESSetPrintLevel((HYPRE_Solver) solver->solver, 0 );
  HYPRE_GMRESSetLogging((HYPRE_Solver) solver->solver, 0 );

  // Symmetric PFMG preconditioner
  int n_pre     = 1;
    int n_post    = 1;
    int rap       = 0;
    int relax     = 1;
    int skip      = 0;
    HYPRE_StructPFMGCreate(MPI_COMM_WORLD, &solver->precond);
    HYPRE_StructPFMGSetMaxIter(solver->precond, 1);
    HYPRE_StructPFMGSetTol(solver->precond, 0.0);
    HYPRE_StructPFMGSetZeroGuess(solver->precond);
    HYPRE_StructPFMGSetRAPType(solver->precond, rap);
    HYPRE_StructPFMGSetRelaxType(solver->precond, relax);
    HYPRE_StructPFMGSetNumPreRelax(solver->precond, n_pre);
    HYPRE_StructPFMGSetNumPostRelax(solver->precond, n_post);
    HYPRE_StructPFMGSetSkipRelax(solver->precond, skip);
    HYPRE_StructPFMGSetPrintLevel(solver->precond, 0);
    HYPRE_StructPFMGSetLogging(solver->precond, 0);
    HYPRE_StructGMRESSetPrecond( solver->solver,
                                 HYPRE_StructPFMGSolve,
                                 HYPRE_StructPFMGSetup,
                                 solver->precond);

  HYPRE_StructGMRESSetup(solver->solver, solver->A, solver->b, solver->x);
  //Additional info allocations
  solver->ilower=ilower;
  solver->iupper=iupper;
  free(values);
  free(trueSize);
  free(subdomain);


  return solver;
}

 void hFree(HypreSolver *solver){
  //free hypre stuff
  HYPRE_StructGMRESDestroy(solver->solver);
  HYPRE_StructPFMGDestroy(solver->precond);
  free(solver->ilower);
  free(solver->iupper);
  free(solver);
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
  //move this to a struct  ^ ^ ^
 	memcpy(nGhostLayers, rho->nGhostLayers, 2*rank*sizeof(*nGhostLayers));

 	gRemoveHalo(rho);
 	gRemoveHalo(phi);

  HYPRE_StructVectorSetBoxValues(solver->b, solver->ilower, solver->iupper, rho->val);

  HYPRE_StructVectorAssemble(solver->b);
  HYPRE_StructVectorAssemble(solver->x);

  // Solve
  HYPRE_StructGMRESSolve(solver->solver, solver->A, solver->b, solver->x);

  HYPRE_StructVectorGetBoxValues(solver->x, solver->ilower, solver->iupper, phi->val);

 	gInsertHalo(rho, nGhostLayers);
 	gInsertHalo(phi, nGhostLayers);
  
  free(nGhostLayers);

 }
