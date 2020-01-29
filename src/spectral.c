/**
 * @file		spectral.c
 * @brief		Spectral solver
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 */
#define _XOPEN_SOURCE 700

#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include "core.h"
#include "spectral.h"

SpectralSolver* sAlloc(const dictionary *ini, const Grid *rho, Grid *phi){

	SpectralSolver *solver = (SpectralSolver *)malloc(sizeof(*solver));

	int nDims = iniGetInt(ini,"grid:nDims");
	long int *trueSize = iniGetLongIntArr(ini,"grid:trueSize",nDims);
	long int size = trueSize[0];
	long int spectralSize = size/2+1;
	solver->spectralSize = spectralSize;
	free(trueSize);

	fftw_complex *spectrum = (fftw_complex *)malloc(size*sizeof(*spectrum));
	double *spectralFactor = (double *)malloc(spectralSize*sizeof(*spectralFactor));

	spectralFactor[0] = 0; // Actually infinity
	for(int n=1; n<spectralSize; n++){

		// Spectral multiplier turning rho into phi
		spectralFactor[n] = size/(2*M_PI*n);
		spectralFactor[n] *= spectralFactor[n];

		// Part of IFFT operations
		spectralFactor[n] /= size;
	}


	spectrum = (fftw_complex *)fftw_malloc(spectralSize*sizeof(fftw_complex));

	// Replacing FFTW_ESTIMATE with FFTW_MEASURE may be beneficial for very
	// large systems. More efficient algorithm by omitting FFTW_PRESERVE_INPUT
	// on c2r-transform.
	solver->fftForward = fftw_plan_dft_r2c_1d(size,rho->val,spectrum,FFTW_ESTIMATE);
	solver->fftInverse = fftw_plan_dft_c2r_1d(size,spectrum,phi->val,FFTW_ESTIMATE||FFTW_PRESERVE_INPUT);

	solver->spectrum = spectrum;
	solver->spectralFactor = spectralFactor;

	return solver;
}

void sFree(SpectralSolver *solver){

	fftw_destroy_plan(solver->fftForward);
	fftw_destroy_plan(solver->fftInverse);
	fftw_cleanup();

	free(solver->spectralFactor);
	free(solver->spectrum);
	free(solver);
}

/**
 * @brief Returns routines necessary for using solver
 * @param[out] solve			&sSolve()
 * @param[out] solverAlloc	&sAlloc()
 * @param[out] solverFree	&sFree()
 */
void sSolver(	void (**solve)(),
				SpectralSolver *(**solverAlloc)(),
				void (**solverFree)()){

	*solve=sSolve;
	*solverAlloc=sAlloc;
	*solverFree=sFree;
}

funPtr sSolver_set(dictionary *ini){

	int nDims = iniGetInt(ini,"grid:nDims");
	if(nDims!=1) msg(ERROR,"sMode only works with grid:nDims=1");

	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", nDims);
	if(nSubdomains[0]!=1) msg(ERROR,"sMode only works with grid:nSubdomains=1");
	free(nSubdomains);

	return sSolver;
}

void sSolve(const SpectralSolver *solver,
	Grid *rho, Grid *phi, const MpiInfo *mpiInfo){

	int rank = rho->rank;
	int *nGhostLayers = (int *)malloc(2*rank*sizeof(*nGhostLayers));
	memcpy(nGhostLayers, rho->nGhostLayers, 2*rank*sizeof(*nGhostLayers));

	gRemoveHalo(rho);
	gRemoveHalo(phi);

	fftw_execute(solver->fftForward);

	// Set DC-component to zero for charge neutrality
	solver->spectrum[0] = 0;

	for(int n=1; n<solver->spectralSize; n++){
		solver->spectrum[n] *= solver->spectralFactor[n];
	}

	fftw_execute(solver->fftInverse);

	gInsertHalo(rho, nGhostLayers);
	gInsertHalo(phi, nGhostLayers);
}

funPtr sMode_set(dictionary *ini){
	int nDims = iniGetInt(ini,"grid:nDims");
	if(nDims!=1) msg(ERROR,"sMode only works with grid:nDims=1");

	int *nSubdomains = iniGetIntArr(ini, "grid:nSubdomains", nDims);
	if(nSubdomains[0]!=1) msg(ERROR,"sMode only works with grid:nSubdomains=1");
	free(nSubdomains);

	return sMode;
}
void sMode(dictionary *ini){


	MpiInfo *mpiInfo = gAllocMpi(ini);
	Grid *phi = gAlloc(ini, SCALAR);
	Grid *rho = gAlloc(ini, SCALAR);

	SpectralSolver *solver = sAlloc(ini, rho, phi);

	int *trueSize = rho->trueSize;
	int *nGhostLayers = rho->nGhostLayers;
	double *rhoValStart = &rho->val[nGhostLayers[1]];
	double *phiValStart = &phi->val[nGhostLayers[1]];
	for(long int j=0; j<trueSize[1]; j++){
		rhoValStart[j] = sin(2*M_PI*j/trueSize[1]);
	}

	adPrint(rhoValStart, trueSize[1]);
	sSolve(solver, rho, phi, mpiInfo);
	adPrint(phiValStart, trueSize[1]);

	sFree(solver);
	gFree(rho);
	gFree(phi);
	gFreeMpi(mpiInfo);

}
