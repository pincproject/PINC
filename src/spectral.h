/**
 * @file		spectral.h
 * @brief		Spectral solver
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */
#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <fftw3.h>

typedef struct {

	fftw_plan fftForward, fftInverse;
	fftw_complex *spectrum;
	double *spectralFactor;
    long int spectralSize;

} SpectralSolver;

funPtr sMode_set(dictionary *ini);
void sMode(dictionary *ini);
SpectralSolver* sAlloc(const dictionary *ini, const Grid *rho, Grid *phi);
void sFree(SpectralSolver *solver);
funPtr sSolve_set(dictionary *ini);
void sSolve(const SpectralSolver *solver, Grid *rho, Grid *phi);

#endif // SPECTRAL_H
