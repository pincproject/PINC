/**
 * @file		spectral.h
 * @brief		Spectral solver
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 */
#ifndef SPECTRAL_H
#define SPECTRAL_H

#include <fftw3.h>

/**
 * @brief Spectral solver
 */
typedef struct {
	fftw_plan fftForward;		///< Forward FFT
	fftw_plan fftInverse;		///< Inverse FFT
	fftw_complex *spectrum;		///< To hold the spectrum
	double *spectralFactor;		///< Multiplicative factor which turns rho into phi
    long int spectralSize;		///< Size of spectrum and spectralFactor
} SpectralSolver;

/**
 * @brief Run mode for spectral solver (currently only testing)
 * @param  ini	Input file
 */
void sMode(dictionary *ini);
funPtr sMode_set(dictionary *ini);

/**
 * @brief Allocates and initializes SpectralSolver
 * @param  ini Input file
 * @param  rho Charge density (source)
 * @param  phi Electric potential (unknown)
 * @return     SpectralSolver
 */
SpectralSolver* sAlloc(const dictionary *ini, const Grid *rho, Grid *phi);

/**
 * @brief Frees SpectralSolver
 * @param solver SpectralSolver
 */
void sFree(SpectralSolver *solver);

/**
 * @brief Solves phi given rho
 * @param solver  SpectralSolver
 * @param rho     Charge density (source)
 * @param phi     Electric potential (unknown)
 * @param mpiInfo MpiInfo
 */
void sSolve(const SpectralSolver *solver, Grid *rho, Grid *phi, const MpiInfo *mpiInfo);

funPtr sSolver_set(dictionary *ini);

#endif // SPECTRAL_H
