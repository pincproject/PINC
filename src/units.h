/**
 * @file		units.h
 * @brief		Handles units and normalization.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *
 * Core module for handling of units and normalization.
 */

#ifndef UNITS_H
#define UNITS_H

/**
 * @brief	Allocates and populates Unit according to ini-file
 * @param	ini		Dictionary to input file
 * @return			Unit
 *
 * The characteristic scales in Units will be populated according to
 * methods:normalization in the input file:
 *
 * - SI:
 *   All physical input parameters are in SI units, and Units is in SI units.
 *
 * - SemiSI:
 *   Similar to SI but for convenience charge is specified in elementary charges
 *   and mass in electron masses. timeStep is specified in in terms of wpe^(-1)
 *   where wpe is the electron plasma frequency. This is just for convenience,
 *   and just during input. These parameters are immediately transformed to SI
 *   units internally, and Units is still in SI units.
 *
 * Remember to call uFree() to free memory.
 */
Units *uAlloc(dictionary *ini);

/**
 * @brief					Frees memory for Units
 * @param[in,out]	units	Pointer to Units to be freed
 */
void uFree(Units *units);

/**
 * @brief					Normalizes input parameters according to Units
 * @param[in,out]	ini		Dictionary to input file
 * @param			units	Units
 *
 * Remember to run this prior to using the ini variable on anything that reads
 * physical values from it. For instance, reading methods:solver prior to
 * running this is not a problem, since it is an algorithmic parameter and not a
 * physical one. Running pAlloc() prior to this function will be catastrophic,
 * since it reads the charge and mass of the species.
 */
void uNormalize(dictionary *ini, const Units *units);

#endif // UNITS_H
