/**
 * @file		units.h
 * @brief		Handles units and normalization.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *
 * Core module for handling of units and normalization.
 */

void parseIndirectInput(dictionary *ini);
Units *uSemiSI(dictionary *ini);
Units *uSI(dictionary *ini);
void uFree(Units *units);
void uNormalize(dictionary *ini, const Units *units);
