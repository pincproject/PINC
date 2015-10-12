/**
 * @file		aux.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC auxiliary functions.
 * @date		11.10.15
 *
 * Small auxiliary functions used throughout the PINC program.
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "pinc.h"

void pincerror(const char* restrict format,...){

	// Retrieve argument list
	va_list args;
	va_start(args,format);

	// Print error message and quit
	fprintf(stderr,"PINC Error: ");
	vfprintf(stderr,format,args);
	fprintf(stderr,"\n");
	exit(EXIT_FAILURE);	

	va_end(args);

}
