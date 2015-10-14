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
#include <string.h>

/******************************************************************************
 * GLOBAL DEFINITIONS
 *****************************************************************************/

void msg(msg_kind kind, const char* restrict format,...){

	// Retrieve argument list
	va_list args;
	va_start(args,format);

	// Set prefix to output and output stream
	char *prefix = malloc(20*sizeof(char));
	FILE *stream;
	switch(kind){
		case STATUS:
			strcpy(prefix,"STATUS: ");
			stream=stdout;
			break;
		case WARNING:
			strcpy(prefix,"WARNING: ");
			stream=stderr;
			break;
		case ERROR:
			strcpy(prefix,"ERROR: ");
			stream=stderr;
			break;
	}

	// Print message
	fprintf(stream,"%s",prefix);
	vfprintf(stream,format,args);
	fprintf(stream,"\n");

	// Stop vararg
	va_end(args);

	// Quit if error
	if(kind==ERROR) exit(EXIT_FAILURE);

}
