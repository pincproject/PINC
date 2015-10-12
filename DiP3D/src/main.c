/**
 * @file	    main.c
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright   University of Oslo, Norway
 * @brief	    PINC main routine.
 * @date        08.10.15
 *
 * Main routine for PINC (Particle-IN-Cell).
 * Replaces old DiP3D main.c file by Wojciech Jacek Miloch.
 */

#include <stdlib.h>
#include <stdio.h>
#include "pinc.h"

int main(int argc, char *argv[]){

	/*
	 * READ INPUT FILE
	 */
    parse_input(argc,argv);


	/*
	 * SUCCESSFUL EXIT
	 */

	printf("PINC completed successfully!\n");
	return 0;

}


