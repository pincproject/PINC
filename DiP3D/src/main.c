/**
 * PINC (DiP3D)
 * Main file
 * Author: Sigvald Marholm
 * University of Oslo, Norway
 * 2015
 */

#include <stdlib.h>
#include <stdio.h>
#include "iniparser.h"

int main(int argc, char *argv[]){

	/**
	 * READ INPUT FILE
	 */

	if(argc!=2){
		fprintf(stderr,"Expects one and only one argument (the input file).\n");
		exit(EXIT_FAILURE);
	}

	dictionary* ini;
	ini = iniparser_load(argv[1]);

	if(ini==NULL){
		fprintf(stderr,"Could not parse %s.\n",argv[1]);
		exit(EXIT_FAILURE);
	}

	for(int i=0;i<argc;i++){
		printf("argv[%i]=%s\n",i,argv[i]);
	}

	return 0;

}


