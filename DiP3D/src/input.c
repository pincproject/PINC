/**
 * @file		input.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC main routine.
 * @date		11.10.15
 *
 * Functions for parsing input to PINC.
 * Replaces old DiP3D input.c file by Wojciech Jacek Miloch.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "iniparser.h"
#include "pinc.h"

/*
 * DECLARING LOCAL FUNCTIONS
 */

/**
 * @brief	Splits a comma-separated list-string to an array of strings.
 * @param	list	Comma-sepaprated list
 * @return	A NULL-terminated array of NULL-terminated strings
 *
 * Example: when str="abc ,def, ghi", list2strarr(str); will return
 * an array arr such that :
 *
 * arr[0]="abc"
 * arr[1]="def"
 * arr[2]="ghi"
 * arr[3]=NULL
 *
 * Note that whitespaces are trimmed away. Remember to free string array. This
 * can be done with free_strarr().
 *
 */
static char** list2strarr(const char* list);

/**
 * @brief Frees dynamically allocated array of string
 * @param	strarr	Pointer to array of strings
 * @return	void
 */
static void free_strarr(char** strarr);

/*
 * DEFINING FUNCTIONS
 */

void parse_input(int argc, char *argv[]){

	int n_specs = 0;	// Used to hold the number of specifications of a part

	// OPEN INPUT FILE

	if(argc!=2)
		pincerror("PINC expects exactly one argument (the input file).");

	dictionary* ini;
	ini = iniparser_load(argv[1]);

	if(ini==NULL) pincerror("Failed to open or parse %s.",argv[1]);

	// PARSE TIME SECTION

	int Nt = iniparser_getint(ini,"time:Nt",0);
	double T = iniparser_getdouble(ini,"time:T",0);
	double dt = iniparser_getdouble(ini,"time:dt",0);

	n_specs = (Nt!=0) + (T!=0) + (dt!=0); // Number of specified time inputs

	printf("Nt=%i, T=%f, dt=%f\n",Nt,T,dt);

	if(n_specs<2) pincerror("[time] in %s is under-determined. "
							"Specify 2 of these: Nt, T and dt.",argv[1]);

	if(n_specs>2) pincerror("[time] in %s is over-determined. "
							"Specify only 2 of these: Nt, T and dt.",argv[1]);

	if(dt==0){
		dt = T/Nt;
		printf("Computed dt to be %f.\n",dt);
	}
	if(T==0){
		T = Nt*dt;
		printf("Computed T to be %f.\n",T);
	}
	if(Nt==0){
		Nt = (int) ceil(T/dt);
		double dt_new = T/Nt;
		printf("Computed Nt to be %i.",Nt);
		if(dt!=dt_new) printf(	" Had to reduce dt from %f to %f to get "
								"integer Nt.",dt,dt_new);
		printf("\n");
		dt=dt_new;

	}

	// PARSE GRID SECTION

//	char *Ng_str = 

	// PARSE PARTICLES SECTION

	// FREE MEMORY

	iniparser_freedict(ini);

}

static char** list2strarr(const char* list){

	// Count elements in list (including NULL-element)
	int count = 2;
	char *temp = (char*)list;
	while(*temp){
		if(*temp==',') count++;
		temp++;
	}

	// Allocate NULL-terminated array of NULL-terminated strings
	char **result = 0;
	result = malloc(count*sizeof(char*));

	if(result){

		count=0;
		char finished=0;
		temp=(char*)list;
		char *start=(char*)list;
		char *stop=(char*)list;

		while(!finished){
			
			if(*temp==',' || *temp=='\0'){
				// New delimeter reached. Set stop of this element.
				stop = temp-1;

				// Trim leading and trailing whitespaces
				while(*start==' ' && start<stop) start++;
				while(*stop==' '  && start<stop) stop--;

				// Length of string (excluding NULL-termination)
				int len = stop-start+1;

				// Copy string
				result[count] = malloc((len+1)*sizeof(char));
				strncpy(result[count],start,len);
				result[count][len]='\0';

				// Set start of next element
				start = temp+1;
				count++;
			}

			// Iterate
			if(*temp=='\0') finished=1;
			temp++;
		}

	}

	result[count]=NULL;

	return result;

}

static void free_strarr(char** strarr){

	for(int i=0;strarr[i];i++) free(strarr[i]);
	free(strarr);

}
