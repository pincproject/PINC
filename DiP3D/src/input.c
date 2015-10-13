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
//static char** list2strarr(const char* list);

/**
 * @brief Frees dynamically allocated array of string
 * @param	strarr	Pointer to array of strings
 * @return	void
 */
//static void free_strarr(char** strarr);

//static int* list2intarr(const char* list, int *count);
static int listparser_getnelements(dictionary* d, const char* key);

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

	if(n_specs<2) pincerror("[time] in %s is under-determined. "
							"Specify 2 of these: Nt, T and dt.",argv[1]);

	if(n_specs>2) pincerror("[time] in %s is over-determined. "
							"Specify only 2 of these: Nt, T and dt.",argv[1]);

	if(dt==0){
		dt = T/Nt;
		printf("T and Nt specified. Computed dt to be %f.\n",dt);
	}
	if(T==0){
		T = Nt*dt;
		printf("Nt and dt specified. Computed T to be %f.\n",T);
	}
	if(Nt==0){
		Nt = (int) ceil(T/dt);
		double dt_new = T/Nt;
		printf("T and dt specified. Computed Nt to be %i.",Nt);
		if(dt!=dt_new) printf(	" Had to reduce dt from %f to %f to get "
								"integer Nt.",dt,dt_new);
		printf("\n");
		dt=dt_new;

	}

	printf("Number of time steps: 	Nt = %i\n",Nt);
	printf("Total simulation time: 	T  = %f Debye lengths\n",T);
	printf("Time step: 		dt = %f Debye lengths\n",dt);

	// PARSE GRID SECTION

	int dim=0;
	int Ng_dim = listparser_getnelements(ini,"grid:Ng");
	int L_dim  = listparser_getnelements(ini,"grid:L");
	int dx_dim = listparser_getnelements(ini,"grid:dx");

	// Number of grid input specifications
	n_specs = (Ng_dim!=0) + (L_dim!=0) + (dx_dim!=0);

	if(n_specs<2) pincerror("[grid] in %s is under-determined. "
							"Specify 2 of these: Ng, L and dx.",argv[1]);

	if(n_specs>2) pincerror("[grid] in %s is over-determined. "
							"Specify only 2 of these: Ng, L and dx.",argv[1]);

	if(Ng_dim==0 && L_dim!=dx_dim)
		pincerror("L and dx have unequal number of elements.");

	if(L_dim==0 && Ng_dim!=dx_dim)
		pincerror("Ng and dx have unequal number of elements.");

	if(dx_dim==0 && Ng_dim!=L_dim)
		pincerror("Ng and L have unequal number of elements.");

	dim = (Ng_dim+L_dim+dx_dim)/2;
	printf("Grid is %i-dimensional.\n",dim);

	int *Ng = malloc(dim*sizeof(int));
	int *L  = malloc(dim*sizeof(int));
	int *dx = malloc(dim*sizeof(int));
		

	char *Ng_str = iniparser_getstring(ini,"grid:Ng","");
	char *L_str  = iniparser_getstring(ini,"grid:L","");
	char *dx_str = iniparser_getstring(ini,"grid:dx","");

	// Number of specified grid inputs
	n_specs = (*Ng_str!="") + (*L_str!="") + (*dx_str!="");

/*
	int Ng_dim = 0, L_dim = 0, dx_dim = 0, dim=0;
	int *Ng = list2intarr(Ng_str,&Ng_dim);
	int *L  = list2intarr(L_str ,&L_dim );
	int *dx = list2intarr(dx_str,&dx_dim);

	printf("Number of dimensions: %i, %i, %i\n",Ng_dim,L_dim,dx_dim);

	if(dx==0){
		if(Ng_dim!=dx_dim)
			pincerror("Ng and dx have unequal number of elements.");
		dim = Ng_dim;
		for(int i=0;i<dim;i++)
		printf("T and Nt specified. Computed dt to be %f.\n",dt);
	}
	if(T==0){
		T = Nt*dt;
		printf("Nt and dt specified. Computed T to be %f.\n",T);
	}
	if(Nt==0){
		Nt = (int) ceil(T/dt);
		double dt_new = T/Nt;
		printf("T and dt specified. Computed Nt to be %i.",Nt);
		if(dt!=dt_new) printf(	" Had to reduce dt from %f to %f to get "
								"integer Nt.",dt,dt_new);
		printf("\n");
		dt=dt_new;

	}
*/
	// PARSE PARTICLES SECTION

	// FREE MEMORY
/*
	free(Ng);
	free(L);
	free(dx);
*/
	iniparser_freedict(ini);

}

static int listparser_getnelements(dictionary* d, const char* key){

	char *list = iniparser_getstring(d,key,"");
	
	if(list[0]=='\0') return 0;	// key not found

	// Count elements
	int count = 1;
	for(int i=0;list[i];i++) count += (list[i]==',');

	return count;

}

/*
static int listcount(const char* list){

	if(list[0]=='\0') return 0;

	int count = 1;
	for(int i=0;list[i];i++) count += (list[i]==',');

}

static int* list2intarr(const char* list, int *count){

	// Count elements in list
	*count = 1;
	for(int i=0;list[i];i++) *count += (list[i]==',');

	int *intarr = 0;
	intarr = malloc(*count*sizeof(int));

	char *temp = (char*)list;

	for(int i=0;i<*count;i++)
		intarr[i] = (int)strtol(temp,&temp,10);

	return intarr;
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
*/
