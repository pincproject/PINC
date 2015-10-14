/**
 * @file		input.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC main routine.
 * @date		13.10.15
 *
 * Functions for parsing input to PINC.
 * Replaces old DiP3D input.c file by Wojciech Jacek Miloch.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "iniparser.h"
#include "pinc.h"

/******************************************************************************
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/

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
static char** list_to_strarr(const char* list);

/**
 * @brief Frees dynamically allocated NULL-terminated array of strings
 * @param	strarr	Pointer to array of strings
 * @return	void
 */
static void free_strarr(char** strarr);

/**
 * @brief Gets number of elements in comma-separated list-string.
 * @param	list	Comma-separated list
 * @return	Number of elements in list. 0 if string is empty.
 */
static int list_getnelements(const char* list);

/**
 * @brief Gets number of elements in comma-separated entry in ini file
 * @param	d		Dictionary to search
 * @param	key		Key string to look for in iniparser dictionary
 * @return	Number of elements in entry. 0 if entry does not exist.
 *
 * This function can be seen as an extension to iniparser in order to handle
 * comma-separated entries. It is syntactically similar to the functions in
 * iniparser.
 */
static int listparser_getnelements(const dictionary* d, const char* key);

/**
 * @brief Get the array of doubles associated to a key in ini file
 * @param		d		Dictionary to search
 * @param		key		Key string to look for in iniparser dictionary
 * @param[out]	result	Pointer to pre-allocated array to store results in
 * @return		void
 *
 * This function can be seen as an extension to iniparser in order to handle
 * comma-separated entries. It is syntactically similar to the functions in
 * iniparser.
 *
 * Array of results must be allocated before being passed to this function. The
 * number of elements to allocate can be obtained by listparser_getnelements().
 *
 * listparser_getallocdouble() is similar to this but allocates memory.
 */
static void listparser_getdouble(	const dictionary *d, const char *key, 
										double *result);

/**
 * @brief Get the array of doubles associated to a key in ini file
 * @param			d		Dictionary to search
 * @param			key		Key string to look for in iniparser dictionary
 * @param[out]		count	Number of elements in returned array
 * @return			Array of doubles
 *
 * This function can be seen as an extension to iniparser in order to handle
 * comma-separated entries. It is syntactically similar to the functions in
 * iniparser.
 *
 * listparser_getdouble() is similar to this but does not allocate memory.
 */
static double* listparser_getallocdouble(	const dictionary *d,
											const char *key, int *count);

/**
 * @brief Get the array of integers associated to a key in ini file
 * @param	d		Dictionary to search
 * @param	key		Key string to look for in iniparser dictionary
 * @param	result	Pointer to pre-allocated array to store results in
 * @return	void
 *
 * This function can be seen as an extension to iniparser in order to handle
 * comma-separated entries. It is syntactically similar to the functions in
 * iniparser.
 *
 * Array of results must be allocated before being passed to this function. The
 * number of elements to allocate can be obtained by listparser_getnelements().
 */
static void listparser_getint(	const dictionary *d, const char *key, 
										int *result);

/**
 * @brief Get the array of integers associated to a key in ini file
 * @param			d		Dictionary to search
 * @param			key		Key string to look for in iniparser dictionary
 * @param[out]		count	Number of elements in returned array
 * @return			Array of integers
 *
 * This function can be seen as an extension to iniparser in order to handle
 * comma-separated entries. It is syntactically similar to the functions in
 * iniparser.
 *
 * listparser_getint() is similar to this but does not allocate memory.
 */

static int* listparser_getallocint(	const dictionary *d, const char *key, 
									int *count);

static void fprintarr(	FILE *stream, const char* restrict format,
						void *arr, int count);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

void parse_input(int argc, char *argv[]){

	int nparams = 0;	// Used to hold the number of specified parameters

	/*
	 * OPEN INPUT FILE
	 */

	if(argc!=2)
		msg(ERROR,"exactly one argument expected (the input file).");

	dictionary* ini;
	ini = iniparser_load(argv[1]);

	if(ini==NULL) msg(ERROR,"Failed to open (or parse) %s.",argv[1]);

	msg(STATUS,"parsing %s.",argv[1]);

	/*
	 * PREPARE OUTPUT DUMP FOR INPUT PARSE
	 */

	// This file can be used to retrieve information on how parsing went.

	char *dumpfilename = iniparser_getstring(ini,"files:parsedump","");
	FILE *dumpfile;
	if(strcmp(dumpfilename,"")==0){
		msg(WARNING,"no parsedump file specified. Using stdout.");
		dumpfile = stdout;
	} else if(strcmp(dumpfilename,"stdout")==0) {
		dumpfile = stdout;
	} else {
		dumpfile = fopen(dumpfilename,"w");
	}

	/*
	 * PARSE [time] SECTION
	 */

	fprintf(dumpfile,"Parsing [time]\n");

	// Load input parameters from [time]
	int Nt = iniparser_getint(ini,"time:Nt",0);
	double T = iniparser_getdouble(ini,"time:T",0);
	double dt = iniparser_getdouble(ini,"time:dt",0);

	// Number of parameters specified
	nparams = (Nt!=0) + (T!=0) + (dt!=0);

	// Check for correct number of input parameters	
	if(nparams<2) msg(ERROR,"[time] in %s is under-determined. "
							"Specify 2 of these: Nt, T and dt.",argv[1]);

	if(nparams>2) msg(ERROR,"[time] in %s is over-determined. "
							"Specify only 2 of these: Nt, T and dt.",argv[1]);

	// Compute non-specified input parameter
	if(dt==0){
		dt = T/Nt;
		fprintf(dumpfile,"Computed dt.\n");
	}
	if(T==0){
		T = Nt*dt;
		fprintf(dumpfile,"Computed T.\n");
	}
	if(Nt==0){
		Nt = (int) ceil(T/dt);
		double dt_new = T/Nt;
		fprintf(dumpfile,"Computed Nt.\n");
		if(dt!=dt_new){
			msg(WARNING,	"had to reduce dt from %f to %f to get "
							"integer Nt.",dt,dt_new);
			fprintf(dumpfile,"Reduced dt to get integer Nt.\n");
		}
		dt=dt_new;

	}
	fprintf(dumpfile,"Nt=%i\n",Nt);
	fprintf(dumpfile,"T=%f Debye lengths (exact: %a)\n",T,T);
	fprintf(dumpfile,"dt=%f Debye lengths (exact: %a)\n",dt,dt);


	/*
	 * PARSE GRID SECTION
	 */

	fprintf(dumpfile,"Parsing [grid]\n");

	// Get dimensions specified by each parameters (0 if unspecified)
	int Ng_dim = listparser_getnelements(ini,"grid:Ng");
	int L_dim  = listparser_getnelements(ini,"grid:L");
	int dx_dim = listparser_getnelements(ini,"grid:dx");

	// Number of grid input parameters specified
	nparams = (Ng_dim!=0) + (L_dim!=0) + (dx_dim!=0);

	// Check for correct number of input parameters
	if(nparams<2) msg(ERROR,"[grid] in %s is under-determined. "
							"Specify 2 of these: Ng, L and dx.",argv[1]);

	if(nparams>2) msg(ERROR,"[grid] in %s is over-determined. "
							"Specify only 2 of these: Ng, L and dx.",argv[1]);

	// Check for equal length of lists
	if(Ng_dim==0 && L_dim!=dx_dim)
		msg(ERROR,"L and dx have unequal number of elements.");

	if(L_dim==0 && Ng_dim!=dx_dim)
		msg(ERROR,"Ng and dx have unequal number of elements.");

	if(dx_dim==0 && Ng_dim!=L_dim)
		msg(ERROR,"Ng and L have unequal number of elements.");

	// Get number of dimensions (one is zero, the other two equals dim)
	int dim = (Ng_dim + L_dim + dx_dim)/2;

	fprintf(dumpfile,"Dimensions=%i\n",dim);

	// Check for valid number of dimensions
	if(dim!=3) msg(ERROR,	"%i dimensions specified in [grid] but only 3D "
							"simulations are supported.",dim);

	// Load specified input parameters
	int 	*Ng = malloc(dim*sizeof(int));
	double	*L  = malloc(dim*sizeof(double));
	double	*dx = malloc(dim*sizeof(double));
	listparser_getint(ini,"grid:Ng",Ng);
	listparser_getdouble(ini,"grid:L",L);
	listparser_getdouble(ini,"grid:dx",dx);

	// Compute non-specified input parameter
	if(dx_dim==0){
		fprintf(dumpfile,"Computed dx.\n");
		for(int i=0;i<dim;i++){
			dx[i] = L[i]/Ng[i];
			//printf("	dx[%i]=%f.\n",i,dx[i]);
		}
	}
	if(L_dim==0){
		fprintf(dumpfile,"Computed L.\n");
		for(int i=0;i<dim;i++){
			L[i] = Ng[i]*dx[i];
			//printf("	L[%i]=%f.\n",i,L[i]);
		}
	}
	if(Ng_dim==0){
		fprintf(dumpfile,"Computed Ng.\n");
		for(int i=0;i<dim;i++){
			Ng[i] = (int) ceil(L[i]/dx[i]);
			double dx_new = L[i]/Ng[i];
			//printf("	Ng[%i]=%i.",i,Ng[i]);
			if(dx[i]!=dx_new){
				msg(WARNING,	"had to reduce dx[%i] from %f to %f to get "
								"integer Ng[%i].",i,dx[i],dx_new,i);
				fprintf(dumpfile,
					"Reduced dx[%i] to get integer Ng[%i].\n",i,i);
			}
			dx[i]=dx_new;
		}
	}

	fprintf(dumpfile,"Ng=");
	fprintarr(dumpfile,"%i",(void*)Ng,dim);
	fprintf(dumpfile,"\n");
	/*
	printf("Ng=%i,%i,%i\n",Ng[0],Ng[1],Ng[2]);
	printf("L=%f,%f,%f\n",L[0],L[1],L[2]);
	printf("dx=%f,%f,%f\n",dx[0],dx[1],dx[2]);
	*/

	// PARSE PARTICLES SECTION

	// Get number of particles per specie (Nps) and number of species (Ns)
	int Ns=0, Ns_check=0;
	int *Nps = listparser_getallocint(ini,"particles:Nps",&Ns);

	// Get the specie charges
	double *Q = listparser_getallocdouble(ini,"particles:Q",&Ns_check);
	if(Ns_check!=Ns) msg(ERROR,"Q and Nps have unequal number of elements.");

	// Get the specie masses
	double *M = listparser_getallocdouble(ini,"particles:M",&Ns_check);
	if(Ns_check!=Ns) msg(ERROR,"M and Nps have unequal number of elements.");

	// CLOSE DUMPFILE

	if(strcmp(dumpfilename,"") && strcmp(dumpfilename,"stdout"))
		fclose(dumpfile);

	// FREE MEMORY

	iniparser_freedict(ini);

	free(Nps);

	free(Ng);
	free(L);
	free(dx);

}

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

static void fprintarr(	FILE *stream, const char* restrict format,
						void *arr, int count){

	barr

	fprintf(stream,"%i",barr[0]);
	for(int i=1;i<count;i++) fprintf(stream,",%i",barr[i]);

}

static int list_getnelements(const char* list){

	if(list[0]=='\0') return 0;	// key not found

	// Count elements
	int count = 1;
	for(int i=0;list[i];i++) count += (list[i]==',');

	return count;

}

static int listparser_getnelements(const dictionary* d, const char* key){

	char *list = iniparser_getstring((dictionary*)d,key,"");
	return list_getnelements(list);

}

static void listparser_getint(	const dictionary *d, const char *key,
								int *result){

	char *list = iniparser_getstring((dictionary*)d,key,"");
	char **strarr = list_to_strarr(list);
	int count = list_getnelements(list);

	for(int i=0;i<count;i++) result[i] = (int)strtol(strarr[i],NULL,10);

	free_strarr(strarr);

}

static void listparser_getdouble(	const dictionary *d, const char *key, 
									double *result){

	char *list = iniparser_getstring((dictionary*)d,key,"");
	char **strarr = list_to_strarr(list);
	int count = list_getnelements(list);

	for(int i=0;i<count;i++) result[i] = strtod(strarr[i],NULL);

	free_strarr(strarr);

}

static int* listparser_getallocint(	const dictionary *d, const char *key, 
									int *count){

	*count = listparser_getnelements(d,key);
	int *result = malloc(*count*sizeof(int));

	listparser_getint(d,key,result);

	return result;

}

static double* listparser_getallocdouble(	const dictionary *d,
											const char *key, int *count){

	*count = listparser_getnelements(d,key);
	double *result = malloc(*count*sizeof(double));

	listparser_getdouble(d,key,result);

	return result;

}


static char** list_to_strarr(const char* list){

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

