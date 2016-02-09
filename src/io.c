/**
 * @file		io.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		PINC input/output handling.
 * @date		13.10.15
 *
 * Functions for dealing with input and output. Dealing with input/output is in
 * PINC largely done by the iniparser library and the HDF5 API. io.c may expand
 * upon this with extra functionality suitable for PINC. Moreover, io.c
 * facilitates standardized means of output printing.
 *
 * This file is _not_ intended to do the reading and writing of other modules'
 * data. They should be somewhat self-contained, merely using this file kind of
 * like a library of supporting functions. This file should not depend on any
 * of those modules in case they are to be replaced in the future.
 */

#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <mpi.h>
#include <hdf5.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "iniparser.h"
#include "pinc.h"

#define BUFFSIZE 128

/******************************************************************************
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/

 /**
  * @brief Makes a directory
  * @param	dir		Directory name
  * @return	0 for success, 1 for failure
  *
  * dir can be a path but ancestors must exist. Directory will have permissions
  * 0775. Used in makePath().
  *
  * @see makePath().
  */
 static int makeDir(const char *dir);

 /**
  * @brief Makes all parent directories of URL path
  * @param	path	Path
  * @return	0 for success, 1 for failure
  *
  * Examples:
  *	path="dir/dir/file"	generates the folder "./dir/dir/"
  *  path="dir/dir/dir/" generates the folder "./dir/dir/dir/"
  *  path="../dir/file" generates the folder "../dir/"
  *	path="/dir/file" generates the folder "/dir/"
  *
  * Already existing folders are left as-is. This function can be used to ensure
  * that the parent directories of its path exists.
  */
 int makePath(const char *path);

/**
 * @brief	Splits a comma-separated list-string to an array of strings.
 * @param	list	Comma-sepaprated list
 * @return	A NULL-terminated array of NULL-terminated strings
 * @see		freeStrArr(), listGetNElements()
 *
 * Example: when str="abc ,def, ghi", listToStrArr(str); will return
 * an array arr such that :
 *
 * arr[0]="abc"
 * arr[1]="def"
 * arr[2]="ghi"
 * arr[3]=NULL
 *
 * Note that whitespaces are trimmed away. Remember to free string array with
 * freeStrArr().
 */
static char** listToStrArr(const char* list);

/**
 * @brief Gets number of elements in comma-separated list-string.
 * @param	list	Comma-separated list
 * @return	Number of elements in list. 0 if string is empty.
 * @see listGetNElements(), iniGetNElements()
 */
static int listGetNElements(const char* list);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

dictionary* iniOpen(int argc, char *argv[]){

	// Sanity check on input arguments
	if(argc<2)
		msg(ERROR,"at least one argument expected (the input file).");

	// Open ini-file
	dictionary *ini = iniparser_load(argv[1]);
	if(ini==NULL) msg(ERROR,"Failed to open %s.",argv[1]);

	for(int i=2;i<argc;i++){
		if(!strcmp(argv[i],"getnp")){	// Just returning number of processes
			int nDims;
			int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",&nDims);
			int np = aiProd(nSubdomains,nDims);
			printf("%i\n",np);
			free(nSubdomains);
			exit(0);
		} else {
			char *value = strstr(argv[i],"=");
			*value = '\0';
			value++;
			iniparser_set(ini,argv[i],value);
		}
	}

	// Start new fmsg()-files (iterate through all files in [msgfiles] section)
	int nKeys = iniparser_getsecnkeys(ini,"msgfiles");
	char **keys = iniparser_getseckeys(ini,"msgfiles");
	for(int i=0;i<nKeys;i++){

		// Get filename corresponding to key
		char *fName = iniparser_getstring(ini,keys[i],""); // don't free

		// Make file empty (unless using stdout or stderr)
		if(strcmp(fName,"")==0){
			msg(WARNING|ONCE,"%s not specified. Using stdout.",keys[i]);
		} else if(strcmp(fName,"stdout")!=0 && strcmp(fName,"stderr")!=0) {

			if(makePath(fName))
				msg(ERROR|ONCE,"Could not open or create path of '%s'",fName);

			// check whether file exists
			FILE *fh = fopen(fName,"r");
			if(fh!=NULL){
				fclose(fh);
				msg(ERROR|ONCE,"'%s' already exists.",fName);
			}
		}

	}

	// Note that keys should be freed but not keys[0] and so on.
	free(keys);

	return ini;

}

void msg(msgKind kind, const char* restrict format,...){

	// Retrieve argument list
	va_list args;
	va_start(args,format);
	const int bufferSize = 132;

	// Set prefix and determine which output to use
	char prefix[8];
	FILE *stream;
	switch(kind&0x0F){
		case STATUS:
			strcpy(prefix,"STATUS");
			stream=stdout;
			break;
		case WARNING:
			strcpy(prefix,"WARNING");
			stream=stderr;
			break;
		case ERROR:
			strcpy(prefix,"ERROR");
			stream=stderr;
			break;
	}

	// Parse and assemble message
	int rank;
	char msg[bufferSize], buffer[bufferSize];
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vsnprintf(msg,bufferSize,format,args);
	snprintf(buffer,bufferSize,"%s (%i): %s",prefix,rank,msg);
	va_end(args);

	// Print message
	if(!(kind&ONCE) || rank==0){
		fprintf(stream,"%s\n",buffer);
	}

	// Quit if error
	if((kind&0x0F)==ERROR) exit(EXIT_FAILURE);

}

void fMsg(dictionary *ini, const char* restrict fNameKey, const char* restrict format, ...){

	// Get filename
	char key[BUFFSIZE] = "msgfiles:";
	strcat(key,fNameKey);
	char *fName = iniparser_getstring(ini,key,"");

	// Have you opened a file that must be closed?
	int fileOpen = 0;

	// Open file (or other stream)
	FILE *file;
	if(strcmp(fName,"stdout")==0)		file = stdout;
	else if(strcmp(fName,"stderr")==0)	file = stderr;
	else if(strcmp(fName,"")==0)		file = stdout;
	else {
		file = fopen(fName,"a");
		fileOpen = 1;
	}

	// Print
	va_list args;
	va_start(args,format);
	vfprintf(file,format,args);
	va_end(args);

	// Close file
	if(fileOpen) fclose(file);

}

/******************************************************************************
 * DEFINING INI PARSING FUNCTIONS (expanding iniparser library)
 *****************************************************************************/

int iniGetNElements(const dictionary* ini, const char* key){

	char *list = iniparser_getstring((dictionary*)ini,key,"");
	return listGetNElements(list);

}


char** iniGetStrArr(const dictionary *ini, const char *key, int *nElements){

	char *list = iniparser_getstring((dictionary*)ini,key,"");	// don't free this
	char **strArr = listToStrArr(list);
	*nElements = listGetNElements(list);

	return strArr;

}

long int* iniGetLongIntArr(const dictionary *ini, const char *key, int *nElements){

	char **strArr = iniGetStrArr(ini,key,nElements);

	long int *result = malloc(*nElements*sizeof(long int));
	for(int i=0;i<*nElements;i++) result[i] = strtol(strArr[i],NULL,10);

	freeStrArr(strArr);

	return result;

}

int* iniGetIntArr(const dictionary *ini, const char *key, int *nElements){

	char **strArr = iniGetStrArr(ini,key,nElements);

	int *result = malloc(*nElements*sizeof(int));
	for(int i=0;i<*nElements;i++) result[i] = (int)strtol(strArr[i],NULL,10);

	freeStrArr(strArr);

	return result;

}

double* iniGetDoubleArr(const dictionary *ini, const char *key, int *nElements){

	char **strArr = iniGetStrArr(ini,key,nElements);

	double *result = malloc(*nElements*sizeof(double));
	for(int i=0;i<*nElements;i++) result[i] = strtod(strArr[i],NULL);

	freeStrArr(strArr);

	return result;

}

int iniAssertEqualNElements(const dictionary *ini, int nKeys, ...){

	va_list args;
	va_start(args,nKeys);

	char buffer[BUFFSIZE] = "";
	char *key = va_arg(args,char*);
	sprintf(buffer,"%s%s ",buffer,key);
	int nElements = iniGetNElements(ini,key);
	int equal = 1;
	for(int i=1;i<nKeys;i++){
		key = va_arg(args,char*);
		sprintf(buffer,"%s%s ",buffer,key);
		int temp = iniGetNElements(ini,key);
		if(temp!=nElements) equal=0;
	}

	va_end(args);

	if(!equal){
		sprintf(buffer,"%smust have equal length.",buffer);
		msg(ERROR,buffer);
	}

	return nElements;

}

/******************************************************************************
 * DEFINING HDF5 FUNCTIONS (expanding HDF5 API)
 *****************************************************************************/

hid_t createH5File(const dictionary *ini, const char *fName, const char *fSubExt){


	// Determine filename
	char *fPrefix = iniparser_getstring((dictionary *)ini,"files:output","");	// don't free

	// Add separator if filename prefix (not just folder) is specified
	char sep[2] = "\0\0";
	char lastchar = fPrefix[strlen(fPrefix)-1];
	if(strcmp(fPrefix,".")==0) sep[0]='/';
	else if(lastchar!='/') sep[0]='_';

	char *fTotName = strCatAlloc(6,fPrefix,sep,fName,".",fSubExt,".h5");

	// Create file with MPI-I/O access
	hid_t pList = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pList,MPI_COMM_WORLD,MPI_INFO_NULL);


	// Make sure parent folder exist
	if(makePath(fName))
		msg(ERROR|ONCE,"Could not open or create folder for '%s'.",fTotName);

	// check whether file exists
	FILE *fh = fopen(fTotName,"r");
	if(fh!=NULL){
		fclose(fh);
		msg(ERROR|ONCE,"'%s' already exists.",fTotName);
	}

	// create file
	hid_t file = H5Fcreate(fTotName,H5F_ACC_EXCL,H5P_DEFAULT,pList);
	H5Pclose(pList);

	free(fTotName);

	return file;

}

/******************************************************************************
 * DEFINING LOCAL LIST PARSING FUNCTIONS
 *****************************************************************************/

static int listGetNElements(const char* list){

	if(list[0]=='\0') return 0;	// key not found

	// Count elements
	int nElements = 1;
	for(int i=0;list[i];i++) nElements += (list[i]==',');

	return nElements;

}

static char** listToStrArr(const char* restrict list){

	// Count elements in list (including NULL-element)
	int nElements = 2;				// first element and NULL
	char *temp = (char*)list;
	while(*temp){
		if(*temp==',') nElements++;	// one more element per delimeter
		temp++;
	}

	// Allocate NULL-terminated array of NULL-terminated strings
	char **result = 0;
	result = malloc(nElements*sizeof(char*));

	if(result){

		nElements=0;
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
				result[nElements] = malloc((len+1)*sizeof(char));
				strncpy(result[nElements],start,len);
				result[nElements][len]='\0';

				// Set start of next element
				start = temp+1;
				nElements++;
			}

			// Iterate
			if(*temp=='\0') finished=1;
			temp++;
		}

	}

	result[nElements]=NULL;

	return result;

}

void freeStrArr(char** strArr){

	for(int i=0;strArr[i];i++) free(strArr[i]);
	free(strArr);

}

/******************************************************************************
 * DEPRECATED FUNCTIONS
 *****************************************************************************/

/*
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
*/

/*
void ini_complete_time(dictionary *ini){

	// Load input parameters from [time]
	int Nt = iniparser_getint(ini,"time:Nt",0);
	double T = iniparser_getdouble(ini,"time:T",0);
	double dt = iniparser_getdouble(ini,"time:dt",0);

	// Number of parameters specified
	int nparams = (Nt!=0) + (T!=0) + (dt!=0);

	// Check for correct number of input parameters
	if(nparams<2) msg(ERROR,"[time] is under-determined. Specify 2 of these: Nt, T and dt.");
	if(nparams>2) msg(ERROR,"[time] is over-determined. Specify only 2 of these: Nt, T and dt.");

	// Compute non-specified input parameter and store in dictionary
	// Hexadecimal (%a) specifier is used since this cause no precision loss when converting to string
	if(dt==0){
		dt = T/Nt;

		char buffer[BUFFSIZE];
		sprintf(buffer,"%a",dt);
		iniparser_set(ini,"time:dt",buffer);

		fmsg(ini,"parsedump","Computed dt.\n");
	}
	if(T==0){
		T = Nt*dt;

		char buffer[BUFFSIZE];
		sprintf(buffer,"%a",T);
		iniparser_set(ini,"time:T",buffer);

		fmsg(ini,"parsedump","Computed T.\n");
	}
	if(Nt==0){
		Nt = (int) ceil(T/dt);
		double dt_new = T/Nt;

		char buffer[BUFFSIZE];
		sprintf(buffer,"%i",Nt);
		iniparser_set(ini,"time:Nt",buffer);
		sprintf(buffer,"%a",dt_new);
		iniparser_set(ini,"time:dt",buffer);

		fmsg(ini,"parsedump","Computed Nt.\n");
		if(dt!=dt_new){
			msg(WARNING,			"had to reduce dt from %f to %f to get integer Nt.",dt,dt_new);
			fmsg(ini,"parsedump",	"Had to reduce dt from %f to %f to get integer Nt.\n",dt,dt_new);
		}

	}
	fmsg(ini,"parsedump","Nt=%i\n",Nt);
	fmsg(ini,"parsedump","T=%f 1/omega_p (exact: %a)\n",T,T);
	fmsg(ini,"parsedump","dt=%f 1/omega_p (exact: %a)\n",dt,dt);

}

void ini_complete_grid(dictionary *ini){

	// Get dimensions specified by each parameters (0 if unspecified)
	int Ng_dim = listparser_getnelements(ini,"grid:Ng");
	int L_dim  = listparser_getnelements(ini,"grid:L");
	int dx_dim = listparser_getnelements(ini,"grid:dx");

	// Number of grid input parameters specified
	int nparams = (Ng_dim!=0) + (L_dim!=0) + (dx_dim!=0);

	// Check for correct number of input parameters
	if(nparams<2) msg(ERROR,"[grid] is under-determined. Specify 2 of these: Ng, L and dx.");
	if(nparams>2) msg(ERROR,"[grid] is over-determined. Specify only 2 of these: Ng, L and dx.");

	// Check for equal length of lists
	if(Ng_dim==0 && L_dim!=dx_dim)
		msg(ERROR,"L and dx have unequal number of elements.");

	if(L_dim==0 && Ng_dim!=dx_dim)
		msg(ERROR,"Ng and dx have unequal number of elements.");

	if(dx_dim==0 && Ng_dim!=L_dim)
		msg(ERROR,"Ng and L have unequal number of elements.");

	// Get number of dimensions (one is zero, the other two equals dim)
	int dim = (Ng_dim + L_dim + dx_dim)/2;

	fmsg(ini,"parsedump","Nd=%i dimensions\n",dim);

	// Check for valid number of dimensions
	if(dim!=3) msg(ERROR,"%i dimensions specified but only 3D simulations are supported.",dim);

	// Load specified input parameters
	int 	*Ng = malloc(dim*sizeof(int));
	double	*L  = malloc(dim*sizeof(double));
	double	*dx = malloc(dim*sizeof(double));
	listparser_getint(ini,"grid:Ng",Ng);
	listparser_getdouble(ini,"grid:L",L);
	listparser_getdouble(ini,"grid:dx",dx);

	// Compute non-specified input parameter and store in dictionary
	// Hexadecimal (%a) specifier is used since this cause no precision loss when converting to string
	if(dx_dim==0){
		fmsg(ini,"parsedump","Computed dx.\n");

		char buffer[BUFFSIZE]="";
		for(int i=0;i<dim;i++){
			dx[i] = L[i]/Ng[i];
			sprintf(buffer,"%s,%a",buffer,dx[i]);
		}

		char *temp=buffer+1;	// Skip first comma
		iniparser_set(ini,"grid:dx",temp);

	}
	if(L_dim==0){
		fmsg(ini,"parsedump","Computed L.\n");

		char buffer[BUFFSIZE]="";
		for(int i=0;i<dim;i++){
			L[i] = Ng[i]*dx[i];
			sprintf(buffer,"%s,%a",buffer,L[i]);
		}

		char *temp=buffer+1;	// Skip first comma
		iniparser_set(ini,"grid:L",temp);
	}
	if(Ng_dim==0){
		fmsg(ini,"parsedump","Computed Ng.\n");

		char buffer_Ng[BUFFSIZE]="";
		char buffer_dx[BUFFSIZE]="";
		for(int i=0;i<dim;i++){
			Ng[i] = (int) ceil(L[i]/dx[i]);
			double dx_new = L[i]/Ng[i];

			sprintf(buffer_Ng,"%s,%i",buffer_Ng,Ng[i]);
			sprintf(buffer_dx,"%s,%a",buffer_dx,dx_new);

			if(dx[i]!=dx_new){
				msg(WARNING,			"had to reduce dx[%i] from %f to %f to get integer Ng[%i].",i,dx[i],dx_new,i);
				fmsg(ini,"parsedump",	"Had to reduce dx[%i] from %f to %f to get integer Ng[%i].\n",i,dx[i],dx_new,i);
			}
		}

		char *temp_Ng=buffer_Ng+1;		// Skip first comma
		char *temp_dx=buffer_dx+1;
		iniparser_set(ini,"grid:Ng",temp_Ng);
		iniparser_set(ini,"grid:dx",temp_dx);

	}

	fmsg(ini,"parsedump","Ng=%f",Ng[0]);
	for(int i=1;i<dim;i++) fmsg(ini,"parsedump",",%f",Ng[i]);
	fmsg(ini,"parsedump","\n");

	fmsg(ini,"parsedump","L=%f",L[0]);
	for(int i=1;i<dim;i++) fmsg(ini,"parsedump",",%f",L[i]);
	fmsg(ini,"parsedump"," (exact: %a",L[0]);
	for(int i=1;i<dim;i++) fmsg(ini,"parsedump",",%a",L[i]);
	fmsg(ini,"parsedump",")\n");

	fmsg(ini,"parsedump","dx=%f",dx[0]);
	for(int i=1;i<dim;i++) fmsg(ini,"parsedump",",%f",dx[i]);
	fmsg(ini,"parsedump"," (exact: %a",dx[0]);
	for(int i=1;i<dim;i++) fmsg(ini,"parsedump",",%a",dx[i]);
	fmsg(ini,"parsedump",")\n");

	free(Ng);
	free(L);
	free(dx);

}
*/

static int makeDir(const char *dir){
    struct stat st;
    int status = 0;

    if(stat(dir, &st) != 0){
        if (mkdir(dir, 0775) != 0 && errno != EEXIST) status = -1;
    } else if(!S_ISDIR(st.st_mode)) {
        errno = ENOTDIR;
        status = -1;
    }

    return(status);
}

int makePath(const char *path){
    char *pp;
    char *sp;
    int status;
    char *copy = strdup(path);

    status = 0;
    pp = copy;
    while (status == 0 && (sp = strchr(pp, '/')) != 0){
        if (sp != pp){
            *sp = '\0';
            status = makeDir(copy);
            *sp = '/';
        }
        pp = sp + 1;
    }

    free(copy);
    return(status);
}
