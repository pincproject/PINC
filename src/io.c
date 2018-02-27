/**
 * @file		io.c
 * @brief		PINC input/output handling.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
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
#include "core.h"

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
 *	path="dir/dir/dir/" generates the folder "./dir/dir/dir/"
 *	path="../dir/file" generates the folder "../dir/"
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

/**
 * @brief Asserts that key exists in ini-file and emits ERROR if not.
 * @param	ini		ini-file dictionary
 * @param	key		Key to check existence of
 * @return			void
 */
void iniAssertExistence(const dictionary *ini, const char* key);

/**
 * @brief Expands/repeats a string array
 * @param	strArr		String array to expand
 * @param	nElements	Elements to expand to
 * @return	Expanded string array
 *
 * Allocates a new string array and repeats the values of strArr until the new
 * array is filled. E.g. if strArr= {"a","b"} and nElements = 5 then the result
 * will be {"a","b","a","b","a"}. Remember to free this newly allocated array
 * using freeStrArr().
 */
char **strArrExpand(char **strArr, int nElements);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

funPtr selectInner(dictionary *ini, const char *key, const char *list,...){

	va_list args;
	va_start(args,list);

	char *value = iniGetStr(ini,key);
	funPtr (*setFunction)() = NULL;

	// "list" is all variadic arguments stringified by select() macro, e.g.
	// "fun1_set,fun2_set". Separate it into a string array.
	char **strArr = listToStrArr(list);

	/*
	 * MATCHING FUNCTION POINTER TO VALUE IN INI-FILE
	 */
	char **strTemp = strArr;
	for(char *str=*strTemp; str!=NULL; str=*(++strTemp) ){

		strtok(str,"_"); // Trim _set part
		funPtr (*fun)() = va_arg(args,funPtr (*)());

		if(!strcmp(str,value)){
			setFunction = fun;
		}
	}

	/*
	 * ERROR HANDLING (Print list of valid values in case of error)
	 */
	if(setFunction==NULL){

		char *valid = malloc(strlen(list+1)*sizeof(*valid));
		valid[0] = '\0';

		char **strTemp = strArr; // _set already trimmed
		for(char *str=*strTemp; str!=NULL; str=*(++strTemp) ){
			strcat(valid," ");
			strcat(valid,str);
		}

		msg(ERROR,"%s=%s invalid. Valid arguments:%s.",key,value,valid);
		free(valid);
	}

	/*
	 * EXECUTING _set()-FUNCTION
	 */

	freeStrArr(strArr);
	free(value);
	va_end(args);

	return setFunction(ini);
}

void msg(msgKind kind, const char* restrict format,...){

	// Retrieve argument list
	va_list args;
	va_start(args,format);
	const int bufferSize = 132;

	// Set prefix and determine which output to use
	char prefix[8];
	FILE *stream = stdout;
	switch(kind&0x0F){
		case STATUS:
			strcpy(prefix, "STATUS");
			break;
		case WARNING:
			strcpy(prefix, "WARNING");
			stream = stderr;
			break;
		case ERROR:
			strcpy(prefix, "ERROR");
			stream = stderr;
			break;
	    case TIMER:
		    strcpy(prefix, "TIMER");
		    break;
	}

	// Parse and assemble message
	int rank;
	char msg[bufferSize], buffer[bufferSize];
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	vsnprintf(msg,bufferSize,format,args);
	if((kind&ALL)){
		snprintf(buffer,bufferSize,"%s (%i): %s",prefix,rank,msg);
	} else {
		snprintf(buffer,bufferSize,"%s: %s",prefix,msg);
	}
	va_end(args);

	// Print message
	if((kind&ALL) || rank==0){
		fprintf(stream,"%s\n",buffer);
	}

	// Quit if error
	if((kind&0x0F)==ERROR) exit(EXIT_FAILURE);

}

void fMsg(dictionary *ini, const char* restrict fNameKey, const char* restrict format, ...){

	// Get filename
	char key[BUFFSIZE] = "msgfiles:";
	strcat(key,fNameKey);
	char *fName = iniGetStr(ini,key);

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
	free(fName);
}

void parseIndirectInput(dictionary *ini){

	/*
	 * APPLIES MULTIPLIERS TO ELEMENTS WITH SUFFICES
	 */

	int nDims = iniGetInt(ini, "grid:nDims");

	double V = (double)gGetGlobalVolume(ini);

	int *L = gGetGlobalSize(ini);
	double *mul = malloc(nDims*sizeof(nDims));
	for(int i=0;i<nDims;i++) mul[i] = 1.0/L[i];

	iniApplySuffix(ini, "population:nParticles", "pc", &V, 1);
	iniApplySuffix(ini, "population:nAlloc", "pc", &V, 1);
	iniApplySuffix(ini, "grid:nEmigrantsAlloc", "pc", &V, 1);
	iniApplySuffix(ini, "grid:stepSize", "tot", mul, nDims);

	free(mul);
}

/******************************************************************************
 * DEFINING NORMALIZATION FUNCTIONS
 *****************************************************************************/

Scales *normalizeSemiSI(dictionary *ini){

	const double elementaryCharge = 1.60217733e-19; // [C]
	const double electronMass = 9.10938188e-31; // [kg]
	const double vacuumPermittivity = 8.854187817e-12; // [F/m]

	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);
	double timeStep = iniGetDouble(ini, "time:timeStep");

	const double tol = 1e-10;
	if(abs(charge[0]+1)>tol)
		msg(ERROR, "Species 0 must have charge -1 with this normalization");
	if(abs(mass[0]-1)>tol)
		msg(ERROR, "Species 0 must have mass 1 with this normalization");

	adScale(charge, nSpecies, elementaryCharge);
	adScale(mass, nSpecies, electronMass);

	double wpe = sqrt(pow(elementaryCharge,2)*density[0]/
			(vacuumPermittivity*electronMass));
	timeStep /= wpe;

	iniSetDoubleArr(ini, "population:charge", charge, nSpecies);
	iniSetDoubleArr(ini, "population:mass", mass, nSpecies);
	iniSetDouble(ini, "time:timeStep", timeStep);

	free(charge);
	free(mass);
	free(density);

	return normalizeSI(ini);
}

void nFree(Scales *scales){
	free(scales);
}

static void nAddDerivedUnits(Scales *scales){

	double nDims = scales->nDims;
	double charge = scales->charge;
	double mass = scales->mass;
	double length = scales->length;
	double time = scales->time;

	scales->velocity = length/time;
	scales->acceleration = length/pow(time,2);
	scales->chargeDensity = charge/pow(length,nDims);
	scales->potential = pow(length/time,2)*mass/charge;
	scales->eField = length*mass/(pow(time,2)*charge);
	scales->bField = mass/(time*charge);
	scales->energy = mass*pow(length/time,2);
}

Scales *normalizeSI(dictionary *ini){
	
	const double vacuumPermittivity = 8.854187817e-12; // [F/m]

	int nDims = iniGetInt(ini, "grid:nDims");
	int nSpecies = iniGetInt(ini, "population:nSpecies");
	double timeStep = iniGetDouble(ini, "time:timeStep");
	double *stepSize = iniGetDoubleArr(ini, "grid:stepSize", nDims);
	long int *nParticles = iniGetLongIntArr(ini, "population:nParticles", nSpecies);
	double *density = iniGetDoubleArr(ini, "population:density", nSpecies);
	double *charge = iniGetDoubleArr(ini, "population:charge", nSpecies);
	double *mass = iniGetDoubleArr(ini, "population:mass", nSpecies);
	double *thermalVelocity = iniGetDoubleArr(ini, "population:thermalVelocity", nSpecies);

	double V  = gGetGlobalVolume(ini)*pow(stepSize[0],nDims);

	double *K = (double*)malloc(nSpecies*sizeof(*K));
	for(int s=0; s<nSpecies; s++){
		K[s] = density[s]*V/nParticles[s];
	}

	double X  = stepSize[0];
	double T  = timeStep;
	double Q  = K[0]*fabs(charge[0]);
	double M  = pow(T*Q,2)/(vacuumPermittivity*pow(X,nDims));
	msg(STATUS, "Characteristic length: %15g m  (%a)",X,X);
	msg(STATUS, "Characteristic time  : %15g s  (%a)",T,T);
	msg(STATUS, "Characteristic charge: %15g C  (%a)",Q,Q);
	msg(STATUS, "Characteristic mass  : %15g kg (%a)",M,M);

	adMul(charge, K, charge, nSpecies);
	adMul(mass,   K, mass,   nSpecies);
	adScale(charge, nSpecies, 1.0/Q);
	adScale(mass,   nSpecies, 1.0/M);
	iniSetDoubleArr(ini, "population:charge", charge, nSpecies);
	iniSetDoubleArr(ini, "population:mass", mass, nSpecies);

	adScale(thermalVelocity, nSpecies, T/X);
	/* adScale(perturbAmplitude, nSpecies, 1/X); */

	iniSetDoubleArr(ini, "population:thermalVelocity", thermalVelocity, nSpecies);


	free(K);
	free(charge);
	free(mass);
	free(density);
	free(thermalVelocity);
	free(stepSize);
	free(nParticles);

	Scales *scales = malloc(sizeof(*scales));
	scales->length = X;
	scales->time = T;
	scales->charge = Q;
	scales->mass = M;
	nAddDerivedUnits(scales);
	return scales;
}

/******************************************************************************
 * DEFINING INI PARSING FUNCTIONS (expanding iniparser library)
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
			int nDims = iniGetInt(ini,"grid:nDims");
			int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
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
			msg(WARNING,"%s not specified. Using stdout.",keys[i]);
		} else if(strcmp(fName,"stdout")!=0 && strcmp(fName,"stderr")!=0) {

			if(makePath(fName))
				msg(ERROR,"Could not open or create path of '%s'",fName);

			// // check whether file exists
			// FILE *fh = fopen(fName,"r");
			// if(fh!=NULL){
			// 	fclose(fh);
			// 	msg(ERROR,"'%s' already exists.",fName);
			// }
		}

	}

	// Note that keys should be freed but not keys[0] and so on.
	free(keys);

	return ini;

}

void iniClose(dictionary *ini){
	iniparser_freedict(ini);
}

void iniAssertExistence(const dictionary *ini, const char* key){

	if(!iniparser_find_entry((dictionary*)ini,key)){
		msg(ERROR,"Key \"%s\" not found in input file",key);
	}
}

int iniGetNElements(const dictionary* ini, const char* key){

	iniAssertExistence(ini,key);
	char *list = iniparser_getstring((dictionary*)ini,key,"");	// don't free
	return listGetNElements(list);

}

int iniGetInt(const dictionary* ini, const char *key){

	iniAssertExistence(ini,key);
	char *res = iniparser_getstring((dictionary*)ini,key,0);	// don't free
	return (int)atof(res); // Parses scientific notation
}

long int iniGetLongInt(const dictionary* ini, const char *key){

	iniAssertExistence(ini,key);
	char *res = iniparser_getstring((dictionary*)ini,key,0);	// don't free
	return (long int)atof(res); // Parses scientific notation

}

double iniGetDouble(const dictionary* ini, const char *key){

	iniAssertExistence(ini,key);
	char *res = iniparser_getstring((dictionary*)ini,key,0);	// don't free
	return atof(res); // Parses scientific notation

}

char* iniGetStr(const dictionary *ini, const char *key){

	iniAssertExistence(ini,key);
	char *temp = iniparser_getstring((dictionary*)ini,key,NULL); // don't free

	// temp points to instance inside dictionary which is ruined if free'd.
	// Creating a copy which must be free'd to make this function's behavior
	// consistent with other ini-functions, e.g. iniGetIntArr().
	int len = strlen(temp);
	char *value = malloc((len+1)*sizeof(*value));
	strcpy(value,temp);

	return value;

}

int* iniGetIntArr(const dictionary *ini, const char *key, int nElements){

	char **strArr = iniGetStrArr(ini,key,nElements); // asserts existence

	int *result = malloc(nElements*sizeof(*result));
	for(int i=0;i<nElements;i++) result[i] = (int)atof(strArr[i]);

	freeStrArr(strArr);

	return result;

}

long int* iniGetLongIntArr(const dictionary *ini, const char *key, int nElements){

	char **strArr = iniGetStrArr(ini,key,nElements); // asserts existence

	long int *result = malloc(nElements*sizeof(*result));
	for(int i=0;i<nElements;i++) result[i] = (long int)atof(strArr[i]);

	freeStrArr(strArr);

	return result;

}

double* iniGetDoubleArr(const dictionary *ini, const char *key, int nElements){

	char **strArr = iniGetStrArr(ini,key,nElements); // asserts existence

	double *result = malloc(nElements*sizeof(double));
	for(int i=0;i<nElements;i++) result[i] = atof(strArr[i]);

	freeStrArr(strArr);

	return result;

}

char** iniGetStrArr(const dictionary *ini, const char *key, int nElements){

	iniAssertExistence(ini,key);

	char *list = iniparser_getstring((dictionary*)ini,key,"");	// don't free

	int nListElements = listGetNElements(list);

	if(nElements<nListElements){
		int nIgnored = nListElements - nElements;
		if(nIgnored==1){
			msg(WARNING, "Ignoring last element in %s",key);
		} else {
			msg(WARNING, "Ignoring last %d elements in %s",nIgnored,key);
		}
	}

	char **strArr = listToStrArr(list);
	char **strArrExpanded = strArrExpand(strArr,nElements);
	freeStrArr(strArr);

	return strArrExpanded;

}

void iniSetInt(dictionary *ini, const char *key, int value){
	iniSetLongInt(ini, key, value);
}

void iniSetLongInt(dictionary *ini, const char *key, long int value){

	const int numSize=32;
	char num[numSize];
	snprintf(num,numSize,"%ld",value);
	iniparser_set(ini, key, num);

}

void iniSetDouble(dictionary *ini, const char *key, double value){

	const int numSize=32;
	char num[numSize];
	snprintf(num,numSize,"%a",value);
	iniparser_set(ini, key, num);

}


void iniSetStr(dictionary *ini, const char *key, const char *value){
	iniparser_set(ini, key, value);
}

void iniSetIntArr(		dictionary *ini, const char *key,
						const int *values, int nElements){

	const int listSize=1024;
	const int numSize=32;
	char num[numSize];
	char list[listSize];

	num[0] = '\0';
	list[0] = '\0';
	for(int i=0; i<nElements; i++){
		snprintf(num,numSize,",%d",values[i]);
		strcat(list,num);
	}
	iniparser_set(ini, key, &list[1]);
}

void iniSetLongIntArr(	dictionary *ini, const char *key,
						const long int *values, int nElements){

	const int listSize=1024;
	const int numSize=32;
	char num[numSize];
	char list[listSize];

	num[0] = '\0';
	list[0] = '\0';
	for(int i=0; i<nElements; i++){
		snprintf(num,numSize,",%ld",values[i]);
		strcat(list,num);
	}
	iniparser_set(ini, key, &list[1]);
}

void iniSetDoubleArr(	dictionary *ini, const char *key,
						const double *values, int nElements){

	const int listSize=1024;
	const int numSize=32;
	char num[numSize];
	char list[listSize];

	num[0] = '\0';
	list[0] = '\0';
	for(int i=0; i<nElements; i++){
		snprintf(num,numSize,",%a",values[i]);
		strcat(list,num);
	}
	iniparser_set(ini, key, &list[1]);

}

void iniApplySuffix(dictionary *ini, const char *key,
					const char *suffix, const double *mul, int mulLen){

	int nElements = iniGetNElements(ini,key);
	char **strArr = iniGetStrArr(ini,key,nElements);

	const int listSize=1024;
	const int numSize=32;
	char num[numSize];
	char list[listSize];
	num[0] = '\0';
	list[0] = '\0';

	double *arr = malloc(nElements*sizeof(*arr));
	for(int i=0;i<nElements;i++){
		double val = atof(strArr[i]); // ignores suffix
		if(strstr(strArr[i],suffix)) val *= mul[i%mulLen];
		snprintf(num,numSize,",%a",val);
		strcat(list,num);
	}
	iniparser_set(ini,key,&list[1]);

	freeStrArr(strArr);

}

/******************************************************************************
 * DEFINING HDF5 FUNCTIONS (expanding HDF5 API)
 *****************************************************************************/

hid_t openH5File(const dictionary *ini, const char *fName, const char *fSubExt){

	// Determine filename
	char *fPrefix = iniGetStr(ini,"files:output");

	// Add separator if filename prefix (not just folder) is specified
	char sep[2] = "\0\0";
	char lastchar = fPrefix[strlen(fPrefix)-1];
	if(strcmp(fPrefix,".")==0) sep[0]='/';
	else if(strlen(fPrefix)>0 && lastchar!='/') sep[0]='_';

	char *fTotName = strCatAlloc(6,fPrefix,sep,fName,".",fSubExt,".h5");

	// Enable MPI-I/O access
	hid_t pList = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pList,MPI_COMM_WORLD,MPI_INFO_NULL);

	// Make sure parent folder exist
	if(makePath(fTotName))
		msg(ERROR,"Could not open or create folder for '%s'.",fTotName);

	hid_t file;	// h5 file handle

	// Open or create file (if it doesn't exist)
	FILE *fh = fopen(fTotName,"r");
	if(fh!=NULL){
		fclose(fh);
		file = H5Fopen(fTotName,H5F_ACC_RDWR,pList);
	} else {
		file = H5Fcreate(fTotName,H5F_ACC_EXCL,H5P_DEFAULT,pList);
	}

	H5Pclose(pList);
	free(fPrefix);
	free(fTotName);

	return file;

}

void setH5Attr(hid_t h5, const char *name, const double *value, int size){

	if(H5Aexists(h5,name)){

		int fNameSize = 128;
		char fName[fNameSize];
		H5Fget_name(h5,fName,fNameSize);

		msg(WARNING,"overwriting attribute \"%s\" in %s",name,fName);

		H5Adelete(h5,name);

	}

	// Create attribute dataspace
	hsize_t attrSize = (hsize_t)size;
	hid_t attrSpace = H5Screate_simple(1,&attrSize,NULL);

	// Create attribute and write data to it
	hid_t attribute = H5Acreate(h5, name, H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, value);

	// Free
    H5Aclose(attribute);
    H5Sclose(attrSpace);

}

void createH5Group(hid_t h5, const char *name){

	// Makes a new editable copy of name. Input may be string literal, which
	// cannot be edited anyway.
	char *str = malloc((strlen(name)+1)*sizeof(*str));
	strcpy(str,name);

	for(char *c=str+1; *c!='\0'; c++){

		if(*c=='/'){
			*c='\0';	// Temporarily ending string prematurely

			// Creates this part of the path if it doesn't already exist
			if(!H5Lexists(h5,str,H5P_DEFAULT)){
				hid_t group = H5Gcreate(h5,str,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
				H5Gclose(group);
			}

			*c='/';		// Changes string back
		}
	}

	free(str);

}


hid_t xyOpenH5(const dictionary *ini, const char *fName){

	return openH5File(ini,fName,"xy");
}

void xyCreateDataset(hid_t h5, const char *name){

	int mpiRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	createH5Group(h5,name);	// Creates parent groups

	const int arrSize=2;

	// Enable chunking of data in order to use extendible (unlimited) datasets
	hsize_t chunkDims[] = {1,2};
	hid_t pList = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(pList, arrSize, chunkDims);

	// Create dataspace for file initially empty but extendable
	hsize_t fileDims[] = {0,2};
	hsize_t fileDimsMax[] = {H5S_UNLIMITED,2};
	hid_t fileSpace = H5Screate_simple(arrSize,fileDims,fileDimsMax);

	// Create dataset in file using mentioned dataspace
	hid_t dataset = H5Dcreate(h5,name,H5T_IEEE_F64LE,fileSpace,H5P_DEFAULT,pList,H5P_DEFAULT);

	H5Sclose(fileSpace);
	H5Dclose(dataset);

}
void xyWrite(hid_t h5, const char* name, double x, double y, MPI_Op op){

	int mpiRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);

	// Reduce data across nodes
	double yReduced;
	MPI_Reduce(&y,&yReduced,1,MPI_DOUBLE,op,0,MPI_COMM_WORLD);

	// Load dataset
	hid_t dataset = H5Dopen(h5,name,H5P_DEFAULT);

	// Extend dataspace in file by one row (must be done on all MPI nodes)
	const int arrSize=2;
	hid_t fileSpace = H5Dget_space(dataset);
	hsize_t fileDims[arrSize];
	H5Sget_simple_extent_dims(fileSpace,fileDims,NULL);
	fileDims[0]++;
	H5Dset_extent(dataset,fileDims);

	// update fileSpace after change
	H5Sclose(fileSpace);
	fileSpace = H5Dget_space(dataset);

	// Write only from MPI rank 0
	if(mpiRank==0){
		// Select hyperslab to write to
		hsize_t offset[] = {fileDims[0]-1,0};
		hsize_t count[] = {1,1};
		hsize_t memDims[] = {1,2};
		H5Sselect_hyperslab(fileSpace,H5S_SELECT_SET,offset,NULL,count,memDims);

		// Write to file
		double data[] = {x,yReduced};
		hid_t memSpace = H5Screate_simple(arrSize,memDims,NULL);
		H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, H5P_DEFAULT, data);
		H5Sclose(memSpace);
	}

	H5Sclose(fileSpace);
	H5Dclose(dataset);

}

void xyCloseH5(hid_t h5){

	H5Fclose(h5);
}

/******************************************************************************
 * DEFINING LOCAL LIST PARSING FUNCTIONS
 *****************************************************************************/

static int listGetNElements(const char* list){

	if(list[0]=='\0') return 0;

	// Count elements
	int nElements = 1;
	for(int i=0;list[i];i++) nElements += (list[i]==',');

	return nElements;

}

static char** listToStrArr(const char* restrict list){

	// Count elements in list
	int nElements = 2;				// first element and NULL
	char *temp = (char*)list;
	while(*temp){
		if(*temp==',') nElements++;
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

int strArrLen(char **strArr){

	int nElements=0;
	for(; *strArr != NULL; strArr++, nElements++){}
	return nElements;

}

char **strArrExpand(char **strArr, int nElements){

	int nElementsOld = strArrLen(strArr);

	char **result = malloc((nElements+1)*sizeof(result));

	for(int i=0;i<nElements;i++){

		int len = strlen(strArr[i%nElementsOld]);
		result[i] = malloc((len+1)*sizeof(char));
		strcpy(result[i],strArr[i%nElementsOld]);

	}

	result[nElements]=NULL;

	return result;

}

void freeStrArr(char** strArr){

	for(int i=0;strArr[i];i++) free(strArr[i]);
	free(strArr);

}

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
