/**
 * @file	    test.c
 * @brief	    PINC unit testing framework
 * @author	    Sigvald Marholm <sigvaldm@fys.uio.no>
 */

#include "test.h"
#include <stdarg.h>
#include <stdio.h>

#define BUFFSIZE 256

int utAssertInner(const char *file, const char *func, int line, int test, const char * restrict format,...){

	if(test) return 0;
	else {

		va_list args;
		va_start(args,format);

		char buffer[BUFFSIZE];
		vsnprintf(buffer,BUFFSIZE,format,args);
		fprintf(stderr,"%s:%i %s: %s\n",file,line,func,buffer);

		va_end(args);

		return 1;

	}
}

void utRun(int (*fun)()){

	static int nTrials = 0;
	static int nSuccess = 0;

	if(fun==NULL){
		if(nSuccess==nTrials)
			printf("SUMMARY: All %i tests passed.\n",nTrials);
		else
			printf("SUMMARY: Passed %i of %i tests. %i failed.\n",nSuccess,nTrials,nTrials-nSuccess);
	} else {
		nTrials++;
		nSuccess += 1-(*fun)();
	}

}

void utSummary(){
	utRun(NULL);
}

dictionary *iniSetDummy(int argc, char **argv){

	static char fileName[128];
	if(argc==0){
		return iniparser_load(fileName);
	} else if(argc==2){
		strcpy(fileName,argv[1]);
		return NULL;
	} else {
		printf("error in iniSetDummy");
		return NULL;
	}
}

dictionary *iniGetDummy(){
	return iniSetDummy(0,0);
}
