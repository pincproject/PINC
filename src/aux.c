/**
 * @file		aux.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Small auxiliary functions
 * @date		29.10.15
 *
 * Small auxiliary functions.
 */

#define _XOPEN_SOURCE 700

#include "pinc.h"
#include <time.h>
#include <mpi.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

/******************************************************************************
 * LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/

static void tFormat(char *str, int len, const TimeSpec *time);

/**
 * @brief Makes a directory
 * @param	dir		Directory name
 * @return	0 for success, 1 for failure
 *
 * dir can be a path but ancestors must exist. Directory will have permissions
 * 0775.
 */
static int makeDir(const char *dir);

/******************************************************************************
 * ARRAY FUNCTIONS
 *****************************************************************************/

int *aiMul(const int *a, const int *b, int nElements){

	int *result = malloc(nElements*sizeof(int));

	for(int i=0;i<nElements;i++){
		result[i]=a[i]*b[i];
	}

	return result;

}

// Shift all elements by value (aiAdd is reserved for adding to arrays)
void aiShift(int *a, long int n, int value){

	for(long int i=0;i<n;i++){
		a[i] += value;
	}

}
void alShift(long int *a, long int n, long int value){

	for(long int i=0;i<n;i++){
		a[i] += value;
	}

}

int aiMax(const int *a, int nElements){

	int max = a[0];

	for(int i=1;i<nElements;i++){
		if(a[i] > max) max = a[i];
	}

	return max;

}

long int alMax(const long int *a, int nElements){

	int max = a[0];

	for(int i=1;i<nElements;i++){
		if(a[i] > max) max = a[i];
	}

	return max;

}

int *aiCumProd(const int *a, int nElements){

	int *result = malloc((nElements+1)*sizeof(*result));
	result[0]=1;

	for(int i=1;i<nElements+1;i++){
		result[i]=result[i-1]*a[i-1];
	}

	return result;

}

long int *ailCumProd(const int *a, int nElements){

	long int *result = malloc((nElements+1)*sizeof(*result));
	result[0]=1;

	for(int i=1;i<nElements+1;i++){
		result[i]=result[i-1]*a[i-1];
	}

	return result;

}

int aiProd(const int *a, int nElements){

	int result = 1;
	for(int i=0;i<nElements;i++){
		result *= a[i];
	}

	return result;
}

int aiDotProd(const int *a, const int *b, int nElements){

	int result = 0;
	for(int i=0;i<nElements;i++){
		result += a[i]*b[i];
	}
	return result;
}

double adProd(const double *a, int nElements){

	double result = 1;
	for(int i=0;i<nElements;i++){
		result *= a[i];
	}

	return result;

}

int aiEq(const int *a, const int *b, long int nElements){
	int equal = 1;
	for(long int i=0;i<nElements;i++) if(a[i]!=b[i]) equal=0;
	return equal;
}

int alEq(const long int *a, const long int *b, long int nElements){
	int equal = 1;
	for(long int i=0;i<nElements;i++) if(a[i]!=b[i]) equal=0;
	return equal;
}

int aEq(const char *a, const char *b, long int nBytes){
	int equal = 1;
	for(long int i=0;i<nBytes;i++) if(a[i]!=b[i]) equal=0;
	return equal;
}

int adEq(const double *a, const double *b, long int nElements, double tol){
	for(long int i=0;i<nElements;i++){
		if(fabs(a[i]-b[i])>tol) return 0;
	}
	return 1;
}

void adSet(double *a, long int n, ...){

	va_list args;
	va_start(args,n);

	for(int i=0;i<n;i++){
		double arg = va_arg(args,double);
		a[i] = arg;
	}

	va_end(args);

}

void aiSet(int *a, long int n, ...){

	va_list args;
	va_start(args,n);

	for(long int i=0;i<n;i++){
		int arg = va_arg(args,int);
		a[i] = arg;
	}

	va_end(args);

}

void alSet(long int *a, long int n, ...){

	va_list args;
	va_start(args,n);

	for(long int i=0;i<n;i++){
		int arg = va_arg(args,long int);
		a[i] = arg;
	}

	va_end(args);

}

void adSetAll(double *a, long int n, double value){

	for(long int i=0;i<n;i++) a[i] = value;

}

void alSetAll(long int *a, long int n, long int value){

	for(long int i=0;i<n;i++) a[i] = value;

}

void adPrintInner(double *a, long int n, char *varName){

	printf("%s = \n  [",varName);
	for(int i=0;i<n-1;i++){
		printf("%f ",a[i]);
	}
	printf("%f]\n",a[n-1]);

}

void aiPrintInner(int *a, long int n, char *varName){

	printf("%s = \n  [",varName);
	for(int i=0;i<n-1;i++){
		printf("%i ",a[i]);
	}
	printf("%i]\n",a[n-1]);

}

void alPrintInner(long int *a, long int n, char *varName){

	printf("%s = \n  [",varName);
	for(int i=0;i<n-1;i++){
		printf("%li ",a[i]);
	}
	printf("%li]\n",a[n-1]);

}

/******************************************************************************
 * TIMING FUNCTIONS
 *****************************************************************************/

static void tFormat(char *str, int len, const TimeSpec *time){

	long int sec = time->tv_sec;
	long int nsec = time->tv_nsec;

	if(sec>=1){
		snprintf(str,len,"%6.2fs ",(double)sec+(double)nsec/1000000000);
	} else if(nsec>1000000) {
		snprintf(str,len,"%6.2fms",(double)nsec/1000000);
	} else if(nsec>1000) {
		snprintf(str,len,"%6.2fus",(double)nsec/1000);
	} else {
		snprintf(str,len,"%6.2fns",(double)nsec);
	}
}

void tMsg(Timer *timer, const char *restrict format, ...){

	if(timer->rank>=0){
		TimeSpec now, diff;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);

		if(format!=NULL){

			// Take difference
			diff.tv_sec = now.tv_sec-timer->previous.tv_sec;
			diff.tv_nsec = now.tv_nsec-timer->previous.tv_nsec;

			// Borrow from tv_sec if necessary
			if(diff.tv_nsec<0){
				diff.tv_sec-=1;
				diff.tv_nsec+=1e9;
			}

			// Format time
			const int strSize = 16;
			char nowStr[strSize];
			char diffStr[strSize];
			tFormat(nowStr,strSize,&now);
			tFormat(diffStr,strSize,&diff);

			// Get message
			const int bufferSize = 132;
			char msg[bufferSize];
			va_list args;
			va_start(args,format);
			vsnprintf(msg,bufferSize,format,args);
			va_end(args);

			// Print
			char buffer[bufferSize];
			snprintf(buffer,bufferSize,"TIMER (%i): [tot=%s, diff=%s] %s",timer->rank,nowStr,diffStr,msg);
			fprintf(stdout,"%s\n",buffer);
		}

		// Store time of this call
		//timer->previous = now;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timer->previous);
	}
}

Timer *tAlloc(int rank){

	// Get rank of this MPI node
	int thisRank;
	MPI_Comm_rank(MPI_COMM_WORLD,&thisRank);

	// Assign struct
	Timer *timer = malloc(sizeof(Timer));
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timer->previous);
	if(rank==thisRank || rank<0){
		timer->rank=thisRank;
	} else {
		timer->rank=-1;	// Deactivate this timer
	}

	return timer;
}

void tFree(Timer *timer){
	free(timer);
}

char *strAllocCat(int n, ...){

	va_list args;
	va_start(args,n);

	int len = 0;

	for(int i=0;i<n;i++){
		char *s = va_arg(args,char*);
		len += strlen(s);
	}

	va_end(args);

	char *result = malloc(len+1);	// Remember '\0'

	va_start(args,n);

	char *s = va_arg(args,char*);
	strcpy(result,s);

	for(int i=1;i<n;i++){
		s = va_arg(args,char*);
		strcat(result,s);
	}

	return result;

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

int makeParentPath(const char *path){
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
