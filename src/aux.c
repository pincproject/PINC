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
#include <math.h>

/******************************************************************************
 * LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/
//
// /**
//  * @brief Automatically formats a TimeSpec into a string
//  * @param[out]		str		Formatted string
//  * @param			len		Allocated length in str
//  * @param			time	TimeSpec to format
//  *
//  * Used in tMsg().
//  */
// static void tFormat(char *str, int len, const TimeSpec *time);
//
// /******************************************************************************
//  * TIMING FUNCTIONS
//  *****************************************************************************/

Timer *tAlloc(){

	Timer *t = malloc(sizeof(*t));
	t->total = 0;

	return t;
}

void tFree(Timer *t){
	free(t);
}

//Linux or mac implementation to getNanoSec (Local function)
#if _POSIX_TIMERS>0 && defined(_POSIX_MONOTONIC_CLOCK)

	unsigned long long int getNanoSec(void){

		struct timespec time;
		clock_gettime(CLOCK_MONOTONIC, &time);

		return time.tv_sec*1e9 + time.tv_nsec;
	}

#elif	defined(__APPLE__)

         unsigned long long int getNanoSec(void){
		//TBD
		//msg(WARNING, "Timer function for macOS is not completed.");

                #include <mach/mach_time.h> // mach_absolute_time

	   //#define BILLION 1000000000L

	   //struct timespec time; // represent the elapsed time, 
	        mach_timebase_info_data_t timebase; // needed to convert mach_absolute time to something meaningful  
	        mach_timebase_info(&timebase); // needed to convert mach_absolute time to something meaningful

		uint64_t clock = mach_absolute_time(); // capture the current time
		uint64_t nanosecs = clock * (uint64_t)timebase.numer / (uint64_t)timebase.denom; // convert to nanoseconds

		//time.tv_sec = nanosecs / BILLION; // convert to timespec
		//time.tv_nsec = nanosecs % BILLION; // convert to timespec
	  
		return nanosecs; //time.tv_sec*1e9 + time.tv_nsec;
	    
	}

#endif


void tStart(Timer *t){
	t->start = getNanoSec();
}

void tStop(Timer *t){
	t->total += getNanoSec() - t->start;
}

void tReset(Timer *t){
	t->total = 0;
}

void tMsg(long long int nanoSec, const char *string){

	if(nanoSec >= 1e9){
		msg(TIMER|ONCE, "%s %6.2fs ", string, (double) nanoSec/1e9);
	} else if(nanoSec>1e6) {
		msg(TIMER|ONCE, "%s %6.2fms ", string, (double) nanoSec/1e6);
	} else if(nanoSec>1e3) {
		msg(TIMER|ONCE, "%s %6.2fus ", string, (double) nanoSec/1e3);
	} else {
		msg(TIMER|ONCE, "%s %6.2fns ", string, (double) nanoSec);
	}

}


/****************************************************************************
 *	Outdated Timer functions (Some functionality is lost in new)
 ****************************************************************************/


// void tMsg(Timer *timer, const char *restrict format, ...){
//
// 	if(timer->rank>=0){
// 		TimeSpec now, diff;
// 		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);
//
// 		if(format!=NULL){
//
// 			// Take difference
// 			diff.tv_sec = now.tv_sec-timer->previous.tv_sec;
// 			diff.tv_nsec = now.tv_nsec-timer->previous.tv_nsec;
//
// 			// Borrow from tv_sec if necessary
// 			if(diff.tv_nsec<0){
// 				diff.tv_sec-=1;
// 				diff.tv_nsec+=1e9;
// 			}
//
// 			// Format time
// 			const int strSize = 16;
// 			char nowStr[strSize];
// 			char diffStr[strSize];
// 			tFormat(nowStr,strSize,&now);
// 			tFormat(diffStr,strSize,&diff);
//
// 			// Get message
// 			const int bufferSize = 132;
// 			char msg[bufferSize];
// 			va_list args;
// 			va_start(args,format);
// 			vsnprintf(msg,bufferSize,format,args);
// 			va_end(args);
//
// 			// Print
// 			char buffer[bufferSize];
// 			snprintf(buffer,bufferSize,"TIMER (%i): [tot=%s, diff=%s] %s",timer->rank,nowStr,diffStr,msg);
// 			fprintf(stdout,"%s\n",buffer);
// 		}
//
// 		// Store time of this call
// 		//timer->previous = now;
// 		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timer->previous);
// 	}
// }

/******************************************************************************
 * STRING FUNCTIONS
 *****************************************************************************/

char *strCatAlloc(int n, ...){

	va_list args;
	va_start(args,n);

	int len = 0;

	for(int i=0;i<n;i++){
		char *s = va_arg(args,char*);
		len += strlen(s);
	}

	va_end(args);

	char *res = malloc(len+1);	// Remember '\0'

	va_start(args,n);

	char *s = va_arg(args,char*);
	strcpy(res,s);

	for(int i=1;i<n;i++){
		s = va_arg(args,char*);
		strcat(res,s);
	}

	return res;

}

/******************************************************************************
 * ARRAY FUNCTIONS
 *****************************************************************************/

void adAdd(const double *a, const double *b, double *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]+b[i];
}

void aiAdd(const int *a, const int *b, int *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]+b[i];
}

void alAdd(const long int *a, const long int *b, long int *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]+b[i];
}

void adMul(const double *a, const double *b, double *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]*b[i];
}

void aiMul(const int *a, const int *b, int *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]*b[i];
}

void alMul(const long int *a, const long int *b, long int *res, long int n){
	for(long int i=0;i<n;i++) res[i]=a[i]*b[i];
}

void adShift(double *a, long int n, double value){
	for(long int i=0;i<n;i++) a[i] += value;
}

void aiShift(int *a, long int n, int value){
	for(long int i=0;i<n;i++) a[i] += value;
}

void alShift(long int *a, long int n, long int value){
	for(long int i=0;i<n;i++) a[i] += value;
}

double adMax(const double *a, long int n){
	double res = a[0];
	for(long int i=1;i<n;i++) if(a[i] > res) res = a[i];
	return res;
}

int aiMax(const int *a, long int n){
	int res = a[0];
	for(long int i=1;i<n;i++) if(a[i] > res) res = a[i];
	return res;
}

long int alMax(const long int *a, long int n){
	long int res = a[0];
	for(long int i=1;i<n;i++) if(a[i] > res) res = a[i];
	return res;
}

double adMin(const double *a, long int n){
	double res = a[0];
	for(long int i=1;i<n;i++) if(a[i] < res) res = a[i];
	return res;
}

int aiMin(const int *a, long int n){
	int res = a[0];
	for(long int i=1;i<n;i++) if(a[i] < res) res = a[i];
	return res;
}

long int alMin(const long int *a, long int n){
	long int res = a[0];
	for(long int i=1;i<n;i++) if(a[i] < res) res = a[i];
	return res;
}

double adExt(const double *a, long int n){
	double min = adMin(a,n);
	double max = adMax(a,n);
	return (max>-min)?max:min;
}

int aiExt(const int *a, long int n){
	int min = aiMin(a,n);
	int max = aiMax(a,n);
	return (max>-min)?max:min;
}

long int alExt(const long int *a, long int n){
	long int min = alMin(a,n);
	long int max = alMax(a,n);
	return (max>-min)?max:min;
}

double adSum(const double *a, long int n){
	double sum = 0;
	for(long int i=0;i<n;i++) sum += a[i];
	return sum;
}

long int aiSum(const int *a, long int n){
	int sum = 0;
	for(long int i=0;i<n;i++) sum += a[i];
	return sum;
}

long int alSum(const long int *a, long int n){
	long int sum = 0;
	for(long int i=0;i<n;i++) sum += a[i];
	return sum;
}

double adAvg(const double *a, long int n){
	return (double)adSum(a,n)/n;
}

double aiAvg(const int *a, long int n){
	return (double)aiSum(a,n)/n;
}

double alAvg(const long int *a, long int n){
	return (double)alSum(a,n)/n;
}

double adProd(const double *a, long int n){
	double res = 1;
	for(long int i=0;i<n;i++) res *= a[i];
	return res;
}

long int aiProd(const int *a, long int n){
	int res = 1;
	for(long int i=0;i<n;i++) res *= a[i];
	return res;
}

long int alProd(const int *a, long int n){
	long int res = 1;
	for(long int i=0;i<n;i++) res *= a[i];
	return res;
}

int adDotProd(const double *a, const double *b, long int n){
	double res = 0;
	for(long int i=0;i<n;i++) res += a[i]*b[i];
	return res;
}

int aiDotProd(const int *a, const int *b, long int n){
	int res = 0;
	for(long int i=0;i<n;i++) res += a[i]*b[i];
	return res;
}

int alDotProd(const long int *a, const long int *b, long int n){
	long int res = 0;
	for(long int i=0;i<n;i++) res += a[i]*b[i];
	return res;
}

int adEq(const double *a, const double *b, long int n, double tol){
	for(long int i=0;i<n;i++) if(fabs(a[i]-b[i])>tol) return 0;
	return 1;
}

int aiEq(const int *a, const int *b, long int n){
	for(long int i=0;i<n;i++) if(a[i]!=b[i]) return 0;
	return 1;
}

int alEq(const long int *a, const long int *b, long int n){
	for(long int i=0;i<n;i++) if(a[i]!=b[i]) return 0;
	return 1;
}

void adCumProd(const double *a, double *res, long int n){
	res[0]=1;
	for(long int i=0;i<n;i++) res[i+1]=res[i]*a[i];
}

void aiCumProd(const int *a, int *res, long int n){
	res[0]=1;
	for(long int i=0;i<n;i++) res[i+1]=res[i]*a[i];
}

void ailCumProd(const int *a, long int *res, long int n){
	res[0]=1;
	for(long int i=0;i<n;i++) res[i+1]=res[i]*a[i];
}

void alCumProd(const long int *a, long int *res, long int n){
	res[0]=1;
	for(long int i=0;i<n;i++) res[i+1]=res[i]*a[i];
}

void adSetAll(double *a, long int n, double value){
	for(long int i=0;i<n;i++) a[i] = value;
}

void aiSetAll(int *a, long int n, int value){
	for(long int i=0;i<n;i++) a[i] = value;
}

void alSetAll(long int *a, long int n, long int value){
	for(long int i=0;i<n;i++) a[i] = value;
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

void adPrintInner(double *a, long int n, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s = \n  [",rank,varName);
	for(int i=0;i<n-1;i++){
		printf("%f ",a[i]);
	}
	printf("%f]\n",a[n-1]);

}

void aiPrintInner(int *a, long int n, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s = \n  [",rank,varName);
	for(int i=0;i<n-1;i++){
		printf("%i ",a[i]);
	}
	printf("%i]\n",a[n-1]);

}

void alPrintInner(long int *a, long int n, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s = \n  [",rank,varName);
	for(int i=0;i<n-1;i++){
		printf("%li ",a[i]);
	}
	printf("%li]\n",a[n-1]);

}

/***************************************************************************
 *			Debug help
 ***************************************************************************/

 void dumpTrueGrid(dictionary *ini, Grid *grid){

 	int *size = grid->size;
 	long int *sizeProd = grid->sizeProd;
 	int nDims = grid->rank -1;

 	msg(STATUS|ONCE, "Dumps grid to parsefile");
 	if(nDims == 3){
 		fMsg(ini,"parsedump", "\nDump of 3D grid: (%dx%dx%d) \n \n",
 		 			size[1], size[2], size[3]);
 		//Cycles trough and prints the grid (not optimized)
 		int p;
 		for(int l = 1; l < size[3]-1; l++){
 			fMsg(ini, "parsedump", "\t\t\t l = %d \n", l);
 			for(int k = size[2] - 2; k > 0; k--){ //y-rows
 				for(int j = 1; j < size[1]-1; j++){ //x-rows
 					p = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
 					fMsg(ini,"parsedump", "%3.1f, ",  grid->val[p]);
 				}
 				fMsg(ini,"parsedump", "\n\n");
 			}
 		}
 	} else if(nDims==2) {
 		fMsg(ini,"parsedump", "2D grid: (%dx%d): \n",
 		 			size[1], size[2]);
 		int p;
 		for(int k = size[2] - 2; k > 0; k--){ //y-rows
 			for(int j = 1; j < size[1]-1; j++){ //x-rows
 				p = j*size[0] + k*size[1];
 				fMsg(ini,"parsedump", "%3.1f \t", grid->val[p]);
 			}
 			fMsg(ini,"parsedump", "\n");
 		}

 	}

 	return;
 }

 void dumpWholeGrid(dictionary *ini, Grid *grid){

    int *size = grid->size;
    long int *sizeProd = grid->sizeProd;
    int nDims = grid->rank -1;

    msg(STATUS|ONCE, "Dumps grid to parsefile");
    if(nDims == 3){
 	   fMsg(ini,"parsedump", "\nDump of 3D grid: (%dx%dx%d) \n \n",
 				   size[1], size[2], size[3]);
 	   //Cycles trough and prints the grid (not optimized)
 	   int p;
 	   for(int l = 0; l < size[3]; l++){
 		   fMsg(ini, "parsedump", "\t\t\t l = %d \n", l);
 		   for(int k = size[2] - 1; k > -1; k--){ //y-rows
 			   for(int j = 0; j < size[1]; j++){ //x-rows
 				   p = j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
 				   fMsg(ini,"parsedump", "%3.1f, ",  grid->val[p]);
 			   }
 			   fMsg(ini,"parsedump", "\n\n");
 		   }
 	   }
    } else if(nDims==2) {
 	   fMsg(ini,"parsedump", "2D grid: (%dx%d): \n",
 				   size[1], size[2]);
 	   int p;
 	   for(int k = size[2] - 1; k > -1; k--){ //y-rows
 		   for(int j = 0; j < size[1]; j++){ //x-rows
 			   p = j*size[0] + k*size[1];
 			   fMsg(ini,"parsedump", "%3.1f \t", grid->val[p]);
 		   }
 		   fMsg(ini,"parsedump", "\n");
 	   }

    }

    return;
 }
