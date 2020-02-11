/**
 * @file		aux.c
 * @brief		Small auxiliary functions
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 *
 * Small auxiliary functions.
 */
#define _XOPEN_SOURCE 700

#include "core.h"
#include <time.h>
#include <mpi.h>
#include <stdarg.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/******************************************************************************
 * LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/

#if _POSIX_TIMERS>0 && defined(_POSIX_MONOTONIC_CLOCK) // Linux

	unsigned long long int getNanoSec(void){

		struct timespec time;
		clock_gettime(CLOCK_MONOTONIC, &time);

		return time.tv_sec*1e9 + time.tv_nsec;
	}

#elif defined(__APPLE__) // Mac

	#include <mach/mach_time.h>

    unsigned long long int getNanoSec(void){

		// needed to convert mach_absolute time to something meaningful
	    mach_timebase_info_data_t timebase;
	    mach_timebase_info(&timebase);

		uint64_t clock = mach_absolute_time();
		uint64_t nanosecs = clock * (uint64_t)timebase.numer / (uint64_t)timebase.denom;

		return nanosecs;
	}

#endif

Timer *tAlloc(){

	Timer *t = malloc(sizeof(*t));
	t->total = 0;

	return t;
}

void tFree(Timer *t){
	free(t);
}


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
		msg(TIMER, "%s %6.2fs ", string, (double) nanoSec/1e9);
	} else if(nanoSec>1e6) {
		msg(TIMER, "%s %6.2fms ", string, (double) nanoSec/1e6);
	} else if(nanoSec>1e3) {
		msg(TIMER, "%s %6.2fus ", string, (double) nanoSec/1e3);
	} else {
		msg(TIMER, "%s %6.2fns ", string, (double) nanoSec);
	}

}

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

void adSub(const double *a, const double *b, double *res, long int n){
	for(long int i=0;i<n;i++) res[i] = a[i] - b[i];
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

void adScale(double *a, long int n, double value){
	for(long int i=0;i<n;i++) a[i] *= value;
}

void aiScale(int *a, long int n, int value){
	for(long int i=0;i<n;i++) a[i] *= value;
}

void alScale(long int *a, long int n, int value){
	for(long int i=0;i<n;i++) a[i] *= value;
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

void adCrossProd(const double *a, const double *b, double *res){
	
	res[0] +=  (a[1]*b[2]-a[2]*b[1]);
	res[1] += -(a[0]*b[2]-a[2]*b[0]);
	res[2] +=  (a[0]*b[1]-a[1]*b[0]);
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

void adNormal(const double *a, const double *b, double *res, long int n){
    
	double length = 1.;

    adCrossProd(a, b, res);
    length = adDotProd(res, res,3);
    length = sqrt(length);

    for(int i=0;i<2;i++) res[i] /= length;
}

//w = v - 2(v . n)n where n=norm, v=ray, w=res, surface points a and b (x,y,z) locations
void adReflect(const double *ray, const double *a, const double *b, double *res){
    
	double *norm = NULL;
	adNormal(a,b,norm,3);
    int ray_dot_n = adDotProd(ray,norm,3);

    for(int i=0;i<3;i++) res[i] = ray[i] - 2 * ray_dot_n * norm[i];
}

void adRotateRandom3D(double *a, const gsl_rng *r){

	double beta = PI/2 * gsl_rng_uniform(r);
	double gamma = PI/2 * gsl_rng_uniform(r);
	int signSB = (gsl_rng_uniform_int(r, 2) == 0) ? -1 : 1;
	int signSG = (gsl_rng_uniform_int(r, 2) == 0) ? -1 : 1;
	int signCB = (gsl_rng_uniform_int(r, 2) == 0) ? -1 : 1;
	int signCG = (gsl_rng_uniform_int(r, 2) == 0) ? -1 : 1;


	//compute sin,cos once
	const double sinB = sin(signSB * beta);
	const double sinG = sin(signSG * gamma);
	const double cosB = cos(signCB * beta);
	const double cosG = cos(signCG * gamma); 

	//msg(STATUS, "cosB: %f, sinB: %f, cosG: %f, sinG: %f", cosB, sinB, cosG, sinG);
	gsl_vector_view aVec = gsl_vector_view_array(a, 3);
	
	double Ry[] = {cosB, 0.0, sinB,
				   0.0, 1.0, 0.0,
				   -sinB,0.0, cosB};

	double Rz[] = {cosG, -sinG, 0.0,
				   sinG, cosG, 0.0,
				   0.0, 0.0, 1.0};
				
	
	gsl_matrix_view Ry_mat = gsl_matrix_view_array(Ry, 3, 3);
	gsl_matrix_view Rz_mat = gsl_matrix_view_array(Rz, 3, 3);
	gsl_matrix *R_tot = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_all(R_tot, 0.0);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &Ry_mat.matrix, &Rz_mat.matrix,
				   0.0, R_tot);
	gsl_vector *rotatedVec = gsl_vector_alloc(3);
	gsl_vector_set_all(rotatedVec, 0.0);
	gsl_blas_dgemv(CblasNoTrans, 1.0, R_tot, &aVec.vector, 0.0, rotatedVec);

	//printf("==================== \n");
	//msg(STATUS, "Before rotation");
	//adPrint(a, 3);
	for(int i = 0; i < 3; i++){
		a[i] = gsl_vector_get(rotatedVec, i);
	}

	//msg(STATUS, "After rotation");
	//adPrint(a,3);
	//printf("==================== \n");

	gsl_vector_free(rotatedVec);
	gsl_matrix_free(R_tot);

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

void adCumSum(const double *a, double *res, long int n){
	res[0]=0;
	for(long int i=0;i<n;i++) res[i+1]=res[i]+a[i];
}

void aiCumSum(const int *a, int *res, long int n){
	res[0]=0;
	for(long int i=0;i<n;i++) res[i+1]=res[i]+a[i];
}

void ailCumSum(const int *a, long int *res, long int n){
	res[0]=0;
	for(long int i=0;i<n;i++) res[i+1]=res[i]+a[i];
}

void alCumSum(const long int *a, long int *res, long int n){
	res[0]=0;
	for(long int i=0;i<n;i++) res[i+1]=res[i]+a[i];
}

void bSetAll(bool *a, long int n, bool value){
	for(long int i=0;i<n;i++) a[i] = value;
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

void adCopy(double *a, double *b, long int n){
	for(long int i=0; i<n; i++) b[i] = a[i];
}

void alCopy(long int *a, long int *b, long int n){
	for(long int i=0; i<n; i++) b[i] = a[i];
}


void adPrintInner(double *a, long int inc, long int end, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s(1:%li:%li) = \n  [",rank,varName,inc,end);
	int i;
	for(i=0;i<end-inc;i+=inc){
		printf("%g ",a[i]);
	}

	printf("%g]\n",a[i]);

}

void aiPrintInner(int *a, long int inc, long int end, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s(1:%li:%li) = \n  [",rank,varName,inc,end);
	int i;
	for(i=0;i<end-1;i+=inc){
		printf("%i ",a[i]);
	}
	printf("%i]\n",a[i]);

}

void alPrintInner(long int *a, long int inc, long int end, char *varName){

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	printf("PRINT(%i): %s(1:%li:%li) = \n  [",rank,varName,inc,end);
	int i;
	for(i=0;i<end-1;i+=inc){
		printf("%li ",a[i]);
	}
	printf("%li]\n",a[i]);

}

/******************************************************************************
 * DEBUG HELP
 *****************************************************************************/

 void fillGridIndexes(Grid *grid){

 	msg(WARNING, "Grid is replaced with indexes");

 	long int *sizeProd = grid->sizeProd;
 	int rank = grid->rank;
 	double *val = grid->val;

 	for(long int g = 0; g < sizeProd[rank]; g++){
 		val[g] = (double) g;
 	}

 	return;
 }

void dumpTrueGrid(dictionary *ini, Grid *grid){

	int *size = grid->size;
	long int *sizeProd = grid->sizeProd;
	int nDims = grid->rank -1;

	msg(STATUS, "Dumps grid to parsefile");
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

    msg(STATUS, "Dumps grid to parsefile");
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