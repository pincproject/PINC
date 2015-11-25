#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * TIMING FUNCTION
 */

 // Time probe. See example above.
 void tprobe(const char * label){

 	static struct timespec previous;	// time of previous tprobe call

 	struct timespec now, diff;
 	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &now);

	if(label==NULL){
		// RESET TIMER
	 	previous = now;
	} else {

	 	// Take difference
	 	diff.tv_sec = now.tv_sec-previous.tv_sec;
	 	diff.tv_nsec = now.tv_nsec-previous.tv_nsec;

	 	// Borrow from tv_sec if necessary
	 	if(diff.tv_nsec<0){
	 		diff.tv_sec-=1;
	 		diff.tv_nsec+=1e9;
	 	}

		long int sec = diff.tv_sec;
		long int nsec = diff.tv_nsec;

	 	// Print
	 	printf("Time: ");
		if(sec>=1){
	 		printf("%.2fs",(double)sec+(double)nsec/1000000000);
	 	} else if(nsec>1000000) {
	 		printf("%.1fms",(double)nsec/1000000);
	 	} else if(nsec>1000) {
	 		printf("%.1fus",(double)nsec/1000);
	 	} else {
	 		printf("%ldns",nsec);
		}
	 	printf(" %s\n",label);

	}

 }

 void printSP(){
 	register const long rsp __asm__ ("rsp");
 	printf("Stack pointer: %p\n", (void*)rsp);
 }

/*
 * METHOD WITH PRE-STORED INDICES
 */

void *getSliceIndices(long int *indices, const int *nGPoints, const long int *nGPointsProd, int nDims, int dim, int offset){
 	for(long int p=0; p<nGPointsProd[nDims]; p++){

 		long int index = p;
 		for(int d=0; d<dim; d++) index /= nGPoints[d];
		index %= nGPoints[dim];

		if(index==offset) *(indices++) = p;

 	}
}

void getFromInd(double *slice, double *val, long int *indices, int nElements){
	for(long int p=0;p<nElements;p++){
		slice[p] = val[indices[p]];
	}
}

/*
 * DIRECT METHOD
 */

double *getSliceInner(double *nextGhost, const double **valp, const long int *mul, const int *points, const long int finalMul){

//	printSP();

	if(*mul==finalMul){
		for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
		*valp += (*mul)*(*points-1);
	} else {
		for(int j=0;j<*points;j++) nextGhost = getSliceInner(nextGhost,valp,mul-1,points-1,finalMul);
	}
	return nextGhost;
}

void getSlice(double *halo, const double *val, const long int *nGPointsProd, const int *nGPoints, int nDims, int d, int offset){
	val += offset*nGPointsProd[d];
	getSliceInner(halo,&val,&nGPointsProd[nDims-1],&nGPoints[nDims-1],nGPointsProd[d]);
}

inline void getSlice2(double *halo, const double *val, const long int *nGPointsProd, const int *nGPoints, int nDims, int d, int offset){
	val += offset*nGPointsProd[d];
	getSliceInner(halo,&val,&nGPointsProd[nDims-1],&nGPoints[nDims-1],nGPointsProd[d]);
}

/*
 * MAIN ROUTINE
 */

int main(int argc, char **argv){

	/*
	 * MAKING ARTIFICIAL DATA
	 */

	int nDims = 4;

	int *nGPoints = malloc(nDims*sizeof(int));
	int *nTGPoints = malloc(nDims*sizeof(int));
	long int *nGPointsProd = malloc((nDims+1)*sizeof(long int));

	nGPoints[0] = 3;
	if(nDims>1) nGPoints[1] = 3;
	if(nDims>2) nGPoints[2] = 3;
	if(nDims>3) nGPoints[3] = 3;

	nGPoints[0] = 256;
	nGPoints[1] = nGPoints[0];
	nGPoints[2] = nGPoints[0];

	nGPointsProd[0] = 1;
	for(int d=1;d<=nDims;d++){
		nGPointsProd[d] = nGPointsProd[d-1] * nGPoints[d-1];
	}

	for(int d=0;d<nDims;d++){
		nTGPoints[d] = nGPoints[d]-2;
	}

	// Interior and exterior (boundary) number of elements
	long int nInt = 1;
	for(int d=0;d<nDims;d++){
		nInt *= nTGPoints[d];
	}
	long int nExt = nGPointsProd[nDims]-nInt;

	// Number of elements on edge at dimension d
	long int *nSlice = malloc(nDims*sizeof(nSlice));
	for(int d=0;d<nDims;d++){
		nSlice[d] = 1;
		for(int dd=0;dd<nDims;dd++){
			if(dd!=d) nSlice[d] *= nGPoints[dd];
		}
	}
	long int nSliceMax = 0;
	for(int d=0;d<nDims;d++){
		if(nSlice[d]>nSliceMax) nSliceMax = nSlice[d];
	}

	double *val = malloc(nGPointsProd[nDims]*sizeof(double));
	for(long int p=0;p<nGPointsProd[nDims];p++) val[p] = p;

	double *slice = malloc(nSliceMax*sizeof(double));
	double *slice2 = malloc(nSliceMax*sizeof(double));


	/*
	 * TESTING ALGORITHMS
	 */

	 int d=1;
	 int offset=2;
	 int N=1000;

	tprobe(NULL);
	long int *indices = malloc(nSliceMax*sizeof(long int));
	getSliceIndices(indices,nGPoints,nGPointsProd,nDims,d,offset);
	tprobe("Get indices once (naively)");

	tprobe(NULL);
	for(int i=0;i<N;i++) getFromInd(slice,val,indices,nSlice[d]);
	tprobe("Get slice using pre-stored indices 1000 times");

	tprobe(NULL);
	for(int i=0;i<N;i++) getSlice(slice2,val,nGPointsProd,nGPoints,nDims,d,offset);
	tprobe("Get slice directly using recursive scheme 1000 times");

//	tprobe(NULL);
//	for(int i=0;i<N;i++) getSlice2(slice2,val,nGPointsProd,nGPoints,nDims,d,offset);
//	tprobe("Get slice directly using recursive scheme 1000 times");

//	tprobe(NULL);
//	for(int i=0;i<N;i++) getSlice2(slice2,val,nGPointsProd,nGPoints,nDims,d,offset);
//	tprobe("Get slice directly using recursive scheme 2 1000 times");

//	factorial(20,0);

	// Compare that they are equal
	int success = 1;
	for(int p=0;p<nSlice[d];p++){
		if(slice[p]!=slice2[p]){
			printf("ERROR: slice[%i]=%f!=%f=slice2[%i]\n",p,slice[p],slice2[p],p);
			success = 0;
		}
	}

	if(success) printf("SUCCESS\n");

}
