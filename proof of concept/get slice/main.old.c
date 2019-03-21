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

/*
 * METHOD WITH PRE-STORED INDICES
 */

long int *getGhostIndices(long int *indices, const int *nGPoints, const long int *nGPointsProd, int nDims){
 	for(long int p=0; p<nGPointsProd[nDims]; p++){

 		long int temp = p;
 		for(int d=0; d<nDims; d++){

 			// Find index d and prepare for next iteration
 			int index = temp % nGPoints[d];
 			temp /= nGPoints[d];

 			if(index==0 || index==nGPoints[d]-1){
 				*(indices++) = p;
 				break;	// avoids double entries of corners etc.
 			}
 		}
 	}
}

int getIndHalo(double *halo, double *val, long int *indices, int nExt){
	for(long int p=0;p<nExt;p++){
		halo[p] = val[indices[p]];
	}
}



/*
 * DIRECT METHOD
 */

double *getGhostsOld(double *nextGhost, const double **valp, const long int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++(*valp));
	if(*mul!=1) for(int k=0;k<*points-2;k++) nextGhost = getGhostsOld(nextGhost,valp,mul-1,points-1);
	else *valp += *points-2;
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++(*valp));
	return nextGhost;
}

void getHaloOld(double *halo, const double *val, const long  int *nGPointsProd, const int *nGPoints,int nDims){
	const double **valp = malloc(sizeof(double*));
	*valp = val-1;
	getGhostsOld(halo,valp,&nGPointsProd[nDims-1],&nGPoints[nDims-1]);
	free(valp);
}

double *getSlice(double *nextGhost, const double **valp, const long int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	if(*mul!=1) for(int k=0;k<*points;k++) nextGhost = getSlice(nextGhost,valp,mul-1,points-1);
	else *valp += *points;
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	return nextGhost;
}

inline void getHalo(double *halo, const double *val, const long  int *nGPointsProd, const int *nTGPoints,int nDims){
	getSlice(halo,&val,&nGPointsProd[nDims-1],&nTGPoints[nDims-1]);
}

/*
 * DIRECT METHOD FOR ONE BOUNDARY EDGE ONLY
 */

double *getGridSlice(double *nextGhost, const double **valp, const long int *mul, const int *points, const long int finalMul){
/*
	for(int j=0;j<mul[-1];j++) *(nextGhost++) = *((*valp)++);
	*valp += (*points+1)*(mul[-1]);
	if(*mul==finalMul) return nextGhost;
	else nextGhost = getSliceEdge(nextGhost,valp,mul-1,points-1,finalMul);
*/

	printf("Call to getGridSlice: *mul=%li, *points=%i\n",*mul,*points);

	if(*mul==finalMul){
		for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
		*valp += (*mul)*(*points-1);
		printf("nextGhost=%f\n",*(nextGhost-1));
	} else {
		for(int j=0;j<*points;j++) nextGhost = getGridSlice(nextGhost,valp,mul-1,points-1,finalMul);
	}
	return nextGhost;
}

inline void getHaloSlice(double *halo, const double *val, const long int *nGPointsProd, const int *nGPoints, int nDims, int d, int offset){
	val += offset*nGPointsProd[d];
	getGridSlice(halo,&val,&nGPointsProd[nDims-1],&nGPoints[nDims-1],nGPointsProd[d]);
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

	nGPoints[0] = 5;
	if(nDims>1) nGPoints[1] = 4;
	if(nDims>2) nGPoints[2] = 3;
	if(nDims>3) nGPoints[3] = 3;

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
	long int *nEdge = malloc(nDims*sizeof(nEdge));
	for(int d=0;d<nDims;d++){
		nEdge[d] = 1;
		for(int dd=0;dd<nDims;dd++){
			if(dd!=d) nEdge[d] *= nGPoints[dd];
		}
	}
	long int nEdgeMax = 0;
	for(int d=0;d<nDims;d++){
		if(nEdge[d]>nEdgeMax) nEdgeMax = nEdge[d];
	}

	double *val = malloc(nGPointsProd[nDims]*sizeof(double));
	for(long int p=0;p<nGPointsProd[nDims];p++) val[p] = p;

	double *halo = malloc(nExt*sizeof(double));
	double *halo2 = malloc(nExt*sizeof(double));
	double *halo3 = malloc(nExt*sizeof(double));

	double *haloSlice = malloc(nEdgeMax*sizeof(double));


	/*
	 * TESTING ALGORITHMS
	 */

	 int d=1;
	 getHaloSlice(haloSlice,val,nGPointsProd,nGPoints,nDims,d,2);

	 for(long int p=0;p<nEdge[d];p++)
	 	printf("haloSlice[%li] = %f\n",p,haloSlice[p]);

/*
	long int N = 1000;

	tprobe(NULL);
	long int *indices = malloc(nExt*sizeof(long int));
	getGhostIndices(indices,nGPoints,nGPointsProd,nDims);
	tprobe("Get indices once");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getIndHalo(halo2,val,indices,nExt);
	tprobe("Get halo using pre-stored indices 1000 times");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getHalo(halo,val,nGPointsProd,nTGPoints,nDims);
	tprobe("Get halo directly using recursive scheme 1000 times");

	// Compare that they are equal
	int success = 1;
	for(int p=0;p<nExt;p++){
		if(halo[p]!=halo2[p]){
			printf("ERROR: halo[%i]=%f!=%f=halo2[%i]\n",p,halo[p],halo2[p],p);
			success = 0;
		}
	}

	if(success) printf("SUCCESS\n");
*/
}
