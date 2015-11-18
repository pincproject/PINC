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

double *getGhosts(double *nextGhost, const double **valp, const long int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++(*valp));
	if(*mul!=1) for(int k=0;k<*points-2;k++) nextGhost = getGhosts(nextGhost,valp,mul-1,points-1);
	else *valp += *points-2;
	for(int j=0;j<*mul;j++) *(nextGhost++) = *(++(*valp));
	return nextGhost;
}

inline void getHalo(double *halo, const double *val, const long  int *nGPointsProd, const int *nGPoints,int nDims){
	const double **valp = malloc(sizeof(double*));
	*valp = val-1;
	getGhosts(halo,valp,&nGPointsProd[nDims-1],&nGPoints[nDims-1]);
	free(valp);
}

double *getGhosts2(double *nextGhost, const double **valp, const long int *mul, const int *points){
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	if(*mul!=1) for(int k=0;k<*points;k++) nextGhost = getGhosts2(nextGhost,valp,mul-1,points-1);
	else *valp += *points;
	for(int j=0;j<*mul;j++) *(nextGhost++) = *((*valp)++);
	return nextGhost;
}

inline void getHalo2(double *halo, const double *val, const long  int *nGPointsProd, const int *nTGPoints,int nDims){
	const double **valp = malloc(sizeof(double*));
	*valp = val;
	getGhosts2(halo,valp,&nGPointsProd[nDims-1],&nTGPoints[nDims-1]);
	free(valp);
}

/*
 * MAIN ROUTINE
 */

int main(int argc, char **argv){

	/*
	 * MAKING ARTIFICIAL DATA
	 */

	int nDims = 3;

	int *nGPoints = malloc(nDims*sizeof(int));
	int *nTGPoints = malloc(nDims*sizeof(int));
	long int *nGPointsProd = malloc((nDims+1)*sizeof(long int));

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

	double *val = malloc(nGPointsProd[nDims]*sizeof(double));
	for(long int p=0;p<nGPointsProd[nDims];p++) val[p] = p;

	// Interior and exterior parts of domain
	long int nInt = (nGPoints[0]-2)*(nGPoints[1]-2)*(nGPoints[2]-2);
	long int nExt = nGPointsProd[nDims]-nInt;

	double *halo = malloc(nExt*sizeof(double));
	double *halo2 = malloc(nExt*sizeof(double));

	/*
	 * TESTING ALGORITHMS
	 */


	long int N = 1000;

	tprobe(NULL);
	long int *indices = malloc(nExt*sizeof(long int));
	getGhostIndices(indices,nGPoints,nGPointsProd,nDims);
	tprobe("Get indices once");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getIndHalo(halo,val,indices,nExt);
	tprobe("Get halo using pre-stored indices 1000 times");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getHalo(halo2,val,nGPointsProd,nGPoints,nDims);
	tprobe("Get halo directly using recursive scheme 1000 times");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getHalo2(halo,val,nGPointsProd,nTGPoints,nDims);
	tprobe("Get halo2 directly using recursive scheme 1000 times");

	tprobe(NULL);
//	for(long int i=0;i<N;i++) getGhostIndices(indices,nGPoints,nGPointsProd,nDims);
	for(long int i=0;i<N;i++) getIndHalo(halo,val,indices,nExt);
	tprobe("Get indices and halo 1000 times");

	// Compare that they are equal
	int success = 1;
	for(int p=0;p<nExt;p++){
		if(halo[p]!=halo2[p]){
			printf("ERROR: halo[%i]=%f!=%f=halo2[%i]\n",p,halo[p],halo2[p],p);
			success = 0;
		}
	}
	if(success) printf("SUCCESS\n");

}
