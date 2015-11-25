#define _XOPEN_SOURCE 700

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "halo.h"

/*
 * TIMING FUNCTION (Ignore this)
 */
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
 * YOUR METHOD
 */



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

	// Interior and boundary parts of domain
	long int nInt = (nGPoints[0]-2)*(nGPoints[1]-2)*(nGPoints[2]-2);
	long int nExt = nGPointsProd[nDims]-nInt;

	double *halo = malloc(nExt*sizeof(double));
	double *halo2 = malloc(nExt*sizeof(double));

	/*
	 * TESTING ALGORITHMS
	 */


	long int N = 1000;

	tprobe(NULL);
	// INSERT YOUR METHOD HERE
	tprobe("Your method");

	tprobe(NULL);
	for(long int i=0;i<N;i++) getHalo(halo2,val,nGPointsProd,nTGPoints,nDims);
	tprobe("Sigvald's method");

	/*
	 * CHECKING CORRECTNESS
	 */

	int success = 1;
	for(int p=0;p<nExt;p++){
		if(halo[p]!=halo2[p]){
			printf("ERROR: halo[%i]=%f!=%f=halo2[%i]\n",p,halo[p],halo2[p],p);
			success = 0;
		}
	}

	if(success) printf("SUCCESS\n");

}
