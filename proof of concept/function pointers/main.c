/*
 * This document is to test whether it is practical to select
 * which function (or "module") to use certain places in
 * the code simply by executing a pointer to a function that
 * can be set to point to the appropriate function during
 * initialization. Functions really just points to a certain
 * area in memory space so maybe it's not even slower?
 */

#include <stdlib.h>
//#include <stdio.h>
#include <math.h>
#include "../stool/stool.h"

void addTwo(double *arr, long int arrsize){

	for(long int i=0;i<arrsize;i++){
		arr[i] += 2;
	}

}

void addThree(double *arr, long int arrsize){

	for(long int i=0;i<arrsize;i++){
		arr[i] += 3;
	}

}

int main(int argc, char **argv){

	int mode=0;	// 0 - use addTwo, 1 - use addThree
	int N=10;	// Elements in array
	long int M=100000000; // Iterations

	if(argv[0][0]=='a') mode=1;

	double* arr;
	arr = (double *)malloc(N*sizeof(double));
	vsetall(arr,N,0);

	// Choose function before loop using function pointer
	void (*myFunc)(double*,long int);
	if(mode) myFunc=&addThree;
	else myFunc=&addTwo;

	tprobe("Start");

	for(long int i=0;i<M;i++){
		addTwo(arr,N);
	}

	tprobe("No choice");

	for(long int i=0;i<M;i++){
		if(mode) addThree(arr,N);
		else addTwo(arr,N);
	}
	mode = 1;

	tprobe("Choice by if sentence");

	for(long int i=0;i<M;i++){
		(*myFunc)(arr,N);
	}

	tprobe("Choice by function pointer");

	vprint(arr,5);
	

	

	free(arr);
	return EXIT_SUCCESS;

}
