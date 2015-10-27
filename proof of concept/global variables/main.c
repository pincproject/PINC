#include <stdlib.h>
//#include <stdio.h>
#include <math.h>
#include "../stool/stool.h"

#define ARRSIZE 1000000

// Script to test performance of global variables which in turn renders it
// unnecessary to pass arguments.

struct ArrStruct{
	double *arr;
};

double arrGlobal[ARRSIZE];
double *arrGDyn;
struct ArrStruct myGStruct;

void incArrLocal(double *arr, long int arrsize){

	for(long int i=0;i<arrsize;i++){
		arr[i]++;
	}

}

void incArrGlobal(){

	for(long int i=0;i<ARRSIZE;i++){
		arrGlobal[i]++;
	}

}

void incArrGDyn(){

	for(long int i=0;i<ARRSIZE;i++){
		arrGDyn[i]++;
	}

}

void incArrLocalUnordered(double *arr, long int arrsize){

	for(int j=0;j<10;j++){
		for(long int i=0;i+j<arrsize;i=i+10){
			arr[i+j]++;
		}
	}

}

void incArrGlobalUnordered(){

	for(int j=0;j<10;j++){
		for(long int i=0;i+j<ARRSIZE;i=i+10){
			arrGlobal[i+j]++;
		}
	}

}

void incArrGDynUnordered(){

	for(int j=0;j<10;j++){
		for(long int i=0;i+j<ARRSIZE;i=i+10){
			arrGDyn[i+j]++;
		}
	}

}

void setArrLocal(double *arr, long int arrsize, int value){

	for(long int i=0;i<arrsize;i++){
		arr[i]=value;
	}

}

void setArrGlobal(int value){

	for(long int i=0;i<ARRSIZE;i++){
		arrGlobal[i]=value;
	}

}

void setArrGDyn(int value){

	for(long int i=0;i<ARRSIZE;i++){
		arrGDyn[i]=value;
	}

}

// Local struct housing array
// Not really a fair comparison since the global arrays in use today are of
// struct type.
void setArrStruct(struct ArrStruct *myStruct, long int arrsize, int value){

	for(long int i=0;i<arrsize;i++){
		myStruct->arr[i]=value;
	}

}

void incArrStruct(struct ArrStruct *myStruct, long int arrsize){

	for(long int i=0;i<arrsize;i++){
		myStruct->arr[i]++;
	}
}

void incArrStructUnordered(struct ArrStruct *myStruct, long int arrsize){

	for(int j=0;j<10;j++){
		for(long int i=0;i+j<arrsize;i=i+10){
			myStruct->arr[i+j]++;
		}
	}

}

void setGArrStruct(int value){

	for(long int i=0;i<ARRSIZE;i++){
		myGStruct.arr[i]=value;
	}

}

void incGArrStruct(){

	for(long int i=0;i<ARRSIZE;i++){
		myGStruct.arr[i]++;
	}
}

void incGArrStructUnordered(){

	for(int j=0;j<10;j++){
		for(long int i=0;i+j<ARRSIZE;i=i+10){
			myGStruct.arr[i+j]++;
		}
	}

}


int main(){

	double arrMain[ARRSIZE];
	double *arrDyn;
	struct ArrStruct myStruct;
	myStruct.arr = (double *)malloc(ARRSIZE*sizeof(double));
	myGStruct.arr = (double *)malloc(ARRSIZE*sizeof(double));
	arrDyn = (double *)malloc(ARRSIZE*sizeof(double));
	arrGDyn = (double *)malloc(ARRSIZE*sizeof(double));
	int N=1000;


	tprobe("Start");

	/*
	 * ORDERED WORKING
	 */

	setArrGlobal(10);
	for(int i=0;i<N;i++) incArrGlobal();

	tprobe("Working on global static array (heap?)");

	setArrGDyn(10);
	for(int i=0;i<N;i++) incArrGDyn();

	tprobe("Working on global dynamic array (heap?)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrLocal(arrMain,ARRSIZE);

	tprobe("Working on local static array (stack)");

	setArrLocal(arrDyn,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrLocal(arrDyn,ARRSIZE);

	tprobe("Working on local dynamic array (heap)");

	setGArrStruct(10);
	for(int i=0;i<N;i++) incGArrStruct();

	tprobe("Working on global struct with dynamic array (heap?)");

	setArrStruct(&myStruct,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrStruct(&myStruct,ARRSIZE);

	tprobe("Working on local struct with dynamic array (heap)");

	vprint(arrGlobal,5);
	vprint(arrGDyn,5);
	vprint(arrMain,5);
	vprint(arrDyn,5);
	vprint(myGStruct.arr,5);
	vprint(myStruct.arr,5);
	tprobe("");

	/*
	 * UNORDERED WORKING (TO GENERATE CACHE MISSES)
	 */

	setArrGlobal(10);
	for(int i=0;i<N;i++) incArrGlobalUnordered();

	tprobe("Unordered working on global static array (heap?)");

	setArrGDyn(10);
	for(int i=0;i<N;i++) incArrGDynUnordered();

	tprobe("Unordered working on global dynamic array (heap)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrLocalUnordered(arrMain,ARRSIZE);

	tprobe("Unordered working on local static array (stack)");

	setArrLocal(arrDyn,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrLocalUnordered(arrDyn,ARRSIZE);

	tprobe("Unordered working on local dynamic array (heap)");

	setGArrStruct(10);
	for(int i=0;i<N;i++) incGArrStructUnordered();

	tprobe("Unordered working on global struct with dynamic array (heap?)");

	setArrStruct(&myStruct,ARRSIZE,10);
	for(int i=0;i<N;i++) incArrStructUnordered(&myStruct,ARRSIZE);

	tprobe("Unordered working on local struct with dynamic array (heap)");

	vprint(arrGlobal,5);
	vprint(arrGDyn,5);
	vprint(arrMain,5);
	vprint(arrDyn,5);
	vprint(myGStruct.arr,5);
	vprint(myStruct.arr,5);
	
	free(arrDyn);
	free(arrGDyn);
	free(myStruct.arr);
	free(myGStruct.arr);


	return EXIT_SUCCESS;
}
