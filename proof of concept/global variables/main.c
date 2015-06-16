#include <stdlib.h>
//#include <stdio.h>
#include <math.h>
#include "../stool/stool.h"

#define ARRSIZE 1000000

double arrGlobal[ARRSIZE];
double *arrGDyn;


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

int main(){

	double arrMain[ARRSIZE];
	double *arrDyn;
	arrDyn = (double *)malloc(ARRSIZE*sizeof(double));
	arrGDyn = (double *)malloc(ARRSIZE*sizeof(double));


	tprobe("Start");

	setArrGlobal(10);
	for(int i=0;i<100;i++) incArrGlobal();

	tprobe("Working on global array (heap?)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<100;i++) incArrLocal(arrMain,ARRSIZE);

	tprobe("Working on local array (stack)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<100;i++) incArrLocal(arrDyn,ARRSIZE);

	tprobe("Working on local dynamic array (heap)");

	setArrGlobal(10);
	for(int i=0;i<100;i++) incArrGlobalUnordered();

	tprobe("Unordered working on global array (heap?)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<100;i++) incArrLocalUnordered(arrMain,ARRSIZE);

	tprobe("Unordered working on local array (stack)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<100;i++) incArrLocalUnordered(arrDyn,ARRSIZE);

	tprobe("Unordered working on local dynamic array (heap)");

	setArrLocal(arrMain,ARRSIZE,10);
	for(int i=0;i<100;i++) incArrLocalUnordered(arrGDyn,ARRSIZE);

	tprobe("Unordered working on global dynamic array (heap)");

	vprint(arrGlobal,15);
	vprint(arrMain,15);
	
	free(arrDyn);


	return EXIT_SUCCESS;
}
