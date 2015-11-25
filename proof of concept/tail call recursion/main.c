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

void printSP2() {
  void* p = NULL;
  printf("Stack pointer: %p\n", (void*)&p);
}

void printSP(){
	register const long rsp __asm__ ("rsp");
	printf("Stack pointer: %p\n", (void*)rsp);
}

long int factorial(long int n){

	if(n==5){
		volatile int c = n+3;
//		printf("here\n");
	}

//	printSP();
//	printSP();
	void * addr = __builtin_frame_address(0);
	printf("Frame address: %p\n",addr);

	if(n==0) return 1;
	else return n*factorial(n-1);
}

int main(int argc, char **argv){

	int n=10;
	long int l = factorial(n);
	printf("%i!=%li\n",n,l);
}
