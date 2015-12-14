#define _XOPEN_SOURCE 700
#define _GNU_SOURCE	// Gives us ffs() to "find first bit set"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>

//#define inline

/*
 * VALUABLE SOURCES:
 * http://stackoverflow.com/questions/757059/position-of-least-significant-bit-that-is-set
 * http://www.jjj.de/fxt/fxtbook.pdf (section 1.16.3 in particular)
 * http://supertech.csail.mit.edu/papers/debruijn.pdf
 */

/*
 * TIMING FUNCTION
 */

// Time probe. See example above.
struct timespec tprobe(const char * label){

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

	return diff;

}

void printSP(){
	register const long rsp __asm__ ("rsp");
	printf("Stack pointer: %p\n", (void*)rsp);
}

void printBinary(int num, int nBits){

	int *bits = malloc(nBits*sizeof(*bits));

	for(int i=0;i<nBits;i++){
		int bit = num%2;
		bits[i] = bit;
		num /= 2;
	}

	for(int i=0;i<nBits;i++)
		printf("%i",bits[nBits-i-1]);

	free(bits);

}

/*
 * GRAY METHOD
 */

#define DEBRUIJN 0x077CB531U
static const int debruijn[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

inline int incGray(int x){
	x ^= 1;
	x ^= (x&-x)<<1;
	return x;
}

inline int toGray(int x){
	int gray = x^(x>>1);
	return gray;
}

inline int getIncGrayMask(int *x){

	int next = incGray(*x);
	int mask = (next^*x)>>1;
	*x = next;
	return mask;
}

inline int getIncGrayMask2(int *x, int i){
	int next = toGray(i);
	int mask = (next^*x);
	*x = next;
	return mask;
}

inline int getBitNo(int mask){

	int ind = 0;
	while(mask>1){
		mask >>= 1;
		ind++;
	}
	return ind;
}

static inline int getBitNoDebruijn(int mask){

	mask &= -mask;
	mask *= DEBRUIJN;
	mask >>= 27;
	return debruijn[mask];

}

inline double getWeight(int signMask, const double *decimal, const double *complement, int nDims){

	double weight = 1;
	for(int d=0;d<nDims;d++){
		int sign = (signMask&(1<<d)) && 1; // 1 represents +, 0 repr. -
		weight *= sign*decimal[d] - sign*complement[d] + complement[d];
	}
	return weight;

}

static inline double interpGray(const double *val, const double *pos, const long int *nGPointsProd, int nDims, double *decimal, double *complement, int nCorners){

	long int p = 0;
	double weight = 1;
	for(int d=0;d<nDims;d++){
		int integer = (int)pos[d];
		decimal[d] = pos[d]-integer;
		complement[d] = 1-decimal[d];
		p += nGPointsProd[d] * integer;
		weight *= complement[d];
	}
	double result = weight*val[p];

	int gray = 0;
	for(int c=1;c<nCorners;c++){

		int prev = gray;
		gray = toGray(c);
		int mask = (prev^gray);
		int d = getBitNoDebruijn(mask);

		int sign = (gray&mask) && 1; // 1 represents +, 0 repr. -
		p += (2*sign-1) *nGPointsProd[d];

		weight = getWeight(gray,decimal,complement,nDims);
		result += weight*val[p];
	}

	return result;
}


static inline double interpCGray(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int nCorners){

	long int p = 0;
	for(int d=0;d<nDims;d++){
		p += nGPointsProd[d] * (int)pos[d];
	}
	double result = val[p];

	int gray = 0;
	for(int c=1;c<nCorners;c++){

		int prev = gray;
		gray = toGray(c);
		int mask = (prev^gray);
		int d = getBitNoDebruijn(mask);

		int sign = (gray&mask) && 1; // 1 represents +, 0 repr. -
		p += (2*sign-1) *nGPointsProd[d];

		result += val[p];
	}

	return 0.125*result;
}

double interpAllGray(const double *val, const double *pos, const long int *nGPointsProd, int nDims, double *decimal, double *complement, int nCorners, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpGray(val,&pos[nDims*i],nGPointsProd,nDims,decimal,complement,nCorners);
	}
	return result;

}



double interpCAllGray(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int nCorners, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpCGray(val,&pos[nDims*i],nGPointsProd,nDims,nCorners);
	}
	return result;

}


/*
 * RECURSIVE METHOD
 */
void inner(double *result, const double *val, const long int *mul, long int p, double factor, double *decimal, double *complement){

	//printSP();

	if(*mul==1){
		*result += *complement*factor*val[p];
		*result +=    *decimal*factor*val[p+1];
	} else {
		inner(result,val,mul-1,p,  *complement*factor,decimal-1,complement-1);
		inner(result,val,mul-1,p+*mul,*decimal*factor,decimal-1,complement-1);
	}

}

double inner2(const double *val, const long int *mul, long int p, double *decimal, double *complement){

	//printSP();

	double result;
	if(*mul==1){
		result  = *complement*val[p];		// stay
		result += *decimal   *val[p+1];		// incr.
	} else {
		result  = *complement*inner2(val,mul-1,p     ,decimal-1,complement-1);	// stay
		result += *decimal   *inner2(val,mul-1,p+*mul,decimal-1,complement-1);	// incr.
	}
	return result;

}

double innerC2(const double *val, const long int *mul, long int p){

	//printSP();

	double result;
	if(*mul==1){
		result  = val[p];		// stay
		result += val[p+1];		// incr.
	} else {
		result  = innerC2(val,mul-1,p     );	// stay
		result += innerC2(val,mul-1,p+*mul);	// incr.
	}
	return result;

}

void innerC(double *result, const double *val, const long int *mul, long int p){

	// printSP();

	if(*mul==1){
		*result += val[p];
		*result += val[p+1];
	} else {
		innerC(result,val,mul-1,p);
		innerC(result,val,mul-1,p+*mul);
	}

}

inline double interpRecursive(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int *integer, double *decimal, double *complement){

	double result = 0;
	long int p = 0;

	for(int d=0;d<nDims;d++){
		integer[d] = (int)pos[d];
		decimal[d] = pos[d]-integer[d];
		complement[d] = 1-decimal[d];
		p += nGPointsProd[d] * integer[d];
	}

	inner(&result,val,&nGPointsProd[nDims-1],p,1,&decimal[nDims-1],&complement[nDims-1]);
	// innerSimpl(&result,val,&nGPointsProd[nDims-1],p);

	return result;

}

inline double interpRecursive2(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int *integer, double *decimal, double *complement){


	long int p = 0;
	for(int d=0;d<nDims;d++){
		integer[d] = (int)pos[d];
		decimal[d] = pos[d]-integer[d];
		complement[d] = 1-decimal[d];
		p += nGPointsProd[d] * integer[d];
	}

	return inner2(val,&nGPointsProd[nDims-1],p,&decimal[nDims-1],&complement[nDims-1]);

}

inline double interpCRecursive2(const double *val, const double *pos, const long int *nGPointsProd, int nDims){


	long int p = 0;
	for(int d=0;d<nDims;d++){
		p += nGPointsProd[d] * (int)pos[d];
	}

	return 0.125*innerC2(val,&nGPointsProd[nDims-1],p);

}

inline double interpCRecursive(const double *val, const double *pos, const long int *nGPointsProd, int nDims){

	double result = 0;
	long int p = 0;

	for(int d=0;d<nDims;d++){
		p += nGPointsProd[d] * (int)pos[d];
	}

	innerC(&result,val,&nGPointsProd[nDims-1],p);

	return 0.125*result;

}


double interpAllRecursive(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int *integer, double *decimal, double *complement, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpRecursive(val,&pos[nDims*i],nGPointsProd,nDims,integer,decimal,complement);
	}
	return result;

}

double interpAllRecursive2(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int *integer, double *decimal, double *complement, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpRecursive2(val,&pos[nDims*i],nGPointsProd,nDims,integer,decimal,complement);
	}
	return result;

}

double interpCAllRecursive(const double *val, const double *pos, const long int *nGPointsProd, int nDims, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpCRecursive(val,&pos[nDims*i],nGPointsProd,nDims);
	}
	return result;

}

double interpCAllRecursive2(const double *val, const double *pos, const long int *nGPointsProd, int nDims, long int nParticles){

	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpCRecursive2(val,&pos[nDims*i],nGPointsProd,nDims);
	}
	return result;

}

/*
 * PURE BINARY INDEXING METHOD
 */

inline double getWeightAndIndex(long int *p, const int mask, const double *decimal, const double *complement, int nDims, const long int *nGPointsProd){

	double weight = 1;
	for(int d=0;d<nDims;d++){
		int sign = (mask&(1<<d)) && 1; // 1 represents +, 0 repr. -
		weight *= sign*decimal[d] - sign*complement[d] + complement[d];
		*p += sign*nGPointsProd[d];
	}
	return weight;
}

inline double interpBinary(const double *val, const double *pos, const long int *nGPointsProd, int nDims, double *decimal, double *complement, int nCorners){

	long int p = 0;
	// double factor = 1;
	for(int d=0;d<nDims;d++){
		int integer = (int)pos[d];
		decimal[d] = pos[d]-integer;
		complement[d] = 1-decimal[d];
		p += nGPointsProd[d] * integer;
	}
	double result = 0;

	int signMask;
	for(int i=0;i<nCorners;i++){

		long int pp = p;
		double weight = getWeightAndIndex(&pp,i,decimal,complement,nDims,nGPointsProd);
		result += weight*val[pp];
	}

	return result;

}

inline void getCIndex(long int *p, const int mask, int nDims, const long int *nGPointsProd){

	for(int d=0;d<nDims;d++){
		int sign = (mask&(1<<d)) && 1; // 1 represents +, 0 repr. -
		*p += sign*nGPointsProd[d];
	}
}

inline double interpCBinary(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int nCorners){

	long int p = 0;
	// double factor = 1;
	for(int d=0;d<nDims;d++){
		p += nGPointsProd[d] * (int)pos[d];
	}
	double result = 0;

	int signMask;
	for(int i=0;i<nCorners;i++){

		long int pp = p;
		getCIndex(&pp,i,nDims,nGPointsProd);
		result += val[pp];
	}

	return 0.125*result;

}

double interpAllBinary(const double *val, const double *pos, const long int *nGPointsProd, int nDims, double *decimal, double *complement, int nCorners, long int nParticles){
	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpBinary(val,&pos[nDims*i],nGPointsProd,nDims,decimal,complement,nCorners);
	}
	return result;
}

double interpCAllBinary(const double *val, const double *pos, const long int *nGPointsProd, int nDims, int nCorners, long int nParticles){
	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpCBinary(val,&pos[nDims*i],nGPointsProd,nDims,nCorners);
	}
	return result;
}

/*
 * FIXED 3D METHOD
 */

inline double interp3D(const double *val, const double *pos, const long int *nGPointsProd){

	// Integer parts of position
	int j = (int) pos[0];
	int k = (int) pos[1];
	int l = (int) pos[2];

	// Decimal parts of position and their complement
	double x = pos[0]-j;
	double y = pos[1]-k;
	double z = pos[2]-l;
	double xcomp = 1-x;
	double ycomp = 1-y;
	double zcomp = 1-z;

	// Index of neighbouring nodes
	long int p 		= j + k*nGPointsProd[1] + l*nGPointsProd[2];
	long int pj 	= p + 1;
	long int pk 	= p + nGPointsProd[1];
	long int pjk 	= pk + 1;
	long int pl 	= p + nGPointsProd[2];
	long int pjl 	= pl + 1;
	long int pkl 	= pl + nGPointsProd[1];
	long int pjkl 	= pkl + 1;

	// Linear interpolation
	double result =	 zcomp*( ycomp*(xcomp*val[p   ]+x*val[pj  ])
							+y    *(xcomp*val[pk  ]+x*val[pjk ]) )
					+z    *( ycomp*(xcomp*val[pl  ]+x*val[pjl ])
							+y    *(xcomp*val[pkl ]+x*val[pjkl]) );

	return result;

}

 inline double interpC3D(const double *val, const double *pos, const long int *nGPointsProd){

	// Integer parts of position
	int j = (int) pos[0];
	int k = (int) pos[1];
	int l = (int) pos[2];

	// Index of neighbouring nodes
	long int p 		= j + k*nGPointsProd[1] + l*nGPointsProd[2];
	long int pj 	= p + 1;
	long int pk 	= p+nGPointsProd[1];
	long int pjk 	= pk + 1;
	long int pl 	= p + nGPointsProd[2];
	long int pjl 	= pl + 1;
	long int pkl 	= pl + nGPointsProd[1];
	long int pjkl 	= pkl + 1;

	// Linear interpolation
	double result =	 0.125*(val[p]+val[pj]+val[pk]+val[pjk]+val[pl]+val[pjl]+val[pkl]+val[pjkl]);
	return result;

}

double interpAll3D(const double *val, const double *pos, const long int *nGPointsProd, int nDims, long int nParticles){
	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interp3D(val,&pos[nDims*i],nGPointsProd);
	}
	return result;
}

double interpCAll3D(const double *val, const double *pos, const long int *nGPointsProd, int nDims, long int nParticles){
	double result = 0;
	for(long int i=0;i<nParticles;i++){
		result += interpC3D(val,&pos[nDims*i],nGPointsProd);
	}
	return result;
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
	long int *nGPointsProd = malloc((nDims+1)*sizeof(long int));

	nGPoints[0] = 5;
	if(nDims>1) nGPoints[1] = 4;
	if(nDims>2) nGPoints[2] = 3;
	if(nDims>3) nGPoints[3] = 3;

	nGPoints[0] = 128;
	nGPoints[1] = 128;
	nGPoints[2] = 128;
	// nGPoints[3] = 129;
	// nGPoints[4] = 129;

	nGPointsProd[0] = 1;
	for(int d=1;d<=nDims;d++){
		nGPointsProd[d] = nGPointsProd[d-1] * nGPoints[d-1];
	}

	double *val = malloc(nGPointsProd[nDims]*sizeof(*val));
	for(long int p=0;p<nGPointsProd[nDims];p++) val[p] = p;

	int *integer = malloc(nDims*sizeof(*integer));
	double *decimal = malloc(nDims*sizeof(*decimal));
	double *complement = malloc(nDims*sizeof(*complement));
	double *weight = malloc((nDims+1)*sizeof(*weight));
	weight[nDims]=1;

	/*
	 * TESTING ALGORITHMS
	 */


	long int nParticles = 1e7;
	nParticles = strtol(argv[1],NULL,10);
	printf("nParticles=%li\n",nParticles);
	double *pos = malloc(nDims*nParticles*sizeof(*pos));
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	for(long int i=0;i<nParticles*nDims;i++){
		pos[i] = gsl_rng_uniform(rng)*(nGPoints[i%nDims]-1);
	}
	gsl_rng_free(rng);
	// pos[0] = 2.5;
	// pos[1] = 0.5;
	// pos[2] = 0.5;

	printf("True interpolation:\n");
	int nCorners = pow(2,nDims);

	double result0 = interpAll3D(val,pos,nGPointsProd,nDims,nParticles);

	tprobe(NULL);
	double result1 = interpAll3D(val,pos,nGPointsProd,nDims,nParticles);
	struct timespec t1 = tprobe("Fixed 3D method");

	tprobe(NULL);
	double result2 = interpAllRecursive(val,pos,nGPointsProd,nDims,integer,decimal,complement,nParticles);
	struct timespec t2 = tprobe("Recursive method");

	tprobe(NULL);
	double result3 = interpAllBinary(val,pos,nGPointsProd,nDims,decimal,complement,nCorners,nParticles);
	struct timespec t3 = tprobe("Binary iterative method");

	tprobe(NULL);
	double result4 = interpAllGray(val,pos,nGPointsProd,nDims,decimal,complement,nCorners,nParticles);
	struct timespec t4 = tprobe("Gray iterative method");

	tprobe(NULL);
	double result5 = interpAllRecursive2(val,pos,nGPointsProd,nDims,integer,decimal,complement,nParticles);
	struct timespec t5 = tprobe("Recursive method improved");

	// tprobe(NULL);
	// double result5 = interpAllGray2(val,pos,nGPointsProd,nDims,integer,decimal,complement,nCorners,weight,nParticles);
	// struct timespec t5 = tprobe("Gray method 2");

	printf("result0=%f\n",result0);
	printf("result1=%f\n",result1);
	printf("result2=%f\n",result2);
	printf("result3=%f\n",result3);
	printf("result4=%f\n",result4);
	printf("result5=%f\n",result5);

	// double ratio;
	// ratio = (double)t2.tv_nsec/t1.tv_nsec;
	// printf("ratio rec/fixed: %.2f\n",ratio);
	// ratio = (double)t3.tv_nsec/t1.tv_nsec;
	// printf("ratio gray/fixed: %.2f\n",ratio);

	printf("Interpolation to center:\n");

	tprobe(NULL);
	result1 = interpCAll3D(val,pos,nGPointsProd,nDims,nParticles);
	t1 = tprobe("Fixed 3D method");

	tprobe(NULL);
	result2 = interpCAllRecursive(val,pos,nGPointsProd,nDims,nParticles);
	t2 = tprobe("Recursive method");

	tprobe(NULL);
	result3 = interpCAllBinary(val,pos,nGPointsProd,nDims,nCorners,nParticles);
	t3 = tprobe("Binary method");

	tprobe(NULL);
	result4 = interpCAllGray(val,pos,nGPointsProd,nDims,nCorners,nParticles);
	t4 = tprobe("Gray method");

	tprobe(NULL);
	result5 = interpCAllRecursive2(val,pos,nGPointsProd,nDims,nParticles);
	t5 = tprobe("Recursive method improved");

	printf("result1=%f\n",result1);
	printf("result2=%f\n",result2);
	printf("result3=%f\n",result3);
	printf("result4=%f\n",result4);
	printf("result5=%f\n",result5);
	//
	// ratio = (double)t2.tv_nsec/t1.tv_nsec;
	// printf("ratio rec/fixed: %.2f\n",ratio);
	// ratio = (double)t3.tv_nsec/t1.tv_nsec;
	// printf("ratio gray/fixed: %.2f\n",ratio);

}
