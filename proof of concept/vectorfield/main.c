#define _XOPEN_SOURCE 700
#define _GNU_SOURCE	// Gives us ffs() to "find first bit set"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define inline

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

/*
 * Alternative A:
 * Store all grid points one component at a time
 */

inline double interpA(const double *val, const double *pos, const long int *nGPointsProd){

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

void accA1(const double *E, const double *pos, double *vel, const long int *nGPointsProd, int nDims, long int nParticles){
	for(int d=0;d<nDims;d++){
		const double *posloc = pos;
		double *velloc = vel+d;
		for(long int i=0;i<nParticles;i++){
			*velloc += interpA(E,posloc,nGPointsProd);
			posloc += nDims;
			velloc += nDims;
		}
		E += nGPointsProd[nDims];
	}
}

void accA1b(const double *E, const double *pos, double *vel, const long int *nGPointsProd, int nDims, long int nParticles){
	for(int d=0;d<nDims;d++){
		for(long int i=0;i<nParticles;i++){
			vel[i*nDims+d] += interpA(&E[d*nGPointsProd[nDims]],&pos[i*nDims],nGPointsProd);
		}
	}
}

void accA1c(const double *E, const double *pos, double *vel, const long int *nGPointsProd, int nDims, long int nParticles){
	int eind = 0;
	for(int d=0;d<nDims;d++){
		int posind = 0;
		int vind = d;
		for(long int i=0;i<nParticles;i++){
			vel[vind] += interpA(&E[eind],&pos[posind],nGPointsProd);
			posind += nDims;
			vind += nDims;
		}
		eind += nGPointsProd[nDims];
	}
}


void accA2(const double *E, const double *pos, double *vel, const long int *nGPointsProd, int nDims, long int nParticles){

	const double *Ex = E;
	const double *Ey = E+nGPointsProd[nDims];
	const double *Ez = Ey+nGPointsProd[nDims];

	for(long int i=0;i<nParticles;i++){

		*(vel++) += interpA(Ex,pos,nGPointsProd);
		*(vel++) += interpA(Ey,pos,nGPointsProd);
		*(vel++) += interpA(Ez,pos,nGPointsProd);

		pos += nDims;
	}

}

inline void interpAb(double *result, const double *val, const double *pos, const long int *sizeProd){

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
	long int p 		= j + k*sizeProd[1] + l*sizeProd[2];
	long int pj 	= p + 1;//sizeProd[1];
	long int pk 	= p + sizeProd[1];
	long int pjk 	= pk + 1;//sizeProd[1];
	long int pl 	= p + sizeProd[2];
	long int pjl 	= pl + 1;//sizeProd[1];
	long int pkl 	= pl + sizeProd[1];
	long int pjkl 	= pkl + 1;//sizeProd[1];

	// Linear interpolation
	for(int g=0;g<3;g++){
		int v = g*sizeProd[3];
		result[g] +=	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
								+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
						+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
								+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );
	}
}

void accAb(const double *E, const double *pos, double *vel, const long int *sizeProd, int nDims, long int nParticles){

	for(long int i=0;i<nParticles;i++){
		interpAb(vel,E,pos,sizeProd);
		vel+=nDims;
		pos+=nDims;
	}

}

/*
 * Alternative B:
 * Store components lumped together one grid node at a time
 */

inline void interpB(double *result, const double *val, const double *pos, const long int *nGPointsProd){

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
	long int p 		= 3*(j + k*nGPointsProd[1] + l*nGPointsProd[2]);
	long int pj 	= p + 3;
	long int pk 	= p + 3*nGPointsProd[1];
	long int pjk 	= pk + 3;
	long int pl 	= p + 3*nGPointsProd[2];
	long int pjl 	= pl + 3;
	long int pkl 	= pl + 3*nGPointsProd[1];
	long int pjkl 	= pkl + 3;

	// Linear interpolation
	for(int v=0;v<3;v++)
		result[v] =	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
							+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
					+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
							+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

inline void interpBb(double *result, const double *val, const double *pos, const long int *sizeProd){

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
	long int p 		= j*3 + k*sizeProd[2] + l*sizeProd[3];
	long int pj 	= p + 3;//sizeProd[1];
	long int pk 	= p + sizeProd[2];
	long int pjk 	= pk + 3;//sizeProd[1];
	long int pl 	= p + sizeProd[3];
	long int pjl 	= pl + 3;//sizeProd[1];
	long int pkl 	= pl + sizeProd[2];
	long int pjkl 	= pkl + 3;//sizeProd[1];

	// Linear interpolation
	for(int v=0;v<3;v++)
		result[v] =	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
							+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
					+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
							+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

inline void interpBc(double *result, const double *val, const double *pos, const long int *sizeProd){

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
	long int p 		= j*3 + k*sizeProd[2] + l*sizeProd[3];
	long int pj 	= p + 3;//sizeProd[1];
	long int pk 	= p + sizeProd[2];
	long int pjk 	= pk + 3;//sizeProd[1];
	long int pl 	= p + sizeProd[3];
	long int pjl 	= pl + 3;//sizeProd[1];
	long int pkl 	= pl + sizeProd[2];
	long int pjkl 	= pkl + 3;//sizeProd[1];

	// Linear interpolation
	for(int v=0;v<3;v++)
		result[v] +=	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
								+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
						+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
								+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

void accB(const double *E, const double *pos, double *vel, const long int *nGPointsProd, int nDims, long int nParticles, double *result){

	for(long int i=0;i<nParticles;i++){
		interpB(result,E,pos,nGPointsProd);
		for(int d=0;d<nDims;d++){
			*(vel++) += result[d];
		}
		pos += nDims;
	}

}

void accBd(const double *E, const double *pos, double *vel, const long int *sizeProd, int nDims, long int nParticles){

	double result[3];
	// for(long int i=0;i<nParticles;i++){
	// 	interpB(result,E,pos,nGPointsProd);
	// 	for(int d=0;d<nDims;d++){
	// 		*(vel++) += result[d];
	// 	}
	// 	pos += nDims;
	// }

	for(long int i=0;i<nParticles;i++){
		interpBb(result,E,pos,sizeProd);
		for(int d=0;d<nDims;d++){
			*(vel++) += result[d];
		}
		pos += nDims;
	}

}


void accBb(const double *E, const double *pos, double *vel, const long int *sizeProd, int nDims, long int nParticles, double *result){

	for(long int i=0;i<nParticles;i++){
		interpBb(result,E,pos,sizeProd);
		for(int d=0;d<nDims;d++){
			*(vel++) += result[d];
		}
		pos += nDims;
	}

}

void accBc(const double *E, const double *pos, double *vel, const long int *sizeProd, int nDims, long int nParticles){

	for(long int i=0;i<nParticles;i++){
		interpBc(vel,E,pos,sizeProd);
		vel += nDims;
		pos += nDims;
	}

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
	nGPoints[0] = 3;
	nGPoints[1] = nGPoints[0];
	nGPoints[2] = nGPoints[0];

	long int *nGPointsProd = malloc((nDims+1)*sizeof(long int));
	nGPointsProd[0] = 1;
	for(int d=1;d<=nDims;d++)
		nGPointsProd[d] = nGPointsProd[d-1] * nGPoints[d-1];

	int rank = nDims+1;
	int *size = malloc(rank*sizeof(*size));
	size[0] = 3;
	size[1] = 3;
	size[2] = size[1];
	size[3] = size[1];

	long int *sizeProd = malloc((rank+1)*sizeof(*sizeProd));
	sizeProd[0] = 1;
	for(int d=1;d<=rank;d++){
		sizeProd[d] = sizeProd[d-1] * size[d-1];
	}

	printf("Assigned auxiliary variables\n");

	double *EA = malloc(nDims*nGPointsProd[nDims]*sizeof(*EA));
	double *EB = malloc(nDims*nGPointsProd[nDims]*sizeof(*EB));
	for(long int p=0;p<nGPointsProd[nDims];p++){
		for(int d=0;d<nDims;d++){
			EA[p+d*nGPointsProd[nDims]] = p+(double)d/10;
			EB[p*nDims+d]               = p+(double)d/10;
		}
	}

	printf("Created fields\n");

	long int nParticles = 1e7;
	double *pos = malloc(nDims*nParticles*sizeof(*pos));
	double *velA1 = malloc(nDims*nParticles*sizeof(*velA1));
	double *velA1b = malloc(nDims*nParticles*sizeof(*velA1b));
	double *velA1c = malloc(nDims*nParticles*sizeof(*velA1c));
	double *velA2 = malloc(nDims*nParticles*sizeof(*velA2));
	double *velB = malloc(nDims*nParticles*sizeof(*velB));
	double *velBb = malloc(nDims*nParticles*sizeof(*velBb));
	double *velBc = malloc(nDims*nParticles*sizeof(*velBc));
	double *velBd = malloc(nDims*nParticles*sizeof(*velBd));


	double *velAb = malloc(nDims*nParticles*sizeof(*velAb));

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	for(long int i=0;i<nParticles*nDims;i++){
		pos[i] = gsl_rng_uniform(rng)*(nGPoints[i%nDims]-1);
		velA1[i] = 0;//gsl_rng_uniform(rng)*2-1;
		velA1b[i] = velA1[i];
		velA1c[i] = velA1[i];
		velA2[i] = velA1[i];
		velB[i] = velA1[i];
		velBb[i] = velA1[i];
		velBc[i] = velA1[i];
		velBd[i] = velA1[i];

		velAb[i] = velA1[i];

	}
	gsl_rng_free(rng);

	printf("Created velocities\n");

	/*
	 * TESTING ALGORITHMS
	 */

	double *buffer = malloc(nDims*sizeof(*buffer));

	// for(int p=0;p<nDims*nGPointsProd[nDims];p++){
	// 	printf("EA(%i)=%4.2f\t",p,EA[p]);
	// 	printf("EB(%i)=%4.2f\n",p,EB[p]);
	// }

	tprobe(NULL);

	accA1(EA,pos,velA1,nGPointsProd,nDims,nParticles);
	tprobe("accA1");
	tprobe(NULL);


	accA1b(EA,pos,velA1b,nGPointsProd,nDims,nParticles);
	tprobe("accA1b");
	tprobe(NULL);


	accA1c(EA,pos,velA1c,nGPointsProd,nDims,nParticles);
	tprobe("accA1c");
	tprobe(NULL);


	accA2(EA,pos,velA2,nGPointsProd,nDims,nParticles);
	tprobe("accA2");
	tprobe(NULL);


	accB(EB,pos,velB,nGPointsProd,nDims,nParticles,buffer);
	tprobe("accB");
	tprobe(NULL);

	accBb(EB,pos,velBb,sizeProd,nDims,nParticles,buffer);
	tprobe("accBb");
	tprobe(NULL);

	accBc(EB,pos,velBc,sizeProd,nDims,nParticles);
	tprobe("accBc");
	tprobe(NULL);

	accAb(EA,pos,velAb,nGPointsProd,nDims,nParticles);
	tprobe("accAb");
	tprobe(NULL);

	accBd(EB,pos,velBd,sizeProd,nDims,nParticles);
	tprobe("accBd");
	tprobe(NULL);

	for(int i=0;i<3&&i<nParticles;i++){
		printf("velA1 [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velA1[3*i],velA1[3*i+1],velA1[3*i+2]);
		printf("velA1b[%i]=(%12.2f,%12.2f,%12.2f)\n",i,velA1b[3*i],velA1b[3*i+1],velA1b[3*i+2]);
		printf("velA1c[%i]=(%12.2f,%12.2f,%12.2f)\n",i,velA1c[3*i],velA1c[3*i+1],velA1c[3*i+2]);
		printf("velA2 [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velA2[3*i],velA2[3*i+1],velA2[3*i+2]);
		printf("velB  [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velB[3*i],velB[3*i+1],velB[3*i+2]);
		printf("velBb [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velBb[3*i],velBb[3*i+1],velBb[3*i+2]);
		printf("velBc [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velBc[3*i],velBc[3*i+1],velBc[3*i+2]);
		printf("velAb [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velAb[3*i],velAb[3*i+1],velAb[3*i+2]);
		printf("velBd [%i]=(%12.2f,%12.2f,%12.2f)\n",i,velBd[3*i],velBd[3*i+1],velBd[3*i+2]);
	}
}
