/**
 * @file		population.c
 * @brief		Particle handling.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#define _XOPEN_SOURCE 700

#include "core.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <hdf5.h>
#include "iniparser.h"




/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

Population *pAlloc(const dictionary *ini,const MpiInfo *mpiInfo){

	// Get MPI info
	int size;//, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int mpiRank = mpiInfo->mpiRank;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

	// Load data
	int nSpecies = iniGetInt(ini,"population:nSpecies");
	int nDims = iniGetInt(ini,"grid:nDims");
	char **boundaries = iniGetStrArr(ini, "grid:boundaries" , 2*nDims);

	// Number of particles to allocate for (for all computing nodes)
	long int *nAllocTotal = iniGetLongIntArr(ini,"population:nAlloc",nSpecies);

	// Determine memory to allocate for this node
	long int *nAlloc = malloc(nSpecies*sizeof(long int));
	for(int s=0;s<nSpecies;s++){
		nAlloc[s] = ceil((double)nAllocTotal[s]/size);
		if(nAlloc[s]*size != nAllocTotal[s])
			msg(WARNING,"increased number of allocated particles from %i to %i"
			 			"to get integer per computing node",
						nAllocTotal[s], nAlloc[s]*size);
	}

	long int *iStart = malloc((nSpecies+1)*sizeof(long int));
	long int *iStop  = malloc(nSpecies*sizeof(long int));

	iStart[0] = 0;
	for(int s=1;s<nSpecies+1;s++) iStart[s]=iStart[s-1]+nAlloc[s-1];
	for(int s=0;s<nSpecies;s++) iStop[s]=iStart[s]; // No particles yet




	int rank = nDims+1;
	bndType *bnd = malloc(2*rank*sizeof(*bnd));

	bnd[0] = NONE; //initialize
	bnd[rank] = NONE; //initialize

	int b = 0;
	for (int d = 1; d<rank;d++){
		//msg(STATUS,"b=%i, d = %i, rank = %i",b,d,rank);
		int dd = d - 1;
		//msg(STATUS,"mpiRank = %i, subdomain[dd] = %i nSubdomainsProd[dd] = %i",mpiRank, subdomain[dd], nSubdomainsProd[dd]);
		int firstElem = mpiRank - subdomain[dd]*nSubdomainsProd[dd];

		// msg(STATUS,"lowerSubdomain: %i",lowerSubdomain);

		int upperSubdomain = firstElem
			+ ((subdomain[dd] + 1)%nSubdomains[dd])*nSubdomainsProd[dd];
		int lowerSubdomain = firstElem
			+ ((subdomain[dd] - 1 + nSubdomains[dd])%nSubdomains[dd])*nSubdomainsProd[dd];

		//printf("rank: %i, lowerSubdomain: %i, upperSubdomain = %i \n",mpiRank,lowerSubdomain,upperSubdomain);

		//lower
		//for(int r=dd; r<dd+1; r++){
		//printf("b = %i \n",b);
		int r=dd+1;
		if(lowerSubdomain>=mpiRank){
				if(		!strcmp(boundaries[b], "PERIODIC"))		bnd[r] = PERIODIC;
				else if(!strcmp(boundaries[b], "DIRICHLET"))	bnd[r] = DIRICHLET;
				else if(!strcmp(boundaries[b], "NEUMANN"))		bnd[r] = NEUMANN;
				else msg(ERROR,"%s invalid value for grid:boundaries",boundaries[b]);

			}else if (lowerSubdomain<mpiRank){
				//printf("YEPS!");
				bnd[r] = PERIODIC;
			} else {
				bnd[r] = NONE; //initialize
				//printf("bnd[%i] = NONE",r);

		}
		//upper
		//for(int r=rank+dd; r<rank+dd+1; r++){
		r = rank+dd+1;
		//printf("r = %i \n",r);
		if(upperSubdomain<=mpiRank){
				if(		!strcmp(boundaries[b+rank-1], "PERIODIC"))		bnd[r] = PERIODIC;
				else if(!strcmp(boundaries[b+rank-1], "DIRICHLET"))	bnd[r] = DIRICHLET;
				else if(!strcmp(boundaries[b+rank-1], "NEUMANN"))		bnd[r] = NEUMANN;
				else msg(ERROR,"%s invalid value for grid:boundaries",boundaries[b]);

			}else if(upperSubdomain>mpiRank){
				bnd[r] = PERIODIC;
			}else {
				bnd[r] = NONE; //initialize
				//printf("bnd[%i] = NONE",r);
			}
			b++;
		//}
	}
	//msg(ERROR,"%d, %d, %d, %d, %d, %d,%d, %d, ",bnd[0], bnd[1],bnd[2],bnd[3],bnd[4],bnd[5],bnd[6],bnd[7]);
	//printf("in pop rank = %i, %d, %d, %d, %d, %d, %d,%d, %d \n",mpiRank,bnd[0], bnd[1],bnd[2],bnd[3],bnd[4],bnd[5],bnd[6],bnd[7]);

	Population *pop = malloc(sizeof(Population));
	pop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->nSpecies = nSpecies;
	pop->nDims = nDims;
	pop->iStart = iStart;
	pop->iStop = iStop;
	pop->objVicinity = malloc(iStart[nSpecies]*sizeof(long int));
	pop->collisions = malloc(iStart[nSpecies]*sizeof(long int)); //malloc(sizeof pop->collisions)
	pop->kinEnergy = malloc((nSpecies+1)*sizeof(double));
	pop->TemperatureX = malloc((nSpecies+1)*sizeof(double));
	pop->TemperatureY = malloc((nSpecies+1)*sizeof(double));
	pop->TemperatureZ = malloc((nSpecies+1)*sizeof(double));
	pop->TemperatureTot = malloc((nSpecies+1)*sizeof(double));
	pop->potEnergy = malloc((nSpecies+1)*sizeof(double));
	pop->charge = iniGetDoubleArr(ini,"population:charge",nSpecies);
	pop->mass = iniGetDoubleArr(ini,"population:mass",nSpecies);
	pop->bnd = bnd;

	free(nAlloc);
	free(nAllocTotal);

	return pop;

}

void pFree(Population *pop){

	free(pop->pos);
	free(pop->vel);
	free(pop->kinEnergy);
	free(pop->potEnergy);
	free(pop->TemperatureX);
	free(pop->TemperatureY);
	free(pop->TemperatureZ);
	free(pop->TemperatureTot);
	free(pop->iStart);
	free(pop->iStop);
	free(pop->objVicinity);
	free(pop->collisions);
	free(pop->charge);
	free(pop->mass);
	free(pop->bnd);
	free(pop);

}

void pPosUniform(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const gsl_rng *rng){

	// Read from ini
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = gGetGlobalSize(ini);

	for(int s=0;s<nSpecies;s++){

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart;
		double *pos = &pop->pos[iStart*nDims];

		// Iterate through all particles to be generated. Same seed on all MPI
		// nodes ensure same particles are generated everywhere.
		for(long int i=0;i<nParticles[s];i++){

			// Generate position for particle i
			for(int d=0;d<nDims;d++){
				pos[d] = L[d]*gsl_rng_uniform_pos(rng);
				//if(pos[d]<1) printf("\n \n pos[d] = %f \n \n ",pos[d]);
			}

			// Count the number of dimensions where the particle resides in
			// the range of this node
			int correctRange = 0;
			for(int d=0;d<nDims;d++)
				correctRange += (subdomain[d] == (int)(posToSubdomain[d]*pos[d]));

			// Iterate only if particle resides in this sub-domain.
			if(correctRange==nDims){
				pos += nDims;
				iStop++;
			}

		}

		if(iStop>pop->iStart[s+1]){
			int allocated = pop->iStart[s+1]-iStart;
			int generated = iStop-iStart;
			msg(ERROR,	"allocated only %i particles of specie %i per node but"
			 			"%i generated", allocated, s, generated);
		}

		pop->iStop[s]=iStop;

	}

	pToLocalFrame(pop,mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);

}

void pPosLattice(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo){

	// Read from ini
	int nDims = pop->nDims;
	int nSpecies = pop->nSpecies;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = gGetGlobalSize(ini);
	long int V = gGetGlobalVolume(ini);

	for(int s=0;s<nSpecies;s++){

		// Particle-particle distance in lattice
		double l = pow(V/(double)nParticles[s],1.0/nDims);

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart;
		double *pos = &pop->pos[iStart*nDims];

		// Iterate through all particles to be generated
		// Generate particles on global frame on all nodes and discard the ones
		// out of range. This is simpler as it resembles pPosUniform()
		for(long int i=0;i<nParticles[s];i++){

			double linearPos = l*i;
			for(int d=0;d<nDims;d++){
				pos[d] = fmod(linearPos,L[d]);
				linearPos /= L[d];
			}

			// Count the number of dimensions where the particle resides in
			// the range of this node
			int correctRange = 0;
			for(int d=0;d<nDims;d++)
				correctRange += (subdomain[d] == (int)(posToSubdomain[d]*pos[d]));

			// Iterate only if particle resides in this sub-domain.
			if(correctRange==nDims){
				pos += nDims;
				iStop++;
			}

		}

		if(iStop>pop->iStart[s+1]){
			int allocated = pop->iStart[s+1]-iStart;
			int generated = iStop-iStart;
			msg(ERROR,	"allocated only %i particles of specie %i per node but"
			 			"%i generated", allocated, s, generated);
		}

		pop->iStop[s]=iStop;

	}

	pToLocalFrame(pop,mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);

}

void pPosPerturb(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo){

	int nDims = pop->nDims;
	int nSpecies = pop->nSpecies;

	int nElements = nDims *nSpecies;
	double *amplitude = iniGetDoubleArr(ini,"population:perturbAmplitude",nElements);
	double *mode = iniGetDoubleArr(ini,"population:perturbMode",nElements);

	int *L = gGetGlobalSize(ini);
	double *pos = pop->pos;


	pToGlobalFrame(pop,mpiInfo);

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];
		for(long int i=iStart;i<iStop;i++){

			for(int d=0;d<nDims;d++){
				double theta = 2.0*M_PI*mode[s*nDims+d]*pos[i*nDims+d]/L[d];
				pos[i*nDims+d] += amplitude[s*nDims+d]*cos(theta);
			}
		}
	}

	pToLocalFrame(pop,mpiInfo);

	free(L);
	free(amplitude);
	free(mode);

}

void pPosDebug(const dictionary *ini, Population *pop){

	int nSpecies = pop->nSpecies;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);

	int mpiRank, mpiSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	for(int s=0;s<nSpecies;s++){
		nParticles[s] /= mpiSize;
	}

	int nDims = pop->nDims;
	long int *nMigrantsResult = malloc(81*sizeof(*nMigrantsResult));
	alSet(nMigrantsResult,81,
			1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,
			1,1,0,1,1,0,1,1,0,3,3,0,0,0,0,4,4,0,1,1,0,1,1,0,1,1,0,
			1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0);

	for(int s=0;s<nSpecies;s++){
		long int iStart = pop->iStart[s];
		pop->iStop[s] = iStart + nParticles[s];
		double *pos = &pop->pos[iStart*nDims];

		for(long int i=0;i<nParticles[s];i++){
			for(int d=0;d<nDims;d++){
				pos[i*nDims+d] = 1000*mpiRank + i + (double)d/10 + (double)s/100;
			}
		}
	}

	free(nParticles);

}

// Verify that position is indeed within local frame and will not produce any
// segmentation faults
void pPosAssertInLocalFrame(const Population *pop, const Grid *grid){

	int *size = grid->size;
	double *pos = pop->pos;

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0; s<nSpecies; s++){

		long int iStart = pop->iStart[s];
		long int iStop  = pop->iStop[s];
		for(int i=iStart; i<iStop; i++){

			for(int d=0; d<nDims; d++){

				if(pos[i*nDims+d]>size[d+1]-1 || pos[i*nDims+d]<0){
					msg(ERROR,	"Particle i=%li (of specie %i) is out of bounds"
					 			"in dimension %i: %f>%i",
								i, s, d, pos[i*nDims+d], size[d+1]-1);
				}
			}
		}
	}
}

void pVelAssertMax(const Population *pop, double max){

	double *vel = pop->vel;

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0; s<nSpecies; s++){

		long int iStart = pop->iStart[s];
		long int iStop  = pop->iStop[s];
		for(int i=iStart; i<iStop; i++){

			for(int d=0;d<nDims;d++){

				if(vel[i*nDims+d]>max){
					msg(ERROR,	"Particle i=%li (of specie %i) travels too"
					 			"fast in dimension %i: %f>%f",
								i, s, d, vel[i*nDims+d], max);
				}
			}
		}
	}
}

void pVelMaxwell(const dictionary *ini, Population *pop, const gsl_rng *rng){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	long int index = 0;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = velThermal[s];
		velTh = velTh;

		for(long int i=iStart;i<iStop;i++){

			double *vel = &pop->vel[i*nDims];
			for(int d=0;d<nDims;d++){
				index = (s*nDims)+d;
				vel[d] = velDrift[index] + gsl_ran_gaussian_ziggurat(rng,velTh);
			}
		}
	}
	free(velDrift);
	free(velThermal);
}

void pPurgeGhost(Population *pop, const Grid *grid){

	//this will delete all particles reciding on Dirichlet boundarys
	int *size = grid->size;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	double *pos = pop->pos;
	double *vel = pop->vel;
	bool cut = false;
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0; s<nSpecies; s++){

		long int iStart = pop->iStart[s];
		long int iStop  = pop->iStop[s];
		for(int i=iStart; i<iStop; i++){
			//printf("iStop= %i\n",iStop);

			//msg(STATUS,"bnd[0+1] = %d, bnd[1+1] = %d, bnd[2+1] = %d \n",bnd[0+1],bnd[1+1],bnd[2+1]);
			//msg(STATUS,"bnd[rank+0+1] = %d, bnd[rank+1+1] = %d, bnd[rank+2+1] = %d \n",bnd[rank+0+1],bnd[rank+1+1],bnd[rank+2+1]);
			for(int d=0; d<nDims; d++){
				//printf("nGhostLayers[d] %i \n",nGhostLayers[d+1]);
				if( (pos[i*nDims+d]>size[d+1]-nGhostLayers[d+1]-1 && bnd[d+rank+1]==DIRICHLET) ){
					//msg(STATUS,"CUT: pos[i*nDims+d] = %f, bnd[d+1]=%d \n",pos[i*nDims+d],bnd[d+1]);
					cut = true;
				}
				if( (pos[i*nDims+d]<nGhostLayers[d+rank+1] && bnd[d+1]==DIRICHLET) ){
					//msg(STATUS,"CUT: pos[i*nDims+d] = %f, bnd[d+rank+1]=%d \n",pos[i*nDims+d],bnd[d+rank+1]);
					cut = true;
				}
			}
			if (cut == true){
				pCut(pop, s, i*nDims, pos, vel);
				cut = false;
				iStop--;
				i--;
				//printf("iStop= %i \n \n",iStop);
			}
		}
	}
}

// void pFillGhost(const dictionary *ini, Population *pop, const gsl_rng *rng, const MpiInfo *mpiInfo){
//
// 	int nSpecies = pop->nSpecies;
// 	int nDims = pop->nDims;
// 	bndType *bnd = pop->bnd;
// 	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
// 	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
// 	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
// 	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);
// 	int *nGhostLayers = iniGetIntArr(ini,"grid:nGhostLayers",2*nDims);
// 	//double timeStep = iniGetDouble(ini, "time:timeStep");
// 	int *L = gGetGlobalSize(ini);
// 	int rank = nDims+1;
//
// 	// Read from mpiInfo
// 	int *subdomain = mpiInfo->subdomain;
// 	double *posToSubdomain = mpiInfo->posToSubdomain;
//
// 	long int index = 0;
// 	double pos[nDims];
// 	double vel[nDims];
// 	int edge[nDims];
//
// 	for(int s=0;s<nSpecies;s++){
//
// 		long int iStart = pop->iStart[s];
// 		long int iStop = pop->iStop[s];
//
// 		double velTh = velThermal[s];
//
// 		for (int d=0;d<nDims;d++){
// 			index = (s*nDims)+d;
//
// 			//compute edge
// 			for (int dd=0;dd<nDims;dd++){
// 				if (velDrift[index]/velDrift[dd] == 1) edge[dd] = 0;
// 				else edge[dd] = 1;
// 			}
//
// 			long int globalSizeProd = (L[0]*L[1]*L[2]); //TODO: make nDimensional
//
// 			long int sliceSize = (globalSizeProd/(L[d]));
// 			long int newParticles = ( (sliceSize*nParticles[s])/globalSizeProd ); // slice * particles per cell
// 			printf("generating %li particles for specie %i, velTh = %f \n \n",newParticles,s,velTh);
// 			//for(long int i=iStart;i<iStop;i++){
// 			//for(long int i=0;i<newParticles;i++){
//
// 				//Lower ghost slice
// 			if(bnd[d+1]==DIRICHLET){
// 				for(long int i=0;i<newParticles;i++){
// 					//generate velocity for particle
// 					for(int d=0;d<nDims;d++){
// 						vel[d] = velDrift[(s*nDims)+d] + gsl_ran_gaussian_ziggurat(rng,velTh);
// 					}
//
// 					// Generate position for particle
// 					for(int dd=0;dd<nDims;dd++){
// 						pos[dd] = 1+(L[dd])*gsl_rng_uniform_pos(rng);
// 					}
// 					pos[d] = nGhostLayers[d+1]*(gsl_rng_uniform_pos(rng)); //in lower ghost
//
// 					int correctRange = 0;
// 					for(int dd=0;dd<nDims;dd++)
// 						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*pos[dd]));
//
// 					// Add only if particle resides in this sub-domain.
// 					if(correctRange==nDims){
// 						pNew(pop,s,pos,vel);
// 					}
// 				}
// 			}
//
// 				//Upper ghost slice
// 			if(bnd[d+rank+1]==DIRICHLET){
// 				for(long int i=0;i<newParticles;i++){
// 					//generate velocity for particle
//
// 					for(int d=0;d<nDims;d++){
// 						vel[d] = velDrift[(s*nDims)+d] + gsl_ran_gaussian_ziggurat(rng,velTh);
// 					}
//
// 					// Generate position for particle
// 					for(int dd=0;dd<nDims;dd++){
// 						pos[dd] = (L[dd])*gsl_rng_uniform_pos(rng);
// 					}
// 					pos[d] = trueSize[d]+nGhostLayers[d+1]*(gsl_rng_uniform_pos(rng)); //in lower ghost
//
// 					int correctRange = 0;
// 					for(int dd=0;dd<nDims;dd++)
// 						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*(pos[dd]-1)));
//
// 					// Add only if particle resides in this sub-domain.
// 					if(correctRange==nDims){
// 						pNew(pop,s,pos,vel);
//
// 					}
// 				}
//
//
//
// 			}
// 		}
// 	}
// 	free(velDrift);
// 	free(velThermal);
// 	free(trueSize);
// 	free(nParticles);
// 	free(nGhostLayers);
//
// 	return;
// }

void pFillGhost(const dictionary *ini, Population *pop, const gsl_rng *rng, const MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	bndType *bnd = pop->bnd;
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);
	int *nGhostLayers = iniGetIntArr(ini,"grid:nGhostLayers",2*nDims);
	//double timeStep = iniGetDouble(ini, "time:timeStep");
	int *L = gGetGlobalSize(ini);
	int rank = nDims+1;

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	long int index = 0;
	double pos[nDims];
	double vel[nDims];
	// int edge[nDims];

	pToGlobalFrame(pop,mpiInfo);

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = velThermal[s];

		for (int d=0;d<nDims;d++){
			index = (s*nDims)+d;

			//compute edge
			// for (int dd=0;dd<nDims;dd++){
			// 	if (velDrift[index]/velDrift[dd] == 1) edge[dd] = 0;
			// 	else edge[dd] = 1;
			// }

			long int globalSizeProd = (L[0]*L[1]*L[2]); //TODO: make nDimensional

			long int sliceSize = (globalSizeProd/(L[d]));
			long int newParticles = ( (sliceSize*nParticles[s])/globalSizeProd ); // slice * particles per cell
			//printf("generating %li particles for specie %i, velTh = %f \n \n",newParticles,s,velTh);
			//for(long int i=iStart;i<iStop;i++){
			//for(long int i=0;i<newParticles;i++){

				//Lower ghost slice
			if(bnd[d+1]==DIRICHLET){
				for(long int i=0;i<newParticles;i++){
					//generate velocity for particle
					for(int d=0;d<nDims;d++){
						vel[d] = velDrift[(s*nDims)+d] + gsl_ran_gaussian_ziggurat(rng,velTh);
					}

					// Generate position for particle
					for(int dd=0;dd<nDims;dd++){
						pos[dd] = (L[dd])*gsl_rng_uniform_pos(rng)-nGhostLayers[d]*0.5;
					}
					pos[d] = nGhostLayers[d]*(gsl_rng_uniform_pos(rng))-nGhostLayers[d]; //in lower ghost
					//printf("nGhostLayers[d+1] = %i\n",nGhostLayers[d]);
					int correctRange = 0;
					for(int dd=0;dd<nDims;dd++){
						//printf("posToSubdomain[dd] = %f, dd = %i \n",posToSubdomain[dd],dd);
						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*pos[dd]));
					}
					// Add only if particle resides in this sub-domain.
					if(correctRange==nDims){
						if((mpiInfo->mpiRank)==4){
							//printf("adding to pos: %f,%f,%f \n",pos[0],pos[1],pos[2]);
						}
						pNew(pop,s,pos,vel);
					}
				}
			}

				//Upper ghost slice
			if(bnd[d+rank+1]==DIRICHLET){
				for(long int i=0;i<newParticles;i++){
					//generate velocity for particle

					for(int d=0;d<nDims;d++){
						vel[d] = velDrift[(s*nDims)+d] + gsl_ran_gaussian_ziggurat(rng,velTh);
					}

					// Generate position for particle
					for(int dd=0;dd<nDims;dd++){
						pos[dd] = (L[dd])*gsl_rng_uniform_pos(rng)-nGhostLayers[d]*0.5;

					}
					pos[d] = L[d]+nGhostLayers[d]*(gsl_rng_uniform_pos(rng))-nGhostLayers[d]; //in lower ghost

					int correctRange = 0;
					for(int dd=0;dd<nDims;dd++)
						correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*(pos[dd])));

					// Add only if particle resides in this sub-domain.
					if(correctRange==nDims){
						pNew(pop,s,pos,vel);

					}
				}



			}
		}
	}
	pToLocalFrame(pop,mpiInfo);

	free(velDrift);
	free(velThermal);
	free(trueSize);
	free(nParticles);
	free(nGhostLayers);
	free(L);

	return;
}

//Depricated: Not the best way to add a particle flux
void pInfluxDrift(const dictionary *ini, Population *pop, const gsl_rng *rng, const MpiInfo *mpiInfo){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);
	int *nGhostLayers = iniGetIntArr(ini,"grid:nGhostLayers",2*nDims);
	//double timeStep = iniGetDouble(ini, "time:timeStep");
	int *L = gGetGlobalSize(ini);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	long int index = 0;
	double pos[nDims];
	double vel[nDims];
	int edge[nDims];

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = velThermal[s];

		for (int d=0;d<nDims;d++){
			index = (s*nDims)+d;

			//compute edge
			for (int dd=0;dd<nDims;dd++){
				if (velDrift[index]/velDrift[dd] == 1) edge[dd] = 0;
				else edge[dd] = 1;
			}


			// long int globalSizeProd = (( (nGhostLayers[0]+L[0]) )
			// 	*( (nGhostLayers[1]+L[1]))
			// 	*( (nGhostLayers[2]+L[2]) )); //TODO: make nDimensional
			//
			// long int sliceSize = (globalSizeProd/(nGhostLayers[d]+L[d]));


			long int globalSizeProd = (L[0]*L[1]*L[2]); //TODO: make nDimensional

			//printf("L[0]= %i, L[1]= %i, L[2]= %i \n",nGhostLayers[0]+L[0],L[1],L[2]);
			long int sliceSize = (globalSizeProd/(L[d]));
			long int newParticles = ( (sliceSize*nParticles[s]/(L[0]*L[1]*L[2]))*(velDrift[index]+(!edge[d])*(velTh)*velDrift[index]));
				//+(!edge[d])*(1.414*nDims*velDrift[index]*(velTh-velDrift[index]))*(velTh) ) );
			//need check if veldrift is < 0
			//printf("%f,%li \n",velDrift[index],sliceSize);
			//printf("generating %li particles for specie %i, velTh = %f \n \n",newParticles,s,velTh);
			//for(long int i=iStart;i<iStop;i++){
			for(long int i=0;i<newParticles;i++){
				//double *vel = &pop->vel[i*nDims];

				for(int d=0;d<nDims;d++){
					//index = (s*nDims)+d;
					vel[d] = velDrift[(s*nDims)+d] + gsl_ran_gaussian_ziggurat(rng,velTh);
					//printf("vel[d] = %f, d= %i\n",vel[d],d);
				}

				// Generate position for particle
				for(int dd=0;dd<nDims;dd++){
					pos[dd] = edge[dd]*(L[dd])*gsl_rng_uniform_pos(rng);

					//printf("pos[%i] = %f\n",dd,pos[dd]);
				}
				pos[d] = 1+velDrift[index]*(gsl_rng_uniform_pos(rng)); // pos in local frame
				//printf("pos[%i] = %f\n",d,pos[d]);
				// Count the number of dimensions where the particle resides in
				// the range of this node
				int correctRange = 0;
				for(int dd=0;dd<nDims;dd++)
					correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*pos[dd]));

				// Add only if particle resides in this sub-domain.
				if(correctRange==nDims){
					//printf("pos[%i] = %f\n",d,pos[d]);
					//printf("%f,%li \n",velDrift[index],sliceSize);
					pNew(pop,s,pos,vel);
				}
			}
		}
	}
	free(velDrift);
	free(velThermal);

	return;
}

void pVelSet(Population *pop, const double *vel){

	int nDims = pop->nDims;
	int nSpecies = pop->nSpecies;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){
			for(int d=0;d<nDims;d++){
				pop->vel[i*nDims+d] = vel[d];
			}
		}
	}
}

void pVelZero(Population *pop){

	int nDims = pop->nDims;
	int nSpecies = pop->nSpecies;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){
			for(int d=0;d<nDims;d++){
				pop->vel[i*nDims+d] = 0;
			}
		}
	}
}

void pVelConstant(const dictionary *ini, Population *pop, double constant1, double constant2){

	//test function. takes only two species
	int nDims = pop->nDims;
	//int nSpecies = pop->nSpecies;
	//double timeStep = iniGetDouble(ini,"time:timeStep");
	//double stepSize = iniGetDouble(ini,"grid:stepSize");

	long int iStart1 = pop->iStart[0];
	long int iStop1 = pop->iStop[0];
	long int iStart2 = pop->iStart[1];
	long int iStop2 = pop->iStop[1];

	for(long int i=iStart1;i<iStop1;i++){
		//for(int d=0;d<nDims;d++){
			//pop->vel[i*nDims] = constant1;
			for(int d=0;d<nDims;d++){
				pop->vel[i*nDims+d] = constant1;
			//	msg(STATUS, "vel %i = %f", d, pop->vel[i*nDims+d]);
		}
		//}
	}
	for(long int i=iStart2;i<iStop2;i++){
		//for(int d=0;d<nDims;d++){
			//pop->vel[i*nDims] = constant2;
			for(int d=0;d<nDims;d++){
				pop->vel[i*nDims+d] = constant2;
			//	msg(STATUS, "vel %i = %f", d, pop->vel[i*nDims+d]);
		}

	}
}

void pNew(Population *pop, int s, const double *pos, const double *vel){

	int nDims = pop->nDims;
	long int *iStart = pop->iStart;
	long int *iStop = pop->iStop;	// New particle added here

	if(iStop[s]>=iStart[s+1])
		msg(WARNING,"Not enough allocated memory to add new particle to specie"
		 			"%i. New particle ignored.",s);
	else {

		long int p = iStop[s]*nDims;
		for(int d=0;d<nDims;d++){
			pop->pos[p+d] = pos[d];
			pop->vel[p+d] = vel[d];
		}


		iStop[s]++;

	}

}

void pCut(Population *pop, int s, long int p, double *pos, double *vel){

	int nDims = pop->nDims;
	long int pLast = (pop->iStop[s]-1)*nDims;
	//printf("Cutting particle p = %li \n",p);
	for(int d=0;d<nDims;d++){
		pos[d] = pop->pos[p+d];
		vel[d] = pop->vel[p+d];
		//printf("pos[%i] = %f \n",d,pos[d]);
		pop->pos[p+d] = pop->pos[pLast+d];
		pop->vel[p+d] = pop->vel[pLast+d];
	}

	pop->iStop[s]--;

}

funPtr pFindCollisionType(Population *pop, Object *obj, long int n){

	msg(WARNING, "Function to determine collision type not yet implemented!");

}

void pBackscatter(Population *pop, const Object *obj, long int n){

	msg(WARNING, "backscatter function not yet implemented!");
	msg(STATUS, "particle %i backscattered an electron \n", n);

}

void pSecondaryElectron(Population *pop, const Object *obj, long int n){

	msg(WARNING, "Secondary Electron function not yet implemented!");
	msg(STATUS, "particle %i created a secondary electron \n", n);

}

void pReflect(Population *pop, const Object *obj, long int n){
	
	double *newVel;
	double *surfNodes;
	double *norm;
	double *intersect;
	double *pos = &pop->pos[3*n];
	double *vel = &pop->vel[3*n];

	surfNodes = oFindNearestSurfaceNodes(pop, obj, n);

	adNormal(surfNodes, surfNodes + 3, norm, n);
	oFindIntersectPoint(pop, n, norm, surfNodes, intersect);
	adReflect(vel, surfNodes, surfNodes+1, newVel);

	msg(WARNING, "Reflection function not yet implemented!");
	msg(STATUS, "particle %i Reflected off an object \n", n);

}

void pPhotoElectron(Population *pop, const Object *obj, long int n){
	msg(WARNING, "photoelectron function not yet implemented!");
}

void pAdhere(Population *pop, const Object *obj, long int n){

	msg(WARNING, "Adhesion function not yet implemented!");
	msg(STATUS, "particle %i adhered to an object \n", n);
}

void pOpenH5(	const dictionary *ini, Population *pop, const Units *units,
	   			const char *fName){

	/*
	 * CREATE FILE
	 */
	hid_t file = openH5File(ini,fName,"pop");

	/*
	 * CREATE GROUPS
	 */
	hid_t group;
	group = H5Gcreate(file,"/pos",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Gclose(group);
	group = H5Gcreate(file,"/vel",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Gclose(group);

	char name[32];	// int is max 5 digits + "/pos/specie " + '\0'

	int nSpecies = pop->nSpecies;
	for(int s=0;s<nSpecies;s++){
		sprintf(name,"/pos/specie %i",s);
		group = H5Gcreate(file,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Gclose(group);

		sprintf(name,"/vel/specie %i",s);
		group = H5Gcreate(file,name,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Gclose(group);
	}

	pop->h5 = file;

	/*
	 * CREATE ATTRIBUTES
	 */

	setH5Attr(file, "Position denormalization factor", &units->length, 1);
	setH5Attr(file, "Velocity denormalization factor", &units->velocity, 1);

}

void pWriteH5(Population *pop, const MpiInfo *mpiInfo, double posN, double velN){

	int mpiRank = mpiInfo->mpiRank;
	int mpiSize = mpiInfo->mpiSize;
	int nSpecies = pop->nSpecies;

	pToGlobalFrame(pop,mpiInfo);

	/*
 	 * HDF5 HYPERSLAB DEFINITION
	 */
	const int arrSize = 2;
	hsize_t fileDims[arrSize];		// Size of data in file
	hsize_t memDims[arrSize];		// Size of data in memory of this MPI node
	hsize_t offset[arrSize];		// At which offset in file to store data of this MPI node
	fileDims[1] = pop->nDims;
	memDims[1] = pop->nDims;
	offset[1] = 0;

	long int *offsetAllSubdomains = malloc((mpiSize+1)*sizeof(long int));
	offsetAllSubdomains[0] = 0;

	for(int s=0;s<nSpecies;s++){

		long int nParticles = pop->iStop[s] - pop->iStart[s];
		MPI_Allgather(	&nParticles,
						1,
						MPI_LONG,
						&offsetAllSubdomains[1],
						1,
						MPI_LONG,
						MPI_COMM_WORLD);

		// Take cumulative sum to actually get offset
		// Last element equals total number of particles on all nodes
		for(int r=1;r<=mpiSize;r++){
			offsetAllSubdomains[r] += offsetAllSubdomains[r-1];
		}

		// Only proceed if there actually are any particles to save, at least
		// on one MPI node. HDF5 may crash otherwise.
		if(offsetAllSubdomains[mpiSize]){

			fileDims[0] = offsetAllSubdomains[mpiSize];
			memDims[0] = nParticles;
			offset[0] = offsetAllSubdomains[mpiRank];

			hid_t memSpace = H5Screate_simple(arrSize,memDims,NULL);
			hid_t fileSpace = H5Screate_simple(arrSize,fileDims,NULL);

			H5Sselect_hyperslab(fileSpace,
								H5S_SELECT_SET,
								offset,
								NULL,
								memDims,
								NULL);

			/*
			 * STORE DATA COLLECTIVELY
			 */
			hid_t pList = H5Pcreate(H5P_DATASET_XFER);
		    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);

			char name[64];
			hid_t dataset;

			sprintf(name,"/pos/specie %i/n=%.1f",s,posN);
			dataset = H5Dcreate(pop->h5,
								name,
								H5T_IEEE_F64LE,
								fileSpace,
								H5P_DEFAULT,
								H5P_DEFAULT,
								H5P_DEFAULT);

			H5Dwrite(	dataset,
						H5T_NATIVE_DOUBLE,
						memSpace,
						fileSpace,
						pList,
						&pop->pos[pop->iStart[s]*pop->nDims]);

			H5Dclose(dataset);

			sprintf(name,"/vel/specie %i/n=%.1f",s,velN);
			dataset = H5Dcreate(pop->h5,
								name,
								H5T_IEEE_F64LE,
								fileSpace,
								H5P_DEFAULT,
								H5P_DEFAULT,
								H5P_DEFAULT);

			H5Dwrite(	dataset,
				 		H5T_NATIVE_DOUBLE,
						memSpace,
						fileSpace,
						pList,
						&pop->vel[pop->iStart[s]*pop->nDims]);

			H5Dclose(dataset);

			H5Pclose(pList);
			H5Sclose(fileSpace);
			H5Sclose(memSpace);

		} else {
			msg(WARNING,"No particles of specie %i to store in .h5-file",s);
		}
	}
 	free(offsetAllSubdomains);

	pToLocalFrame(pop,mpiInfo);
}

void pCloseH5(Population *pop){
	H5Fclose(pop->h5);
}


void pCreateEnergyDatasets(hid_t xy, Population *pop){

	char name[64];
	int nSpecies = pop->nSpecies;

	sprintf(name,"/energy/potential/total");
	xyCreateDataset(xy,name);

	sprintf(name,"/energy/kinetic/total");
	xyCreateDataset(xy,name);

	for(int s=0;s<nSpecies;s++){
		sprintf(name,"/energy/potential/specie %i",s);
		xyCreateDataset(xy,name);

		sprintf(name,"/energy/kinetic/specie %i",s);
		xyCreateDataset(xy,name);
	}
}

void pCreateTemperatureDatasets(hid_t xy, Population *pop){

	char name[64];
	int nSpecies = pop->nSpecies;

	for(int s=0;s<nSpecies;s++){
		sprintf(name,"/energy/TemperatureTot/specie %i",s);
		xyCreateDataset(xy,name);

		sprintf(name,"/energy/TemperatureX/specie %i",s);
		xyCreateDataset(xy,name);

		sprintf(name,"/energy/TemperatureY/specie %i",s);
		xyCreateDataset(xy,name);

		sprintf(name,"/energy/TemperatureZ/specie %i",s);
		xyCreateDataset(xy,name);
	}
}

void pWriteEnergy(hid_t xy, Population *pop, double x,Units *units){

	char name[64];
	int nSpecies = pop->nSpecies;

	double denorm = units->energy;

	sprintf(name,"/energy/potential/total");
	xyWrite(xy,name,x,denorm*pop->potEnergy[nSpecies],MPI_SUM);

	sprintf(name,"/energy/kinetic/total");
	xyWrite(xy,name,x,denorm*pop->kinEnergy[nSpecies],MPI_SUM);

	for(int s=0; s<nSpecies; s++){

		sprintf(name,"/energy/potential/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->potEnergy[s],MPI_SUM);

		sprintf(name,"/energy/kinetic/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->kinEnergy[s],MPI_SUM);
	}

}

void pWriteTemperature(hid_t xy, Population *pop, double x,Units *units,dictionary *ini){

	// Temperature instead of kinetic energy
	double k_b = 1.38064852e-23;

	char name[64];
	int nSpecies = pop->nSpecies;
	double denorm = units->energy;
	int nDims = pop->nDims;
	int *nSubdomains = iniGetIntArr(ini,"grid:nSubdomains",nDims);
	int nProcs = nSubdomains[0]*nSubdomains[1]*nSubdomains[2];
	//msg(STATUS,"nsudims = %i",nSubdomains[2]);
	//long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",nSpecies);

	for(int s=0; s<nSpecies; s++){
		//msg(STATUS,"s = %i, Tx =%f, Tt = %f",s,pop->TemperatureX[s],pop->TemperatureTot[s]);
		sprintf(name,"/energy/TemperatureTot/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->TemperatureTot[s]/(nProcs*units->weights[s]*k_b),MPI_SUM);
		//msg(STATUS, " nprocs = %i temp = %f",nProcs,(denorm*pop->TemperatureX[s]/(nProcs*units->weights[s]*k_b)));
		// sprintf(name,"/energy/kinetic/specie %i",s);
		// xyWrite(xy,name,x,denorm*pop->kinEnergy[s]/(nParticles[s]*units->weights[s]*k_b),MPI_SUM);

		sprintf(name,"/energy/TemperatureX/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->TemperatureX[s]/(nProcs*units->weights[s]*k_b),MPI_SUM);

		sprintf(name,"/energy/TemperatureY/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->TemperatureY[s]/(nProcs*units->weights[s]*k_b),MPI_SUM);

		sprintf(name,"/energy/TemperatureZ/specie %i",s);
		xyWrite(xy,name,x,denorm*pop->TemperatureZ[s]/(nProcs*units->weights[s]*k_b),MPI_SUM);

	}

}

void pSumKinEnergy(Population *pop){

	int nSpecies = pop->nSpecies;

	pop->kinEnergy[nSpecies] = 0;
	for(int s=0; s<nSpecies; s++){
		pop->kinEnergy[nSpecies] += pop->kinEnergy[s];
	}

}

void pSumPotEnergy(Population *pop){

	int nSpecies = pop->nSpecies;

	pop->potEnergy[nSpecies] = 0;
	for(int s=0; s<nSpecies; s++){
		pop->potEnergy[nSpecies] += pop->potEnergy[s];
	}

}


/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

void pToLocalFrame(Population *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[i*nDims];
			for(int d=0;d<nDims;d++) pos[d] -= offset[d];
		}
	}
}

void pToGlobalFrame(Population *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[i*nDims];
			for(int d=0;d<nDims;d++) pos[d] += offset[d];
		}
	}
}
