/**
 * @file		particles.c
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Particle handling.
 * @date		26.10.15
 *
 * Functions for handling particles: initialization and finalization of
 * particle structs, reading and writing of data and so on.
 */

#include "pinc.h"
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <hdf5.h>
#include "iniparser.h"

Population *allocPopulation(const dictionary *ini){

	// Sanity check
	iniAssertEqualNElements(ini,4,"population:nParticles","population:nAlloc","population:q","population:m");

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);	// Presumes sanity check on nSubdomains by allocGrid()
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data
	int nSpecies;
	long int *nAllocTotal = iniGetLongIntArr(ini,"population:nAlloc",&nSpecies);	// This is for all computing nodes
	int nDims = iniGetNElements(ini,"grid:nTGPoints");
	if(nDims==0) msg(ERROR,"grid:nTGPoints not found");

	// Determine memory to allocate for this node
	long int *nAlloc = malloc(nSpecies*sizeof(long int));
	for(int s=0;s<nSpecies;s++){
		nAlloc[s] = ceil((double)nAllocTotal[s]/size);
		if(nAlloc[s]*size!=nAllocTotal[s])
			msg(WARNING|ONCE,"increased number of allocated particles from %i to %i to get integer per computing node",
				nAllocTotal[s],nAlloc[s]*size);
	}

	// Determine iStart and iStop
	long int *iStart = malloc((nSpecies+1)*sizeof(long int));
	long int *iStop  = malloc(nSpecies*sizeof(long int));

	iStart[0] = 0;
	for(int s=1;s<nSpecies+1;s++) iStart[s]=iStart[s-1]+nAlloc[s-1];
	for(int s=0;s<nSpecies;s++) iStop[s]=iStart[s]-1; // No particles yet, set stop before iStart.

	free(nAlloc);
	free(nAllocTotal);

	// Store in struct
	Population *pop = malloc(sizeof(Population));
	pop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->q = iniGetDoubleArr(ini,"population:q",&nSpecies);
	pop->m = iniGetDoubleArr(ini,"population:m",&nSpecies);
	pop->nSpecies = nSpecies;
	pop->nDims = nDims;
	pop->iStart = iStart;
	pop->iStop = iStop;
	pop->energy = NULL;	// Assume unused for now

	return pop;

}

void freePopulation(Population *pop){

	free(pop->pos);
	free(pop->vel);
	free(pop->q);
	free(pop->m);
	free(pop->energy);
	free(pop->iStart);
	free(pop->iStop);
	free(pop);

}

void posUniform(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const gsl_rng *rng){

	// Read from ini
	int nSpecies, nDims;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",&nSpecies);
	int *nTGPoints = iniGetIntArr(ini,"grid:nTGPoints",&nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*nTGPoints[d]-1;

	for(int s=0;s<nSpecies;s++){

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart-1;
		double *pos = &pop->pos[iStart*nDims];

		// Iterate through all particles to be generated
		// Same seed on all MPI nodes ensure same particles are generated everywhere.
		for(long int i=0;i<nParticles[s];i++){	// Generate nParticles of this specie

			// Generate position for particle i
			for(int d=0;d<nDims;d++) pos[d] = L[d]*gsl_rng_uniform_pos(rng);

			// Count the number of dimensions where the particle resides in the range of this node
			int correctRange = 0;
			for(int d=0;d<nDims;d++) correctRange += (subdomain[d] == (int)(posToSubdomain[d]*pos[d]));
//			msg(STATUS,"pos of particle %i: %f,%f,%f",i,pos[0],pos[1],pos[2]);

			// Iterate only if particle resides in this sub-domain.
			if(correctRange==nDims){
				pos += nDims;
				iStop++;
			}

		}

		if(iStop>=pop->iStart[s+1]){
			int allocated = pop->iStart[s+1]-iStart;
			int generated = iStop-iStart+1;
			msg(ERROR,"allocated only %i particles of specie %i per node but %i generated",allocated,s,generated);
		}

		pop->iStop[s]=iStop;

	}

	// Transform to local reference frame
	int *offset = mpiInfo->offset;
	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<=iStop;i++){

			double *pos = &pop->pos[i*nDims];
			for(int d=0;d<nDims;d++) pos[d] -= offset[d];
		}
	}

	free(L);
	free(nParticles);
	free(nTGPoints);

}

void posDebug(const dictionary *ini, Population *pop){

	int nSpecies;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",&nSpecies);

	int mpiRank, mpiSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	for(int s=0;s<nSpecies;s++){
		nParticles[s] /= mpiSize;
	}

	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){
		long int iStart = pop->iStart[s];
		pop->iStop[s] = iStart + nParticles[s] - 1;
		double *pos = &pop->pos[iStart*nDims];

		for(long int i=0;i<nParticles[s];i++){
			for(int d=0;d<nDims;d++){
				pos[i*nDims+d] = 1000*mpiRank + i + (double)d/10 + (double)s/100;
			}
		}
	}

	free(nParticles);

}

void velMaxwell(const dictionary *ini, Population *pop, const gsl_rng *rng){

	iniAssertEqualNElements(ini,3,"population:temperature","population:drift","population:nParticles");

	int nSpecies;
	double *temp = iniGetDoubleArr(ini,"population:temperature",&nSpecies);
	double *velDrift = iniGetDoubleArr(ini,"population:drift",&nSpecies);

	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = sqrt(temp[s]/temp[0]);

		for(long int i=iStart;i<=iStop;i++){

			double *vel = &pop->vel[i*nDims];
			for(int d=0;d<nDims;d++){
				vel[d] = velDrift[s] + gsl_ran_gaussian_ziggurat(rng,velTh);
			}
		}
	}
	free(temp);
	free(velDrift);
}

void createPopulationH5(const dictionary *ini, Population *pop){

	// Get MPI rank and size
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	// Determine filename
	char *fName = iniparser_getstring((dictionary *)ini,"files:output","");	// don't free
	char *fTotName = strAllocCat(fName,".pop.h5");

	// Create file with MPI-I/O access
	hid_t pList = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(pList,MPI_COMM_WORLD,MPI_INFO_NULL);
	hid_t file = createH5File(fTotName,H5P_DEFAULT,pList);
	H5Pclose(pList);

	free(fTotName);

	// Create groups
	hid_t group;
	group = H5Gcreate(file,"/pos",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Gclose(group);
	group = H5Gcreate(file,"/vel",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Gclose(group);

	char name[18];	// int is max 5 digits + "/pos/specie " + '\0'

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
}

void writePopulationH5(Population *pop, double posN, double velN){

	int mpiRank, mpiSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	int nSpecies = pop->nSpecies;

	const int arrSize = 2;
	hsize_t fileDims[arrSize];		// Dimensions of array in file
	hsize_t memDims[arrSize];		// Dimensions of array in memory
	hsize_t offset[arrSize];		// Where in file to store memory
	fileDims[1] = pop->nDims;
	memDims[1] = pop->nDims;
	offset[1] = 0;

	long int *nParticlesAllSubdomains = malloc((mpiSize+1)*sizeof(long int));
	nParticlesAllSubdomains[0] = 0;

	for(int s=0;s<nSpecies;s++){

		/*
		 * DETERMINE NUMBER OF PARTICLES TO STORE AND AT WHICH OFFSET
		 */
		long int nParticlesThisSubdomain = pop->iStop[s] - pop->iStart[s] + 1;
		MPI_Allgather(&nParticlesThisSubdomain,1,MPI_LONG,&nParticlesAllSubdomains[1],1,MPI_LONG,MPI_COMM_WORLD);

		// Cumulative sum
		for(int r=1;r<=mpiSize;r++){
			nParticlesAllSubdomains[r] += nParticlesAllSubdomains[r-1];
		}

		fileDims[0] = nParticlesAllSubdomains[mpiSize];
		memDims[0] = nParticlesThisSubdomain;
		offset[0] = nParticlesAllSubdomains[mpiRank];

		hid_t memSpace = H5Screate_simple(arrSize,memDims,NULL);
		hid_t fileSpace = H5Screate_simple(arrSize,fileDims,NULL);

		H5Sselect_hyperslab(fileSpace,H5S_SELECT_SET,offset,NULL,memDims,NULL);

		/*
		 * STORE DATA COLLECTIVELY
		 */
		hid_t pList = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);


		char name[64];
		hid_t dataset;

		sprintf(name,"/pos/specie %i/n=%.1f",s,posN);
		dataset = H5Dcreate(pop->h5,name,H5T_IEEE_F64LE,fileSpace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, &pop->pos[pop->iStart[s]*pop->nDims]);
		H5Dclose(dataset);

		sprintf(name,"/vel/specie %i/n=%.1f",s,velN);
		dataset = H5Dcreate(pop->h5,name,H5T_IEEE_F64LE,fileSpace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
		H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, &pop->vel[pop->iStart[s]*pop->nDims]);
		H5Dclose(dataset);

		H5Sclose(fileSpace);
		H5Sclose(memSpace);
		H5Pclose(pList);

	}
	free(nParticlesAllSubdomains);
}
