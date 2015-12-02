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

/******************************************************************************
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/

 /**
  * @brief Transforms particles to local reference frame
  * @param	pop			Population of particles
  * @param	mpiInfo		MPI information about the reference frames
  * @return	void
  * @see toGlobalFrame()
  */
void toLocalFrame(Population *pop, const MpiInfo *mpiInfo);

/**
 * @brief Transforms particles to global reference frame
 * @param	pop			Population of particles
 * @param	mpiInfo		MPI information about the reference frames
 * @return	void
 * @see toLocalFrame()
 *
 * For parallelization by means of configuration space decomposition, the
 * the particles' positions are usually specified with respect to a local
 * reference frame to that subdomain in order to ease computation. Some
 * operations may require the position in global reference frame (e.g. when
 * storing to file) for which purpose it can be transformed using this function.
 */
void toGlobalFrame(Population *pop, const MpiInfo *mpiInfo);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

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

	toLocalFrame(pop,mpiInfo);

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

void createPopulationH5(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const char *fName){

	/*
	 * CREATE FILE
	 */

	hid_t file = createH5File(ini,fName,"pop");

	/*
	 * CREATE GROUPS
	 */
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

	/*
	 * CREATE ATTRIBUTES
	 */

	int nDims;
	double *dr = iniGetDoubleArr(ini,"grid:dr",&nDims);
	double *attrData = malloc(nDims*sizeof(*attrData));
	hsize_t attrSize;
    hid_t attrSpace;
    hid_t attribute;

	attrSize = (hsize_t)nDims;
	attrSpace = H5Screate_simple(1,&attrSize,NULL);

	attribute = H5Acreate(file, "Position denormalization factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, dr);
    H5Aclose(attribute);

	double debye = iniparser_getdouble((dictionary *)ini,"grid:debye",0);
	for(int d=0;d<nDims;d++) attrData[d]=debye;
	attribute = H5Acreate(file, "Position dimensionalizing factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, attrData);
    H5Aclose(attribute);


	double dt = iniparser_getdouble((dictionary *)ini,"time:dt",0);
	for(int d=0;d<nDims;d++) attrData[d]=dr[d]/dt;
	attribute = H5Acreate(file, "Velocity denormalization factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, attrData);
	H5Aclose(attribute);

	double *T = iniGetDoubleArr(ini,"population:temperature",&nDims);
	double *m = iniGetDoubleArr(ini,"population:m",&nDims);
	double vThermal = sqrt(BOLTZMANN*T[0]/(m[0]*ELECTRON_MASS));
	for(int d=0;d<nDims;d++) attrData[d] = vThermal;
	attribute = H5Acreate(file, "Velocity dimensionalizing factor", H5T_IEEE_F64LE, attrSpace, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attribute, H5T_NATIVE_DOUBLE, attrData);
	H5Aclose(attribute);

	free(T);
	free(m);
	free(dr);
	free(attrData);

	H5Sclose(attrSpace);

}

void writePopulationH5(Population *pop, const MpiInfo *mpiInfo, double posN, double velN){

	int mpiRank = mpiInfo->mpiRank;
	int mpiSize = mpiInfo->mpiSize;
	int nSpecies = pop->nSpecies;

	toGlobalFrame(pop,mpiInfo);

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

		long int nParticles = pop->iStop[s] - pop->iStart[s] + 1;
		MPI_Allgather(&nParticles,1,MPI_LONG,&offsetAllSubdomains[1],1,MPI_LONG,MPI_COMM_WORLD);

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

			H5Pclose(pList);
			H5Sclose(fileSpace);
			H5Sclose(memSpace);

		} else {
			msg(WARNING|ONCE,"No particles to store in .h5-file");
		}
	}
 	free(offsetAllSubdomains);

	toLocalFrame(pop,mpiInfo);
}

void closePopulationH5(Population *pop){
	H5Fclose(pop->h5);
}

/******************************************************************************
 * DEFINING LOCAL FUNCTIONS
 *****************************************************************************/

void toLocalFrame(Population *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<=iStop;i++){

			double *pos = &pop->pos[i*nDims];
			for(int d=0;d<nDims;d++) pos[d] -= offset[d];
		}
	}
}

void toGlobalFrame(Population *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<=iStop;i++){

			double *pos = &pop->pos[i*nDims];
			for(int d=0;d<nDims;d++) pos[d] += offset[d];
		}
	}
}
