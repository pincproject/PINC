/**
 * @file		population.c
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
 * @brief Computes specie-specific re-normalization factors for the E-field
 * @param 	q			Charge [elementary charges] as given in input file
 * @param	m			Mass [electron masses] as given in input file
 * @param	nSpecies	Number of species (elements in q and m)
 * @return 				Pointer to allocated array of factors
 *
 * The E-field is normally normalized using the charge and mass of specie 0,
 * such that it becomes the discretized equivalent of acceleration of that
 * specie. Henceforth, adding the electric field to the velocity without any
 * form of multiplication readily yields the accelerated velocity. Multiplying
 * the E-field by renormE[0], where renormE is the return value of this
 * function, re-normalizes the E-field such that it represents the discretized
 * acceleration for specie 1. Subsequent multiplication by renormE[1] (such that
 * it has been multiplied by renormE[0]*renormE[1]) re-normalizes the E-field
 * for specie 2 and so on. The last multiplicative factor, renorm[nSpecies-1],
 * renormalizes it back to the original value, as normalized with respect to
 * specie 0.
 *
 * E is normalized as follows:
 * \f[
 *	\bar{E_x} = K_1\frac{(\Delta t)^2}{\Delta x}\frac{q_0}{m_0}E_x
 * \f]
 * and similarly for the other components. \f$K_1\f$ is a factor which may be
 * set unequal to 1 to speed up certain calculations.
 *
 * Remember to free the return value using free().
 */
double *computeRenormE(const double *q, const double *m, int nSpecies);

/**
 * @brief Computes specie-specific re-normalization factors for the charge density
 * @param 	q			Charge [elementary charges] as given in input file
 * @param	m			Mass [electron masses] as given in input file
 * @param	nSpecies	Number of species (elements in q and m)
 * @param	cellVolume	The volume of a cell in the mesh (product of 'stepSize')
 * @param	timeStep	The time step
 * @param	factor		Multiplicative factor used in normalizing rho
 * @return 				Pointer to allocated array of factors
 *
 * When distributing charges onto the grid to compute the charge density rho,
 * rather than having to multiply each contribution by the charge of the
 * particle, rho can temporarily have a specie-specific normalization allowing
 * the contribution simply to be added without having to multiply the weighting
 * coefficient by any charge.
 *
 * Starting with a grid value rho equal to zero this can be assumed to be
 * normalized for specie 0 and all the contributions for that specie can be
 * added. Multiplying the whole grid by renormRho[0] the re-normalizes it for
 * specie 1 such that it can be added. Doing this for all species, and
 * mulitplying by renormRho[nSpecies-1] at the end, results in a comletely
 * assembled charge density. The charge density in PINC is normalized
 * using the charge and mass of specie 0 as is other quantities.
 *
 * Rho is normalized as follows:
 * \f[
 *	\bar{\rho} = K_1K_2(\Delta t)^2\left(\frac{q_0}{m_0}\right)\rho
 * \f]
 * where \f$K_1K_2\f$ is a multiplicative factor which can be set un-equal to 1
 * to speed up certain calculations.
 *
 * Remember to free the return value using free().
 */
double *computeRenormRho(const double *q, const double *m, int nSpecies, double cellVolume, double timeStep, double factor);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

Population *pAlloc(const dictionary *ini){

	// Sanity check
	iniAssertEqualNElements(ini,4,"population:nParticles","population:nAlloc","population:q","population:m");

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);	// Presumes sanity check on nSubdomains by allocGrid()
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data
	int nSpecies;
	long int *nAllocTotal = iniGetLongIntArr(ini,"population:nAlloc",&nSpecies);	// This is for all computing nodes
	int nDims = iniGetNElements(ini,"grid:trueSize");
	if(nDims==0) msg(ERROR,"grid:trueSize not found");

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
	for(int s=0;s<nSpecies;s++) iStop[s]=iStart[s]; // No particles yet

	free(nAlloc);
	free(nAllocTotal);

	// Store in struct
	Population *pop = malloc(sizeof(Population));
	pop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->nSpecies = nSpecies;
	pop->nDims = nDims;
	pop->iStart = iStart;
	pop->iStop = iStop;
	pop->energy = NULL;	// Assume unused for now

	// Default normalization factors
	pSetSpecieNorm(ini,pop,1,1);

	return pop;

}

void pFree(Population *pop){

	free(pop->pos);
	free(pop->vel);
	free(pop->renormE);
	free(pop->renormRho);
	free(pop->energy);
	free(pop->iStart);
	free(pop->iStop);
	free(pop);

}

void pPosUniform(const dictionary *ini, Population *pop, const MpiInfo *mpiInfo, const gsl_rng *rng){

	// Read from ini
	int nSpecies, nDims;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",&nSpecies);
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",&nDims);

	// Read from mpiInfo
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	double *posToSubdomain = mpiInfo->posToSubdomain;

	// Compute normalized length of global reference frame
	int *L = malloc(nDims*sizeof(int));
	for(int d=0;d<nDims;d++) L[d] = nSubdomains[d]*trueSize[d]-1;

	for(int s=0;s<nSpecies;s++){

		// Start on first particle of this specie
		long int iStart = pop->iStart[s];
		long int iStop = iStart;
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

		if(iStop>pop->iStart[s+1]){
			int allocated = pop->iStart[s+1]-iStart;
			int generated = iStop-iStart;
			msg(ERROR,"allocated only %i particles of specie %i per node but %i generated",allocated,s,generated);
		}

		pop->iStop[s]=iStop;

	}

	pToLocalFrame(pop,mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);

}

void pPosDebug(const dictionary *ini, Population *pop){

	int nSpecies;
	long int *nParticles = iniGetLongIntArr(ini,"population:nParticles",&nSpecies);

	int mpiRank, mpiSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpiSize);

	for(int s=0;s<nSpecies;s++){
		nParticles[s] /= mpiSize;
	}

	int nDims = pop->nDims;	long int *nMigrantsResult = malloc(81*sizeof(*nMigrantsResult));
	alSet(nMigrantsResult,81,	1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,
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

void pVelMaxwell(const dictionary *ini, Population *pop, const gsl_rng *rng){

	iniAssertEqualNElements(ini,3,"population:temperature","population:drift","population:nParticles");

	int nSpecies;
	double *temp = iniGetDoubleArr(ini,"population:temperature",&nSpecies);
	double *velDrift = iniGetDoubleArr(ini,"population:drift",&nSpecies);

	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = sqrt(temp[s]/temp[0]);

		for(long int i=iStart;i<iStop;i++){

			double *vel = &pop->vel[i*nDims];
			for(int d=0;d<nDims;d++){
				vel[d] = velDrift[s] + gsl_ran_gaussian_ziggurat(rng,velTh);
			}
		}
	}
	free(temp);
	free(velDrift);
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

void pNew(Population *pop, int s, const double *pos, const double *vel){

	int nDims = pop->nDims;
	long int *iStart = pop->iStart;
	long int *iStop = pop->iStop;	// New particle added here

	if(iStop[s]>=iStart[s+1])
		msg(WARNING,"Not enough allocated memory to add new particle to specie %i. New particle ignored.",s);
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

	for(int d=0;d<nDims;d++){
		pos[d] = pop->pos[p+d];
		vel[d] = pop->vel[p+d];
		pop->pos[p+d] = pop->pos[pLast+d];
		pop->vel[p+d] = pop->vel[pLast+d];
	}

	pop->iStop[s]--;

}

void pCreateH5(const dictionary *ini, Population *pop, const char *fName){

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

	pToLocalFrame(pop,mpiInfo);
}

void pCloseH5(Population *pop){
	H5Fclose(pop->h5);
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

double *computeRenormE(const double *q, const double *m, int nSpecies){

	double *renormE = malloc(nSpecies*sizeof(*renormE));

	for(int s=0;s<nSpecies-1;s++){
		renormE[s] = (q[s+1]/m[s+1])/(q[s]/m[s]);
	}
	renormE[nSpecies-1] = (q[0]/m[0])/(q[nSpecies-1]/m[nSpecies-1]);

	return renormE;
}

double *computeRenormRho(const double *q, const double *m, int nSpecies, double cellVolume, double timeStep, double factor){

	double *qBar = malloc(nSpecies*sizeof(*qBar));
	double *renormRho = malloc(nSpecies*sizeof(*renormRho));

	for(int s=0;s<nSpecies;s++){
		qBar[s] = factor*(pow(timeStep,2)/cellVolume)*(q[0]/m[0])*q[s];
	}

	for(int s=0;s<nSpecies-1;s++){
		renormRho[s] = qBar[s]/qBar[s+1];
	}
	renormRho[nSpecies-1] = qBar[nSpecies-1];

	free(qBar);
	return renormRho;

}

void pSetSpecieNorm(const dictionary *ini, Population *pop, double timeStepMul, double factor){

	int nSpecies, nDims;
	double *q = iniGetDoubleArr(ini,"population:q",&nSpecies);
	double *m = iniGetDoubleArr(ini,"population:m",&nSpecies);
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",&nDims);
	double timeStep = iniparser_getdouble((dictionary *)ini,"time:timeStep",0.0);
	double cellVolume = adProd(stepSize,nDims);

	timeStep *= timeStepMul;

	pop->renormE = computeRenormE(q,m,nSpecies);
	pop->renormRho = computeRenormRho(q,m,nSpecies,cellVolume,timeStep,factor);

	free(q);
	free(m);
	free(stepSize);

}
