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
 * DECLARING LOCAL FUNCTIONS
 *****************************************************************************/

/**
 * @brief	Sets normalization parameters in Population
 * @param	ini				Dictionary to input file
 * @param	pop[in,out]		Population
 *
 * Normalizes charge and mass and sets specie-specific renormalization
 * parameters in Population.
 *
 */
static void pSetNormParams(const dictionary *ini, Population *pop);

/******************************************************************************
 * DEFINING GLOBAL FUNCTIONS
 *****************************************************************************/

Population *pAlloc(const dictionary *ini){

	// Get MPI info
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	// Load data
	int nSpecies = iniGetInt(ini,"population:nSpecies");
	int nDims = iniGetInt(ini,"grid:nDims");

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

	free(nAlloc);
	free(nAllocTotal);

	Population *pop = malloc(sizeof(Population));
	pop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	pop->nSpecies = nSpecies;
	pop->nDims = nDims;
	pop->iStart = iStart;
	pop->iStop = iStop;
	pop->kinEnergy = malloc((nSpecies+1)*sizeof(double));
	pop->potEnergy = malloc((nSpecies+1)*sizeof(double));

	pSetNormParams(ini,pop);

	return pop;

}

void pFree(Population *pop){

	free(pop->pos);
	free(pop->vel);
	free(pop->renormE);
	free(pop->renormRho);
	free(pop->kinEnergy);
	free(pop->potEnergy);
	free(pop->iStart);
	free(pop->iStop);
	free(pop->charge);
	free(pop->mass);
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
			for(int d=0;d<nDims;d++) pos[d] = L[d]*gsl_rng_uniform_pos(rng);

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
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",nDims);
	double *amplitude = iniGetDoubleArr(ini,"population:perturbAmplitude",nElements);
	double *mode = iniGetDoubleArr(ini,"population:perturbMode",nElements);

	for(int e = 0; e < nElements; e++) {
		if (nDims > 1)	amplitude[e] /= stepSize[e%(nDims-1)];
		else amplitude[e] /= stepSize[0];
	}
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
	double *velDrift = iniGetDoubleArr(ini,"population:drift",nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"population:thermalVelocity",nSpecies);
	double timeStep = iniGetDouble(ini,"time:timeStep");
	double stepSize = iniGetDouble(ini,"grid:stepSize");

	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = (timeStep/stepSize)*velThermal[s];

		for(long int i=iStart;i<iStop;i++){

			double *vel = &pop->vel[i*nDims];
			for(int d=0;d<nDims;d++){
				vel[d] = velDrift[s] + gsl_ran_gaussian_ziggurat(rng,velTh);
			}
		}
	}
	free(velDrift);
	free(velThermal);
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

	for(int d=0;d<nDims;d++){
		pos[d] = pop->pos[p+d];
		vel[d] = pop->vel[p+d];
		pop->pos[p+d] = pop->pos[pLast+d];
		pop->vel[p+d] = pop->vel[pLast+d];
	}

	pop->iStop[s]--;

}

void pOpenH5(const dictionary *ini, Population *pop, const char *fName){

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

	int nDims = pop->nDims;
	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",nDims);
	double timeStep = iniGetDouble(ini,"time:timeStep");
	double debye = 0; // TBD: Does not work

	double vThermal = 0; // TBD: Does not work

	double *attrData = malloc(nDims*sizeof(*attrData));

	setH5Attr(file, "Position denormalization factor", stepSize, nDims);

	adSetAll(attrData,nDims,debye);
	setH5Attr(file, "Position dimensionalizing factor", attrData, nDims);

	for(int d=0;d<nDims;d++) attrData[d]=stepSize[d]/timeStep;
	setH5Attr(file, "Velocity denormalization factor", attrData, nDims);

	adSetAll(attrData,nDims,vThermal);
	setH5Attr(file, "Velocity dimensionalizing factor", attrData, nDims);

	free(stepSize);
	free(attrData);

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

	char name[32];
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

void pWriteEnergy(hid_t xy, Population *pop, double x){

	char name[32];
	int nSpecies = pop->nSpecies;

	sprintf(name,"/energy/potential/total");
	xyWrite(xy,name,x,pop->potEnergy[nSpecies],MPI_SUM);

	sprintf(name,"/energy/kinetic/total");
	xyWrite(xy,name,x,pop->kinEnergy[nSpecies],MPI_SUM);

	for(int s=0; s<nSpecies; s++){

		sprintf(name,"/energy/potential/specie %i",s);
		xyWrite(xy,name,x,pop->potEnergy[s],MPI_SUM);

		sprintf(name,"/energy/kinetic/specie %i",s);
		xyWrite(xy,name,x,pop->kinEnergy[s],MPI_SUM);
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

static void pSetNormParams(const dictionary *ini, Population *pop){

	int nSpecies = pop->nSpecies;
	int nDims = pop->nDims;
	double *charge = iniGetDoubleArr(ini,"population:charge",nSpecies);
	double *mass = iniGetDoubleArr(ini,"population:mass",nSpecies);
	double *multiplicity = iniGetDoubleArr(ini,"population:multiplicity",nSpecies);

	double *stepSize = iniGetDoubleArr(ini,"grid:stepSize",nDims);
	double timeStep = iniGetDouble(ini,"time:timeStep");
	double cellVolume = adProd(stepSize,nDims);

	// Multiply charge and mass by multiplicity
	adMul(charge,multiplicity,charge,nSpecies);
	adMul(mass,multiplicity,mass,nSpecies);

	/*
	 * Normalizing charge and mass (used in energy computations)
	 */
	double *chargeBar = malloc(nSpecies*sizeof(*chargeBar));
	double *massBar = malloc(nSpecies*sizeof(*massBar));
	for(int s=0;s<nSpecies;s++){
		chargeBar[s] = (pow(timeStep,2)/cellVolume)*(charge[0]/mass[0])*charge[s];
		// massBar[s] 	= (pow(timeStep,4)/cellVolume)*pow(charge[0]/mass[0],2)*mass[s];
		massBar[s]	 = pow(timeStep,2)*mass[s];
	}
	pop->charge=chargeBar;
	pop->mass=massBar;

	/*
	 * Computing renormalization factors (used in update equations)
	 */
	double *renormE = malloc(nSpecies*sizeof(*renormE));
	double *renormRho = malloc(nSpecies*sizeof(*renormRho));
	for(int s=0;s<nSpecies-1;s++){
		renormE[s] = (charge[s+1]/mass[s+1])/(charge[s]/mass[s]);
		renormRho[s] = chargeBar[s]/chargeBar[s+1];
		// renormRho[s] = charge[s]/charge[s+1];
	}
	renormE[nSpecies-1] = (charge[0]/mass[0])/(charge[nSpecies-1]/mass[nSpecies-1]);
	// renormRho[nSpecies-1] = chargeBar[nSpecies-1];
	renormRho[nSpecies-1] = pow(timeStep,2)*(charge[0]/mass[0])*charge[nSpecies-1];

	pop->renormE = renormE;
	pop->renormRho = renormRho;

	free(charge);
	free(mass);
	free(multiplicity);
	free(stepSize);

}
