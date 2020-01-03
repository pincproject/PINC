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
#include "pusher.h"



/******************************************************************************
* DEFINING LOCAL FUNCTIONS
*****************************************************************************/

void gDirichletRho(Grid *grid, const int boundary,  const  MpiInfo *mpiInfo){

	//msg(STATUS, "Hello from Dirichlet");

	//Load data
	int rank = grid->rank;
	int *size = grid->size;
	double *bndSlice = grid->bndSlice;

	//Compute dimensions and size of slice
	int d = boundary%rank;
	int offset = 1 + (boundary>rank)*(size[d]-3);

	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}
	//msg(STATUS,"offset. eg. index to set slice in perp. direction %i",offset);
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset); //edge before halo
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset - 1 + (boundary>rank)*2); //halo
	setSlice(&bndSlice[boundary*nSliceMax], grid, d, offset + 1 - (boundary>rank)*2); //edge before edge

	//adPrint(&bndSlice[(boundary)*nSliceMax], nSliceMax);

	return;

}

void nuGBnd(Grid *grid, const MpiInfo *mpiInfo){

	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;

	//If periodic neutralize phi

	bool periodic =  true;
	for(int d = 1; d < rank; d++){
		//msg(STATUS,"d = %i",d);
		if(bnd[d] != PERIODIC){
			//msg(STATUS,"bnd[d] != PERIODIC, d = %i",d);
			periodic = false;
		}
	}
	for(int d = rank+1; d < 2*rank; d++){
		//msg(STATUS,"d = %i",d);
		if(bnd[d] != PERIODIC){
			//msg(STATUS,"bnd[d] != PERIODIC, d = %i",d);
			periodic = false;
		}
	}
	if(periodic == true){
		//printf("PERIODIC cond, rank %i\n",mpiInfo->mpiRank);
		//gPeriodic(grid, mpiInfo);

	}
	//gHaloOp(setSlice, grid,mpiInfo,TOHALO);
	//Lower edge
	for(int d = 1; d < rank; d++){
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				//msg(STATUS,"bnd[d] = DIRICHLET, giving d = %i, rank = %i",d,rank);
				gDirichletRho(grid, d, mpiInfo);
				//msg(STATUS,"bnd[d] = DIRICHLET, d = %i",d);
			}
			else if(bnd[d] == NEUMANN){
				//msg(STATUS,"bnd[d] = NEUMANN, d = %i",d);
				gNeumann(grid, d, mpiInfo);
			}
		}
	}

	//Higher edge
	for(int d = rank+1; d < 2*rank; d++){
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET) gDirichletRho(grid, d, mpiInfo);
			if(bnd[d] == NEUMANN)	gNeumann(grid, d, mpiInfo);
		}
	}
	//printf("after boundary cond, rank %i\n",mpiInfo->mpiRank);
	//msg(STATUS,"phi size = %i",grid->sizeProd[4]);
	//if (mpiInfo->mpiRank == 7){
	//adPrint(grid->val,grid->sizeProd[4] );
	//}
	//MPI_Barrier(MPI_COMM_WORLD);
	//exit(1);
	return;
}


void neNewparticle(NeutralPopulation *pop, int s, const double *pos, const double *vel){

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

void neInjectParticles(int slicePos,int dim,int multiplyDens,const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo){

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	bndType *bnd = pop->bnd;
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	double *velDrift = iniGetDoubleArr(ini,"collisions:neutralDrift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"collisions:thermalVelocityNeutrals",nSpecies);
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

	nePToGlobalFrame(pop,mpiInfo);

	for(int s=0;s<nSpecies;s++){

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		double velTh = velThermal[s];

		int d=dim;

		long int globalSizeProd = (L[0]*L[1]*L[2]); //TODO: make nDimensional

		long int sliceSize = (globalSizeProd/(L[d]));
		long int newParticles = multiplyDens*( (sliceSize*nParticles[s])/globalSizeProd ); // slice * particles per cell
		//printf("generating %li particles for specie %i, velTh = %f \n \n",newParticles,s,velTh);
		//for(long int i=iStart;i<iStop;i++){
		//for(long int i=0;i<newParticles;i++){


		for(long int i=0;i<newParticles;i++){
			//generate velocity for particle
			for(int dd=0;dd<nDims;dd++){
				vel[dd] = velDrift[(s*nDims)+dd] + gsl_ran_gaussian_ziggurat(rng,velTh);
				//printf("vel = %f, drift = %f \n",vel[dd],velDrift[(s*nDims)+dd] );
			}

			// Generate position for particle
			for(int dd=0;dd<nDims;dd++){
				pos[dd] = (L[dd])*gsl_rng_uniform_pos(rng)-nGhostLayers[d];
			}
			pos[d] = slicePos+(gsl_rng_uniform_pos(rng)); //in lower ghost
			//printf("nGhostLayers[d+1] = %i\n",nGhostLayers[d]);
			int correctRange = 0;
			for(int dd=0;dd<nDims;dd++){
				//printf("posToSubdomain[dd] = %f, dd = %i \n",posToSubdomain[dd],dd);
				correctRange += (subdomain[dd] == (int)(posToSubdomain[dd]*pos[dd]));
			}
			// Add only if particle resides in this sub-domain.
			if(correctRange==nDims){
				neNewparticle(pop,s,pos,vel);
			}
		}
	}
	nePToLocalFrame(pop,mpiInfo);

	free(velDrift);
	free(velThermal);
	free(trueSize);
	free(nParticles);
	free(nGhostLayers);
	free(L);

	return;
}

void neScatterParticle(NeutralPopulation *pop, int s, long int p, double *pos, double *vel){

    // neutral object coll
	int nDims = pop->nDims;
	long int pLast = (pop->iStop[s]-1)*nDims;
	double *newvel = vel;
	//printf("Cutting particle p = %li \n",p);
	//Decide in  whitch dimension we crossed boundary
	double oldPos[nDims];
	double oldVel[nDims];
	int traversed[nDims];
	double delta[nDims];



	for(int d=0;d<nDims;d++){
		oldVel[d] = vel[d]; // using oldVel temporarily.
		oldPos[d] = pos[d]-oldVel[d];
		//printf("\n");
		//printf("pos[d] = %f, oldPos[d]=%f \n",pos[d],oldPos[d] );
		//printf("traversed[d] = %i \n",((int)(pos[d])-(int)(oldPos[d])));
		if( (int)(pos[d])-(int)(oldPos[d])!=0){
			traversed[d]=(int)(pos[d])-(int)(oldPos[d]);
			delta[d] = pos[d]-(int)pos[d];
			//printf("delta[%i]= %f \n",d,delta[d]);
			//printf("new pos will be = %f \n",((int)pos[d]-delta[d]) );
			pos[d] = ((int)pos[d]-delta[d]);
			vel[d] = - oldVel[d];
			//printf("collided in dim %i \n",d);
			if(traversed[d]==-1){ //hotfix
				pos[d] += 2;
			}
		}
	}
	// for(int d=0;d<nDims;d++){
	// 	//pos[d] = -pop->pos[p+d];
	// 	//printf("vel[%i] = %f \n",d,newvel[d]);
	// 	newvel[d] = -pop->vel[d];
	// 	pos[d] += newvel[d]; // move back
	// 	//printf("vel[%i] = %f \n \n",d,newvel[d]);
	// 	//printf("pos[%i] = %f \n",d,pos[d]);
	// 	//pop->pos[p+d] = pop->pos[pLast+d];
	// 	//pop->vel[p+d] = pop->vel[pLast+d];
	// }


	//pop->iStop[s]--;
	//exit(0);


}



void neCutParticle(NeutralPopulation *pop, int s, long int p, double *pos, double *vel){

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


void nePToLocalFrame(NeutralPopulation *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpeciesNeutral;
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

void nePToGlobalFrame(NeutralPopulation *pop, const MpiInfo *mpiInfo){

	int *offset = mpiInfo->offset;
	int nSpecies = pop->nSpeciesNeutral;
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



/******************************************************************************
* DEFINING GLOBAL FUNCTIONS
*****************************************************************************/

NeutralPopulation *pNeutralAlloc(const dictionary *ini,const MpiInfo *mpiInfo){

	// Get MPI info
	int size;//, rank;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int mpiRank = mpiInfo->mpiRank;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;

	// Load data
	int nSpecies= iniGetInt(ini,"collisions:nSpeciesNeutral");
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
	//printf(" alloc vals = %li \n",((long int)nDims*iStart[nSpecies]));
	double *mass = iniGetDoubleArr(ini, "collisions:neutralMass", nSpecies);
	double *density = iniGetDoubleArr(ini, "collisions:numberDensityNeutrals", nSpecies);
	double rho0 = density[0];//*mass[0];

    double NvelThermal = iniGetDouble(ini,"collisions:thermalVelocityNeutrals");
	double stiffnessC =  (NvelThermal); //((7.)*(NvelThermal*NvelThermal))/(mass[0]*rho0); // Monaghan (1994)
	//msg(STATUS,"mass[0]= %f, rho0 = %f, density  = %f, NvelThermal = %f, stiffnessC = %f \n",mass[0],rho0,density [0],NvelThermal,stiffnessC);
	//exit(0);
	NeutralPopulation *neutralPop = malloc(sizeof(*neutralPop));
	neutralPop->pos = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	neutralPop->vel = malloc((long int)nDims*iStart[nSpecies]*sizeof(double));
	neutralPop->nSpeciesNeutral = nSpecies;
	neutralPop->nDims = nDims;
	neutralPop->iStart = iStart;
	neutralPop->iStop = iStop;
	//neutralPop->kinEnergy = malloc((nSpecies+1)*sizeof(double));
	//neutralPop->potEnergy = malloc((nSpecies+1)*sizeof(double));
	neutralPop->mass = mass;//iniGetDoubleArr(ini,"collisions:neutralMass",nSpecies);
	neutralPop->bnd = bnd;
	neutralPop->stiffnessConstant = stiffnessC; //0.00001
	neutralPop->rho0 = rho0;
	//printf("neutralPop->mass = %f \n",neutralPop->mass[0] );

	free(nAlloc);
	free(nAllocTotal);
	free(density);
	free(boundaries);
	return neutralPop;

}

void pNeutralFree(NeutralPopulation *pop){

	free(pop->pos);
	free(pop->vel);
	//free(pop->kinEnergy);
	//free(pop->potEnergy);
	free(pop->iStart);
	free(pop->iStop);
	free(pop->mass);
	free(pop->bnd);
	free(pop);


}


void nePosLattice(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo){

	// Read from ini
	int nDims = pop->nDims;
	int nSpecies = pop->nSpeciesNeutral;
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

	nePToLocalFrame(pop,mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);

}

void nePosUniform(const dictionary *ini, NeutralPopulation *pop, const MpiInfo *mpiInfo, const gsl_rng *rng){

	// Read from ini
	int nSpecies = pop->nSpeciesNeutral;
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

	nePToLocalFrame(pop,mpiInfo);

	free(L);
	free(nParticles);
	free(trueSize);

}





void neVelMaxwell(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng){

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	double *velDrift = iniGetDoubleArr(ini,"collisions:neutralDrift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"collisions:thermalVelocityNeutrals",nSpecies);
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

void nePurgeGhost(NeutralPopulation *pop, const Grid *grid){

	//this will delete all particles reciding on Dirichlet boundarys
	int *size = grid->size;
	int *nGhostLayers = grid->nGhostLayers;
	int rank = grid->rank;
	bndType *bnd = grid->bnd;
	double *pos = pop->pos;
	double *vel = pop->vel;
	bool cut = false;
	int nSpecies = pop->nSpeciesNeutral;
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
				if( (pos[i*nDims+d]>size[d+1]-nGhostLayers[d+1]-1 && (bnd[d+rank+1]==DIRICHLET || bnd[d+rank+1]==NEUMANN)) ){
					//msg(STATUS,"CUT: pos[i*nDims+d] = %f, bnd[d+1]=%d \n",pos[i*nDims+d],bnd[d+1]);
					cut = true;
				}
				if( (pos[i*nDims+d]<nGhostLayers[d+rank+1] && (bnd[d+1]==DIRICHLET || bnd[d+1]==NEUMANN)) ){
					//msg(STATUS,"CUT: pos[i*nDims+d] = %f, bnd[d+rank+1]=%d \n",pos[i*nDims+d],bnd[d+rank+1]);
					cut = true;
				}
			}
			if (cut == true){
				neCutParticle(pop, s, i*nDims, pos, vel);
				cut = false;
				iStop--;
				i--;
				//printf("iStop= %i \n \n",iStop);
			}
		}
	}
}

void neFillGhost(const dictionary *ini, NeutralPopulation *pop, const gsl_rng *rng, const MpiInfo *mpiInfo){

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;
	bndType *bnd = pop->bnd;
	int *trueSize = iniGetIntArr(ini,"grid:trueSize",nDims);
	double *velDrift = iniGetDoubleArr(ini,"collisions:neutralDrift",nDims*nSpecies);
	double *velThermal = iniGetDoubleArr(ini,"collisions:thermalVelocityNeutrals",nSpecies);
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

	nePToGlobalFrame(pop,mpiInfo);

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
			if(bnd[d+1]==DIRICHLET || bnd[d+1]==NEUMANN){
				for(long int i=0;i<newParticles;i++){
					//generate velocity for particle
					for(int dd=0;dd<nDims;dd++){
						vel[dd] = velDrift[(s*nDims)+dd] + gsl_ran_gaussian_ziggurat(rng,velTh);
						//printf("vel = %f, drift = %f \n",vel[dd],velDrift[(s*nDims)+dd] );

					}
					//printf(" \n");

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
						neNewparticle(pop,s,pos,vel);
					}
				}
			}

			//Upper ghost slice
			if(bnd[d+rank+1]==DIRICHLET || bnd[d+rank+1]==NEUMANN){
				for(long int i=0;i<newParticles;i++){
					//generate velocity for particle

					for(int dd=0;dd<nDims;dd++){
						vel[dd] = velDrift[(s*nDims)+dd] + gsl_ran_gaussian_ziggurat(rng,velTh);
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
						neNewparticle(pop,s,pos,vel);

					}
				}



			}
		}
	}
	nePToLocalFrame(pop,mpiInfo);

	free(velDrift);
	free(velThermal);
	free(trueSize);
	free(nParticles);
	free(nGhostLayers);
	free(L);

	return;
}



//#########################################
// Distributer
// ########################################



funPtr NeutralDistr3D1_set(dictionary *ini){
	//puSanity(ini,"puDistr3D1",3,1);
	return NeutralDistr3D1;
}
void NeutralDistr3D1(const NeutralPopulation *pop, Grid *rho){

	gZero(rho);
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;
	//adPrint(val,rho->sizeProd[2]);
	int nSpecies = pop->nSpeciesNeutral;

	for(int s=0;s<nSpecies;s++){

		//gMul(rho, 1.0/pop->mass[s]);

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[3*i];

			// Integer parts of position
			int j = (int) pos[0];
			int k = (int) pos[1];
			int l = (int) pos[2];

			// Decimal (cell-referenced) parts of position and their complement
			double x = pos[0]-j;
			double y = pos[1]-k;
			double z = pos[2]-l;
			double xcomp = 1-x;
			double ycomp = 1-y;
			double zcomp = 1-z;

			//printf("")
			// Index of neighbouring nodes
			long int p 		= j + k*sizeProd[2] + l*sizeProd[3];
			long int pj 	= p + 1; //sizeProd[1];
			long int pk 	= p + sizeProd[2];
			long int pjk 	= pk + 1; //sizeProd[1];
			long int pl 	= p + sizeProd[3];
			long int pjl 	= pl + 1; //sizeProd[1];
			long int pkl 	= pl + sizeProd[2];
			long int pjkl 	= pkl + 1; //sizeProd[1];


			// if(p>=sizeProd[4]){
			// 	msg(ERROR,"Particle %i at (%f,%f,%f) out-of-bounds, tried to access node %li",i,pos[0],pos[1],pos[2],pjkl);
			// }
			//12294
			//printf("p = %li\n",sizeProd[4]);
			//printf("val[p] = %f\n",val[p]);
			//MPI_Barrier(MPI_COMM_WORLD);
			val[p] 		+= xcomp*ycomp*zcomp;
			val[pj]		+= x    *ycomp*zcomp;
			val[pk]		+= xcomp*y    *zcomp;
			val[pjk]	+= x    *y    *zcomp;
			val[pl]     += xcomp*ycomp*z    ;
			val[pjl]	+= x    *ycomp*z    ;
			val[pkl]	+= xcomp*y    *z    ;
			val[pjkl]	+= x    *y    *z    ;
			//if(val[p]>100. || val[pj]>100. || val[pk]>100. || val[pjk]>100. || val[pl]>100. || val[pjl]>100. || val[pkl]>100. || val[pjkl]>100.){
            	//printf("val = %f, %f, %f, %f, %f, %f, %f, %f \n",val[p],val[pj],val[pk],val[pjk],val[pl],val[pjl],val[pkl],val[pjkl]);
			//}
		}

		//gMul(rho, pop->mass[s]);
		//adPrint(val,rho->sizeProd[4]);
		//adPrint(rho->val,rho->sizeProd[4]);
		//exit(0);
	}


}


funPtr NeutralDistr3D1Vector_set(dictionary *ini){
	//puSanity(ini,"puDistr3D1",3,1);
	return NeutralDistr3D1Vector;
}
void NeutralDistr3D1Vector(const NeutralPopulation *pop, Grid *rho){

	gZero(rho);
	double *val = rho->val;
	long int *sizeProd = rho->sizeProd;
	//adPrint(val,rho->sizeProd[2]);
	int nSpecies = pop->nSpeciesNeutral;
	int nDims = pop->nDims;

	for(int s=0;s<nSpecies;s++){

		//gMul(rho, 1.0/pop->mass[s]);

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		//gMul(rho, 1.0/(iStart-iStop));

		for(int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[3*i];
			double *vel = &pop->vel[3*i];

			// Integer parts of position
			int j = (int) pos[0];
			int k = (int) pos[1];
			int l = (int) pos[2];

			// Decimal (cell-referenced) parts of position and their complement
			double x = pos[0]-j;
			double y = pos[1]-k;
			double z = pos[2]-l;
			double xcomp = 1-x;
			double ycomp = 1-y;
			double zcomp = 1-z;

			//printf("sizeProd[1] = %i \n",sizeProd[1]);
			// Index of neighbouring nodes
			long int p 		= j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
			long int pj 	= p + sizeProd[1];
			long int pk 	= p + sizeProd[2];
			long int pjk 	= pk + sizeProd[1];
			long int pl 	= p + sizeProd[3];
			long int pjl 	= pl + sizeProd[1];
			long int pkl 	= pl + sizeProd[2];
			long int pjkl 	= pkl + sizeProd[1];


			// if(p>=sizeProd[4]){
			// 	msg(ERROR,"Particle %i at (%f,%f,%f) out-of-bounds, tried to access node %li",i,pos[0],pos[1],pos[2],pjkl);
			// }
			//12294
			//printf("p = %li\n",sizeProd[4]);
			//printf("val[p] = %f\n",val[p]);
			//MPI_Barrier(MPI_COMM_WORLD);
			for (int d=0;d<nDims;d++){
				  double K = vel[d]/12.;
    			val[p+d] 		+= K*xcomp*ycomp*zcomp;
	    		val[pj+d]		+= K*x    *ycomp*zcomp;
		    	val[pk+d]		+= K*xcomp*y    *zcomp;
    			val[pjk+d]	+= K*x    *y    *zcomp;
		    	val[pl+d]   += K*xcomp*ycomp*z    ;
			    val[pjl+d]	+= K*x    *ycomp*z    ;
			    val[pkl+d]	+= K*xcomp*y    *z    ;
			    val[pjkl+d]	+= K*x    *y    *z    ;
					//printf("val[%i] = %f \n",p+d,val[p+d]);
			}
			//if(val[p]>100. || val[pj]>100. || val[pk]>100. || val[pjk]>100. || val[pl]>100. || val[pjl]>100. || val[pkl]>100. || val[pjkl]>100.){
            	//printf("val = %f, %f, %f, %f, %f, %f, %f, %f \n",val[p],val[pj],val[pk],val[pjk],val[pl],val[pjl],val[pkl],val[pjkl]);
			//}
		}

		//gMul(rho, (iStart-iStop));
		//gMul(rho, pop->mass[s]);
		//adPrint(val,rho->sizeProd[4]);
		//adPrint(rho->val,rho->sizeProd[4]);
		//exit(0);
	}


}


//#############################
// accelerator
//#############################

static inline void neInterp3D1(	double *result, const double *pos,
	const double *val, const long int *sizeProd){

		// Integer parts of position
		int j = (int) pos[0];
		int k = (int) pos[1];
		int l = (int) pos[2];

		// Decimal (cell-referenced) parts of position and their complement
		double x = pos[0]-j;
		double y = pos[1]-k;
		double z = pos[2]-l;
		double xcomp = 1-x;
		double ycomp = 1-y;
		double zcomp = 1-z;

		// Index of neighbouring nodes
		long int p 		= j*3 + k*sizeProd[2] + l*sizeProd[3];
		long int pj 	= p + 3; //sizeProd[1];
		long int pk 	= p + sizeProd[2];
		long int pjk 	= pk + 3; //sizeProd[1];
		long int pl 	= p + sizeProd[3];
		long int pjl 	= pl + 3; //sizeProd[1];
		long int pkl 	= pl + sizeProd[2];
		long int pjkl 	= pkl + 3; //sizeProd[1];

		// Linear interpolation
		for(int v=0;v<3;v++)
		result[v] =	zcomp*(	 ycomp*(xcomp*val[p   +v]+x*val[pj  +v])
		+y    *(xcomp*val[pk  +v]+x*val[pjk +v]) )
		+z    *( ycomp*(xcomp*val[pl  +v]+x*val[pjl +v])
		+y    *(xcomp*val[pkl +v]+x*val[pjkl+v]) );

}

static inline void neInterp3D1scalar(	double *result, const double *pos,
	const double *val, const long int *sizeProd){

		result[0] = 0;
		// Integer parts of position
		int j = (int) pos[0];
		int k = (int) pos[1];
		int l = (int) pos[2];

		// Decimal (cell-referenced) parts of position and their complement
		double x = pos[0]-j;
		double y = pos[1]-k;
		double z = pos[2]-l;
		double xcomp = 1-x;
		double ycomp = 1-y;
		double zcomp = 1-z;

		// Index of neighbouring nodes
		long int p 		= j*sizeProd[1] + k*sizeProd[2] + l*sizeProd[3];
		long int pj 	= p + sizeProd[1]; //sizeProd[1];
		long int pk 	= p + sizeProd[2];
		long int pjk 	= pk + sizeProd[1]; //sizeProd[1];
		long int pl 	= p + sizeProd[3];
		long int pjl 	= pl + sizeProd[1]; //sizeProd[1];
		long int pkl 	= pl + sizeProd[2];
		long int pjkl 	= pkl + sizeProd[1]; //sizeProd[1];

		// Linear interpolation
		//for(int v=0;v<3;v++)
		result[0] =	zcomp*(	 ycomp*(xcomp*val[p]+x*val[pj])
		+y    *(xcomp*val[pk]+x*val[pjk]) )
		+z    *( ycomp*(xcomp*val[pl]+x*val[pjl])
		+y    *(xcomp*val[pkl]+x*val[pjkl]) );

}


funPtr neAcc3D1_set(dictionary *ini){
	//puSanity(ini,"puAcc3D1",3,1);
	return neAcc3D1;
}
void neAcc3D1(NeutralPopulation *pop, Grid *Pgrad,Grid *divBulkV){

	int nSpecies = pop->nSpeciesNeutral;
	int nDims = 3; // pop->nDims; // hard-coding allows compiler to replace by value
	double *pos = pop->pos;
	double *vel = pop->vel;

	long int *sizeProd = Pgrad->sizeProd;
	double *val = Pgrad->val;

	for(int s=0;s<nSpecies;s++){

		gMul(Pgrad, 1./pop->mass[s]);

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p+=nDims){
			double dv[3];
			double divergence[1];
			neInterp3D1(dv,&pos[p],val,sizeProd);
			neInterp3D1scalar(divergence,&pos[p],val,sizeProd);
			for(int d=0;d<nDims;d++) vel[p+d] += (dv[d]+vel[p+d]*divergence[0]);
		}

		gMul(Pgrad, pop->mass[s]);
	}
}



//###########################
// Mover
//###########################

void neMove(NeutralPopulation *pop){

	int nSpecies = pop->nSpeciesNeutral;

	int nDims = pop->nDims;
	double *pos = pop->pos;
	double *vel = pop->vel;

	for(int s=0; s<nSpecies; s++){

		long int pStart = pop->iStart[s]*nDims;
		long int pStop = pop->iStop[s]*nDims;

		for(long int p=pStart;p<pStop;p++){
			pos[p] += vel[p];
		}
	}
}


/******************************************************************************
 *	FINITE DIFFERENCE
 *****************************************************************************/

 void divFinDiff1st(const Grid *scalar, Grid *field){

	// Performs first order centered finite difference on field and returns a scalar

	gZero(scalar);
	int rank = scalar->rank;
	// int *size = scalar->size;
	long int *sizeProd = scalar->sizeProd;
	long int *fieldSizeProd = field->sizeProd;

	double *scalarVal = scalar->val;
	double *fieldVal = field->val;

 	// Scalar indices
	long int fNext, fPrev;
	long int s;
	//int fNext = fieldSizeProd[1];

	long int start = alSum(&sizeProd[1], rank-1 );
	long int end = sizeProd[rank]-start;


	// Centered Finite difference
	for(int d = 1; d < rank; d++){
		fNext = start*fieldSizeProd[1] + fieldSizeProd[d];
		fPrev = start*fieldSizeProd[1] - fieldSizeProd[d];
		s = start;



		for(int g = start; g < end; g++){
			scalarVal[s] += 0.5*(fieldVal[fNext] - fieldVal[fPrev]);
			fNext+=rank-1;
			fPrev+=rank-1;
			s ++; // fNext;
			//printf("node: %li, using %li, and %li \n",s,fNext,fPrev);
			//printf("node: %li, using %li, and %li \n",s,fNext,fPrev);
			//printf("fieldSizeProd = %li, scalarsizeProd = %li \n",fieldSizeProd[4],sizeProd[4]);
		}
	}

	return;
}


//#############################
// Migration
//#############################


funPtr neExtractEmigrants3DOpen_set(const dictionary *ini){
	int nDims = iniGetInt(ini, "grid:nDims");
	if(nDims!=3) msg(ERROR, "neExtractEmigrants3DOpen requires grid:nDims=3");
	return neExtractEmigrants3DOpen;
}
void neExtractEmigrants3DOpen(NeutralPopulation *pop, MpiInfo *mpiInfo){

	//TODO: Needs more testing
	int nSpecies = pop->nSpeciesNeutral;
	double *pos = pop->pos;
	double *vel = pop->vel;
	double *thresholds = mpiInfo->thresholds;
	const int neighborhoodCenter = 13;
	long int *nEmigrants = mpiInfo->nEmigrants;
	int nNeighbors = mpiInfo->nNeighbors;
	int *trueSize = mpiInfo->trueSize;
	int *nSubdomainsProd = mpiInfo->nSubdomainsProd;
	int *nSubdomains = malloc(3*sizeof(*nSubdomains));
	int *subdomain = mpiInfo->subdomain;
	bndType *bnd = pop->bnd;
	//int rank = mpiInfo->mpiRank;

	//printf("mpirank = %i, bnd[1] = %d,bnd[2] = %d,bnd[3] = %d, bnd[5] = %d,bnd[6] = %d,bnd[7] = %d \n",rank,bnd[1],bnd[2],bnd[3],bnd[5],bnd[6],bnd[7]);
	double dummyPos[3];

	int *offset = mpiInfo->offset;
	//int nDims = pop->nDims;


	nSubdomains[0] = nSubdomainsProd[1];
	nSubdomains[1] = nSubdomainsProd[2]/nSubdomainsProd[1];
	nSubdomains[2] = (nSubdomainsProd[3]/nSubdomainsProd[2]);

	//msg(STATUS, "sssss %i, %i, %i", nSubdomains[0], nSubdomains[1], nSubdomains[2]);
	//exit(0);

	//msg(STATUS,"neighbours = %i",mpiInfo->nNeighbors);

	// By using the dummy to hold data we won't lose track of the beginning of
	// the arrays when incrementing the pointer
	double **emigrants = mpiInfo->emigrantsDummy;
	for(int ne=0;ne<nNeighbors;ne++){
		emigrants[ne] = mpiInfo->emigrants[ne];

	}
	alSetAll(nEmigrants,nSpecies*nNeighbors,0);

	double lx = thresholds[0];
	double ly = thresholds[1];
	double lz = thresholds[2];
	double ux = thresholds[3];
	double uy = thresholds[4];
	double uz = thresholds[5];

	//adPrint(thresholds,6);

	for(int s=0;s<nSpecies;s++){

		long int pStart = pop->iStart[s]*3;
		long int pStop = pop->iStop[s]*3;
		long int removedupp = 0; //debug
		long int removedlow = 0; //debug
		long int exhanged = 0; //debug
		for(long int p=pStart;p<pStop;p+=3){

			for(int d=0;d<3;d++) {
				//Why is offset size -1 ? ... -1 ?
				//printf("%i \n",offset[d]);
				dummyPos[d] = pop->pos[p+d] + offset[d];
				//if(pop->pos[p+d]<1) printf("\n \n pop->pos[p+d] = %f \n \n ",pop->pos[p+d]);
			}


			// not offset but GLOBAL Size!
			//MpiInfo->subdomain;				///< MPI node (nDims elements)
			//MpiInfo->nSubdomains;


			if ( (dummyPos[0] > trueSize[0]*(nSubdomains[0]) && bnd[5]!=PERIODIC)
				|| ( dummyPos[1] > trueSize[1]*(nSubdomains[1]) && bnd[6]!=PERIODIC)
				|| (  dummyPos[2] > trueSize[2]*(nSubdomains[2]) && bnd[7]!=PERIODIC)   ){

					// if (s==0){
					// 	printf("removed upper \n");
					// 	printf("trueSize[1]*(nSubdomains[0]) = %i \n",trueSize[0]*(nSubdomains[0]));
					// 	printf("bnd[5] = %d,bnd[6] = %d,bnd[7] = %d \n",bnd[5],bnd[6],bnd[7]);
					// 	printf("global: %f, %f, %f \n", dummyPos[0], dummyPos[1], dummyPos[2]);
					// 	printf("local: %f, %f, %f \n", pop->pos[p+0], pop->pos[p+1],  pop->pos[p+2]);
					// 	printf("ofsett: %i, %i, %i \n \n", offset[0], offset[1], offset[2]);
					// }
				//printf("too large \n");
				//msg(STATUS,"%i, %i, %i \n",trueSize[0],trueSize[1],trueSize[2]);
				//printf("global: %f, %f, %f \n",dummyPos[0], dummyPos[1], dummyPos[2]);
				//printf("Local: %f, %f, %f \n",pop->pos[p+0], pop->pos[p+1], pop->pos[p+2]);
				//printf("Boundary %f, %f, %f \n \n",ux+trueSize[0]*(nSubdomains[0]-1),uy+trueSize[1]*(nSubdomains[1]-1),uz+trueSize[2]*(nSubdomains[2]-1) );
				//msg(STATUS, " removing ");
				// Remove particle out of bounds particle

				//msg(STATUS,"dummyPos[0] = %f, thresholds[0]+pos = %f",dummyPos[0],thresholds[0]+pop->pos[p*3*+d]);
				//msg(STATUS,"dummyPos[1] = %f, thresholds[1]+pos = %f",dummyPos[1],thresholds[1]+pop->pos[p*3*+d]);
				//msg(STATUS,"dummyPos[2] = %f, thresholds[2]+pos = %f",dummyPos[2],thresholds[2]+pop->pos[p*3*+d]);
				//printf(" \n");
				removedupp += 1; //debug
				pos[p]   = pos[pStop-3];
				pos[p+1] = pos[pStop-2];
				pos[p+2] = pos[pStop-1];
				vel[p]   = vel[pStop-3];
				vel[p+1] = vel[pStop-2];
				vel[p+2] = vel[pStop-1];

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;

			}
			else if ( (dummyPos[0] < -1. && bnd[1]!=PERIODIC)
				|| ( dummyPos[1] < -1. && bnd[2]!=PERIODIC)
				|| ( dummyPos[2] < -1. && bnd[3]!=PERIODIC) ){
				//printf("%f\n",dummyPos[0]);
				//msg(STATUS, " removing ");
				// Remove particle out of bounds particle

				// if (s==0){
				// 	printf("removed lower \n");
				// 	printf("bnd[1] = %d,bnd[2] = %d,bnd[3] = %d \n",bnd[1],bnd[2],bnd[3]);
				// 	printf("global: %f, %f, %f \n", dummyPos[0], dummyPos[1], dummyPos[2]);
				// 	printf("local: %f, %f, %f \n", pop->pos[p+0], pop->pos[p+1],  pop->pos[p+2]);
				// 	printf("ofsett: %i, %i, %i \n \n", offset[0], offset[1], offset[2]);
				// }
				//printf("too small \n");
				//printf("%f, %f, %f \n",dummyPos[0], dummyPos[1], dummyPos[2]);
				//printf("%f, %f, %f \n \n",pop->pos[p+0], pop->pos[p+1], pop->pos[p+2]);

				removedlow += 1; //debug
				pos[p]   = pos[pStop-3];
				pos[p+1] = pos[pStop-2];
				pos[p+2] = pos[pStop-1];
				vel[p]   = vel[pStop-3];
				vel[p+1] = vel[pStop-2];
				vel[p+2] = vel[pStop-1];

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;

			} else if ( (pos[p] >= ux-1 && bnd[5]==PERIODIC)
				|| ( pos[p+1] >= uy-1 && bnd[6]==PERIODIC)
				|| (  pos[p+2] >= uz-1  && bnd[7]==PERIODIC)
				|| (pos[p] < lx+1 && bnd[1]==PERIODIC)
				|| ( pos[p+1] < ly+1 && bnd[2]==PERIODIC)
				|| ( pos[p+2] < lz+1 && bnd[3]==PERIODIC) ){
					//Should look for a better implementation ^

			double x = pos[p];
			double y = pos[p+1];
			double z = pos[p+2];
			int nx = - (x<lx) + (x>=ux);
			int ny = - (y<ly) + (y>=uy);
			int nz = - (z<lz) + (z>=uz);
			int ne = neighborhoodCenter + nx + 3*ny + 9*nz;

			// if (s==1){
			// 	if (x<lx) msg(STATUS,"exhanged particle backward, s = %i \n",s);
			// 	if (y<ly) msg(STATUS,"exhanged particle backward, s = %i \n",s);
			// 	if (z<lz) msg(STATUS,"exhanged particle backward, s = %i \n",s);
			// }
			// if(p==371*3)
			// 	msg(STATUS,"x1: %f",x);

			if(ne!=neighborhoodCenter){
				//msg(STATUS, "exchanged");
				//msg(STATUS,"ne = %i",ne);
				//if(mpiInfo->mpiRank == 0){
					//printf("pos = %f, %f, %f \n",x,y,z);
					//printf("sending to %i\n \n",ne/3);
				//}
				// if (s==1){
				// 	printf("exhcanged \n" );
				// 	printf("global: %f, %f, %f \n", dummyPos[0], dummyPos[1], dummyPos[2]);
				// 	printf("local: %f, %f, %f \n", pop->pos[p+0], pop->pos[p+1],  pop->pos[p+2]);
				// 	printf("offset: %i, %i, %i \n \n", offset[0]+1, offset[1]+1, offset[2]+1);
				// }

				exhanged += 1;
				*(emigrants[ne]++) = x;
				*(emigrants[ne]++) = y;
				*(emigrants[ne]++) = z;
				*(emigrants[ne]++) = vel[p];
				*(emigrants[ne]++) = vel[p+1];
				*(emigrants[ne]++) = vel[p+2];
				nEmigrants[ne*nSpecies+s]++;
				//if(mpiInfo->mpiRank == 2);
				//printf("nEmigrants[ne*nSpecies+s] = %li\n",nEmigrants[ne*nSpecies+s]);

				pos[p]   = pos[pStop-3];
				pos[p+1] = pos[pStop-2];
				pos[p+2] = pos[pStop-1];
				vel[p]   = vel[pStop-3];
				vel[p+1] = vel[pStop-2];
				vel[p+2] = vel[pStop-1];

				// if(p==371*3)
				// 	msg(STATUS,"x2: %f",pos[p]);

				pStop -= 3;
				p -= 3;
				pop->iStop[s]--;

				}

			}
		}//printf("removedupp = %li,removedlow = %li, exhanged = %li s = %i, for rank = %i \n", removedupp,removedlow,exhanged,s,mpiInfo->mpiRank);

		//msg(STATUS,"pRange: %li, iStop: %li",pStart-pStop,pop->iStop[s]);
	}
	free(nSubdomains);
}


static inline void neShiftImmigrants(MpiInfo *mpiInfo, Grid *grid, int ne,int nSpecies){

	double *immigrants = mpiInfo->immigrants;
	//int nSpecies = 1;//mpiInfo->nSpecies;
	long int nImmigrantsTotal = alSum(&mpiInfo->nImmigrants[ne*nSpecies],nSpecies);
	int nDims = mpiInfo->nDims;

	for(int d=0;d<nDims;d++){
		int n = ne%3-1;
		ne /=3;

		double shift = n*grid->trueSize[d+1];

		for(int i=0;i<nImmigrantsTotal;i++){

			//msg(STATUS," immigrants[d+2*nDims*i] = %f ",immigrants[d+2*nDims*i]);

			immigrants[d+2*nDims*i] += shift;

			double pos = immigrants[d+2*nDims*i];
			//msg(STATUS," ne = %i, n = %i ",ne,n);
			//msg(STATUS," shift = %f, pos = %f \n",shift,pos);

			if(pos>grid->trueSize[d+1]+2){
			 msg(ERROR,"particle %i skipped two domains, pos = %f",i,pos);
			}

		}

	}


}

static inline void neImportParticles(NeutralPopulation *pop, double *particles, long int *nParticles, int nSpecies){

	int nDims = pop->nDims;
	long int *iStop = pop->iStop;


	for(int s=0;s<nSpecies;s++){
		//printf("nParticles[s] = %li \n",nParticles[s]);
		double *pos = &pop->pos[nDims*iStop[s]];
		double *vel = &pop->vel[nDims*iStop[s]];

		for(int i=0;i<nParticles[s];i++){
			for(int d=0;d<nDims;d++) *(pos++) = *(particles++);
			for(int d=0;d<nDims;d++){
				*(vel++) = *(particles++);
				//printf("d = %i, pos = %f \n",d,particles[-1]);
			}
			//printf("\n");
		}

		iStop[s] += nParticles[s];
	}


}

static inline void neExchangeMigrants(NeutralPopulation *pop, MpiInfo *mpiInfo, Grid *grid){

	int nSpecies = pop->nSpeciesNeutral;
	int nNeighbors = mpiInfo->nNeighbors;
	long int nImmigrantsAlloc = mpiInfo->nImmigrantsAlloc;
	int nDims = mpiInfo->nDims;
	double **emigrants = mpiInfo->emigrants;
	double *immigrants = mpiInfo->immigrants;
	long int *nImmigrants = mpiInfo->nImmigrants;
	MPI_Request *send = mpiInfo->send;

	//printf("nImmigrantsAlloc = %li \n",nImmigrantsAlloc);

	for(int ne=0;ne<nNeighbors;ne++){
		if(ne!=mpiInfo->neighborhoodCenter){
			int rank = puNeighborToRank(mpiInfo,ne);
			int reciprocal = puNeighborToReciprocal(ne,nDims);
			long int *nEmigrants  = &mpiInfo->nEmigrants[nSpecies*ne];
			long int length = alSum(nEmigrants,nSpecies)*2*nDims;
			//printf("length = %li \n",length);
			//printf("nEmigrants = %i, nSpecies = %li \n",nEmigrants[0],nSpecies);
			MPI_Isend(emigrants[ne],length,MPI_DOUBLE,rank,reciprocal,MPI_COMM_WORLD,&send[ne]);
		}
	}

	// Since "immigrants" is reused for every receive operation MPI_Irecv cannot
	// be used. However, in order to receive and process whichever comes first
	// MPI_ANY_SOURCE is used.
	for(int a=0;a<nNeighbors-1;a++){

		MPI_Status status;
		MPI_Recv(immigrants,nImmigrantsAlloc,MPI_DOUBLE,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		int ne = status.MPI_TAG;	// Which neighbor it is from equals the tag

		// adPrint(mpiInfo->immigrants,6);
		neShiftImmigrants(mpiInfo,grid,ne,nSpecies);
		// adPrint(mpiInfo->immigrants,6);
		neImportParticles(pop,immigrants,&nImmigrants[ne*nSpecies],nSpecies);

	}

	MPI_Waitall(nNeighbors,send,MPI_STATUS_IGNORE);
}

// TODO: Add fault-handling in case of too small Population struct
static inline void neExchangeNMigrants(NeutralPopulation *pop, MpiInfo *mpiInfo){

	int nSpecies = pop->nSpeciesNeutral;
	int nNeighbors = mpiInfo->nNeighbors;
	MPI_Request *send = mpiInfo->send;
	MPI_Request *recv = mpiInfo->recv;

	// Send number of emigrants and receive number of immigrants.
	// Order of reception is not necessarily as determined by the for-loop since
	// we're using non-blocking send/receive.
	for(int ne=0;ne<nNeighbors;ne++){
		if(ne!=mpiInfo->neighborhoodCenter){
			int rank = puNeighborToRank(mpiInfo,ne);
			int reciprocal = puNeighborToReciprocal(ne,mpiInfo->nDims);
			//msg(STATUS,"ne =%i, rank = %i, reciprocal = %i",ne,rank,reciprocal);
			long int *nEmigrants  = &mpiInfo->nEmigrants[nSpecies*ne];
			long int *nImmigrants = &mpiInfo->nImmigrants[nSpecies*ne];
			// We get segfault on large simulation systems, as a safety measure we can
			// try to use a different tag here than in exchangeMigrants()
			// hence the "2*nNeighbors"
			MPI_Isend(nEmigrants ,nSpecies,MPI_LONG,rank,reciprocal+2*nNeighbors,MPI_COMM_WORLD,&send[ne]);
			MPI_Irecv(nImmigrants,nSpecies,MPI_LONG,rank,ne+2*nNeighbors        ,MPI_COMM_WORLD,&recv[ne]);
		}
	}

	MPI_Waitall(nNeighbors,send,MPI_STATUS_IGNORE);
	MPI_Waitall(nNeighbors,recv,MPI_STATUS_IGNORE);

}

void neMigrate(NeutralPopulation *pop, MpiInfo *mpiInfo, Grid *grid){

	MPI_Barrier(MPI_COMM_WORLD);
	neExchangeNMigrants(pop,mpiInfo);
	MPI_Barrier(MPI_COMM_WORLD); //debug
	neExchangeMigrants(pop,mpiInfo,grid);
	MPI_Barrier(MPI_COMM_WORLD);

}



//#############################
// Pressure solver
//#############################

void nePressureSolve3D(Grid *rhoNeutral,Grid *P,NeutralPopulation *pop, const MpiInfo *mpiInfo){

	gHaloOp(setSlice, rhoNeutral, mpiInfo, TOHALO);
	gZero(P);
	//gZero(rhoNeutral);

    int *nGhostLayers = rhoNeutral->nGhostLayers;
	int *trueSize=rhoNeutral->trueSize;
	long int *sizeProd =  rhoNeutral->sizeProd;
	double stiffC = pop->stiffnessConstant;
	double rho0 = pop->rho0;
	//Seperate values
	double *PVal = P->val;
	double *rhoVal = rhoNeutral->val;
	double val = 0;

    //aiPrint(nGhostLayers,8);
	int index = 0;
	//printf("sizeProd[1] = %i, %i, %i, %i \n",sizeProd[1],sizeProd[2],sizeProd[3],sizeProd[4]);
	for(int k=0;k<trueSize[3]+nGhostLayers[7];k++){
		for(int j=0;j<trueSize[2]+nGhostLayers[6];j++){
			for(int i=0;i<trueSize[1]+nGhostLayers[5];i++){
				index = i+j*sizeProd[2]+k*sizeProd[3];
				val = stiffC*(pow( (rhoVal[index]/rho0),7) - 1); //rho0
				PVal[index] = val;
				//if (val >1000){
				//printf("PVal[index] = %f, rhoval = %f \n",PVal[index],rhoVal[index]);
				//}
				//printf(" %f \n",rhoVal[index]);
				//printf("index = %li, sizeprod = %li \n",index,sizeProd[4]);
			}

		}
	}
	//for(int i=0;i<sizeProd[4];i++){
	//	PVal[i] = stiffC*(pow( (rhoVal[i]/rho0),7) - 1);
	//printf("index = %li, sizeprod = %li \n",index,sizeProd[4]);
	//printf("PVal = %f \n",PVal[i] );
	//}
	gHaloOp(setSlice, P, mpiInfo, TOHALO);
	nuGBnd(P, mpiInfo);
}



void neSetBndSlices(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo){

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *nGhostLayers = grid->nGhostLayers;
	//double *bndSolution = grid->bndSolution;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.
	int nSpecies = iniGetInt(ini,"collisions:nSpeciesNeutral");
	double *velDrift = iniGetDoubleArr(ini,"collisions:neutralDrift",nDims*nSpecies);
	double *B = iniGetDoubleArr(ini,"fields:BExt",nDims);

	double veld[nDims];

	for (int s=0;s<nSpecies;s++){
		for (int d = 0;d<nDims;d++){
			//msg(STATUS,"d = %i, d+s*nDims = %i",d,d+s*nDims);
			veld[d] += (1./nSpecies)*velDrift[d+s*nDims];
		}
	}

	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veld[0],veld[1],veld[2]);

	//double B[3] = {1., 0., 0.};
	//double veld[3] = {0., 1., 1.};
	double veldCrossB[3] = {0., 0., 0.};
	adCrossProd(veld, B, veldCrossB);

	//printf("B[0] = %f, B[1] = %f, B[2] = %f \n",B[0],B[1],B[2]);
	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veldCrossB[0],veldCrossB[1],veldCrossB[2]);
	//printf("veldCrossB = %f,%f,%f",veldCrossB[0],veldCrossB[1],veldCrossB[2]);


	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//double constant1 = 0.; //Dirichlet
	double constant2 = 0.; // solution to bnd cond = 0.

	// The dirichlet part in this function ended up realy ugly. I am sorry for
	// that, but it works.

	// indices are indices of slice perp to d
	long int indices[rank];
	long int edge[rank]; //tells us what dimensions are perp (slice)

	//Lower edge
	for(int d = 1; d < rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-1));
			//printf("q = %li, indices[q] = %li \n",q,indices[q]);
		}
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = nGhostLayers[q+1];
					} else{
						// start indices at minimum local indice
						indices[q] = 0;//
					}
				}
				for(int s = 0; s < nSliceMax; s++){

					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = nGhostLayers[q+1];
						}
					}

					bndSlice[s + (nSliceMax * d)] = 0;
					for(int dd = 0;dd<rank-1;dd++){ //dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + (nSliceMax * d)] += 0;

					}
					// counter to increment only in the slice dims, and not the
					// dim perp to slice
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 domain
						if(indices[dd]<(size[dd+1])-1 && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//nGhostLayers[dd+1];
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}
			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){

					// initialize.
					bndSlice[s + (nSliceMax * d)] = constant2;

					//Solution to equation. constant for now

				}
			}
		}
	}

	//Upper edge
	for(int d = rank+1; d < 2*rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-rank-1));
			//printf("d = %i, edge[q] = %li \n",d,edge[q]);
		}
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
						//printf("nSubdomains = %i\n",nSubdomains[q]);
					} else{
						// start indices at minimum
						indices[q] = 0;
					}
				}
				for(int s = 0; s < nSliceMax; s++){
					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
							//printf("nSubdomains = %i\n",nSubdomains[q]);
						}
					}
					bndSlice[s + (nSliceMax * (d))] = 0;
					for(int dd = 0;dd<rank-1;dd++){
						bndSlice[s + (nSliceMax * (d))] +=  0;
					}
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 local domain
						if(indices[dd]<(size[dd+1]-1) && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//(size[dd+1]-nGhostLayers[dd]);
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}

			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = constant2;
				}
			}
		}
	}

	//msg(STATUS,"nSliceMax = %li",nSliceMax);
	//adPrint(&bndSlice[nSliceMax], nSliceMax*(rank));
	free(velDrift);
	//free(B);
	return;
}


void neSetBndSlicesRho(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo){

	int rank = grid->rank;
	int *size = grid->size;
	bndType *bnd = grid->bnd;
	double *bndSlice = grid->bndSlice;
	int *nGhostLayers = grid->nGhostLayers;
	//double *bndSolution = grid->bndSolution;
	int *subdomain = mpiInfo->subdomain;
	int *nSubdomains = mpiInfo->nSubdomains;
	int nDims = mpiInfo->nDims;

	// using ini is slow, but setting boundary slices is likely only done once.




	int nSpecies = iniGetInt(ini,"collisions:nSpeciesNeutral");
	double *velDrift = iniGetDoubleArr(ini,"collisions:neutralDrift",nDims*nSpecies);
	double *B = iniGetDoubleArr(ini,"fields:BExt",nDims);

	double veld[nDims];

	double *nt = iniGetDoubleArr(ini,"collisions:numberDensityNeutrals",nSpecies);
	double density =0;
	//needed for neutral bnd
	for (int s=0;s<nSpecies;s++){
		density += 0;//nt[s]/(2*nSpecies);
	}

	printf("density = %f",density);

	for (int s=0;s<nSpecies;s++){
		for (int d = 0;d<nDims;d++){
			//msg(STATUS,"d = %i, d+s*nDims = %i",d,d+s*nDims);
			veld[d] += (1./nSpecies)*velDrift[d+s*nDims];
		}
	}

	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veld[0],veld[1],veld[2]);

	//double B[3] = {1., 0., 0.};
	//double veld[3] = {0., 1., 1.};
	double veldCrossB[3] = {0., 0., 0.};
	adCrossProd(veld, B, veldCrossB);

	//printf("B[0] = %f, B[1] = %f, B[2] = %f \n",B[0],B[1],B[2]);
	//printf("veld[0] = %f, veld[1] = %f, veld[2] = %f \n",veldCrossB[0],veldCrossB[1],veldCrossB[2]);
	//printf("veldCrossB = %f,%f,%f",veldCrossB[0],veldCrossB[1],veldCrossB[2]);


	//Number of elements in slice
	long int nSliceMax = 0;
	for(int d=1;d<rank;d++){
		long int nSlice = 1;
		for(int dd=0;dd<rank;dd++){
			if(dd!=d) nSlice *= size[dd];
		}
		if(nSlice>nSliceMax) nSliceMax = nSlice;
	}

	//double constant1 = 0.; //Dirichlet
	double constant2 = 0.; // solution to bnd cond = 0.

	// The dirichlet part in this function ended up realy ugly. I am sorry for
	// that, but it works.

	// indices are indices of slice perp to d
	long int indices[rank];
	long int edge[rank]; //tells us what dimensions are perp (slice)

	//Lower edge
	for(int d = 1; d < rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-1));
			//printf("q = %li, indices[q] = %li \n",q,indices[q]);
		}
		if(subdomain[d-1] == 0){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = nGhostLayers[q+1];
					} else{
						// start indices at minimum local indice
						indices[q] = 0;//
					}
				}
				for(int s = 0; s < nSliceMax; s++){

					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = nGhostLayers[q+1];
						}
					}

					bndSlice[s + (nSliceMax * d)] = density;
					for(int dd = 0;dd<rank-1;dd++){ //dot prod of VxB and indices
						// grid indices (i,j,k) are computed locally, so we need to
						// cast them to global frame in the dot product
						bndSlice[s + (nSliceMax * d)] = density;

					}
					// counter to increment only in the slice dims, and not the
					// dim perp to slice
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 domain
						if(indices[dd]<(size[dd+1])-1 && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//nGhostLayers[dd+1];
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}
			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){

					// initialize.
					bndSlice[s + (nSliceMax * d)] = density;

					//Solution to equation. constant for now

				}
			}
		}
	}

	//Upper edge
	for(int d = rank+1; d < 2*rank; d++){
		for(long int q = 0;q<rank-1;q++){
			//initiallize
			//indices[q] = size[q+1]*subdomain[q];
			edge[q] = (q!=(d-rank-1));
			//printf("d = %i, edge[q] = %li \n",d,edge[q]);
		}
		if(subdomain[d-rank-1]==nSubdomains[d-rank-1]-1){
			if(bnd[d] == DIRICHLET){
				for(long int q = 0;q<rank-1;q++){
					//set dim perp to slice indice to fixed val
					if(!edge[q]){
						indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
						//printf("nSubdomains = %i\n",nSubdomains[q]);
					} else{
						// start indices at minimum
						indices[q] = 0;
					}
				}
				for(int s = 0; s < nSliceMax; s++){
					for(long int q = 0;q<rank-1;q++){
						//set dim perp to slice indice to fixed val
						if(!edge[q]){
							indices[q] = (size[q+1]-2*nGhostLayers[q+1]); //-nGhostLayers[q]
							//printf("nSubdomains = %i\n",nSubdomains[q]);
						}
					}
					bndSlice[s + (nSliceMax * (d))] = density;
					for(int dd = 0;dd<rank-1;dd++){
						bndSlice[s + (nSliceMax * (d))] =  density;
					}
					bool incremented = false;
					for(int dd = 0;dd<rank;dd++){
						//runs up to 34-1 for 32x32x32 local domain
						if(indices[dd]<(size[dd+1]-1) && edge[dd]==1 && incremented == false){
							indices[dd]++;
							incremented = true;
						}else if(incremented == false){
							// reset
							indices[dd] = 0;//(size[dd+1]-nGhostLayers[dd]);
						}
					}
					//printf("indices[0] = %li,indices[1] = %li,indices[2] = %li \n",indices[0],indices[1],indices[2]);
				}
			}

			if(bnd[d] == NEUMANN){
				for(int s = 0; s < nSliceMax; s++){
					bndSlice[s + (nSliceMax * d)] = density;
				}
			}
		}
	}

	//msg(STATUS,"nSliceMax = %li",nSliceMax);
	//adPrint(&bndSlice[nSliceMax], nSliceMax*(rank));
	free(velDrift);
	//free(B);
	return;
}




//#########################################
// Object functions
//#########################################


void nuObjectpurge(NeutralPopulation *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo) {


	//int rank = mpiInfo->mpiRank;
	int size = mpiInfo->mpiSize;

	double *val = rhoObj->val;
	long int *sizeProd = rhoObj->sizeProd;
	long int nDims = pop->nDims;

	int nSpecies = pop->nSpeciesNeutral;

	long int *lookupIntOff = obj->lookupInteriorOffset;
	long int *lookupSurfOff = obj->lookupSurfaceOffset;


	int cutNumber = 0;
	for(int s=0;s<nSpecies;s++) {

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[i*nDims];
			double *vel = &pop->vel[i*nDims];

			// Integer parts of position
			int j = (int) pos[0];
			int k = (int) pos[1];
			int l = (int) pos[2];

			long int p = j + k*sizeProd[2] + l*sizeProd[3];
			long int pIndex = i*nDims; //j + k*sizeProd[2] + l*sizeProd[3];
			//msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
			// Check whether p is one of the object nodes and collect the charge if so.
			for (long int a=0; a<obj->nObjects; a++) {
				for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
					if ((obj->lookupInterior[b])==p) {
						//msg(STATUS,"p, pIndex: %li,%li, %li",p,(pIndex-iStart*nDims),(iStop-iStart));
						//msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
						//msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);
						neCutParticle(pop, s, pIndex, pop->pos, pop->vel);
						cutNumber += 1;
						//msg(STATUS,"iStop = %li",iStop);
						iStop--;

					}
				}

			}

		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &cutNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	//printf("cutNumber = %i \n",cutNumber);
	cutNumber = 0;


}


void nuObjectCollide(NeutralPopulation *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo) {


	//int rank = mpiInfo->mpiRank;
	int size = mpiInfo->mpiSize;

	double *val = rhoObj->val;
	long int *sizeProd = rhoObj->sizeProd;
	long int nDims = pop->nDims;

	int nSpecies = pop->nSpeciesNeutral;

	long int *lookupIntOff = obj->lookupInteriorOffset;
	long int *lookupSurfOff = obj->lookupSurfaceOffset;


	int cutNumber = 0;
	for(int s=0;s<nSpecies;s++) {

		long int iStart = pop->iStart[s];
		long int iStop = pop->iStop[s];

		for(long int i=iStart;i<iStop;i++){

			double *pos = &pop->pos[i*nDims];
			double *vel = &pop->vel[i*nDims];

			// Integer parts of position
			int j = (int) pos[0];
			int k = (int) pos[1];
			int l = (int) pos[2];

			long int p = j + k*sizeProd[2] + l*sizeProd[3];
			long int pIndex = i*nDims; //j + k*sizeProd[2] + l*sizeProd[3];
			//msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
			// Check whether p is one of the object nodes and collect the charge if so.
			for (long int a=0; a<obj->nObjects; a++) {
				for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
					if ((obj->lookupInterior[b])==p) {
						//msg(STATUS,"j,p, pIndex: %i, %li,%li, %li",j,p,(pIndex-iStart*nDims),(iStop-iStart));
						//msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
						//msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);

						//printf("Before: pos = %f,%f,%f \n",pos[0],pos[1],pos[2]);
						//printf("Before: vel = %f,%f,%f \n",vel[0],vel[1],vel[2]);

						neScatterParticle(pop, s, pIndex, pos, vel);

						//printf("After: pos = %f,%f,%f \n",pos[0],pos[1],pos[2]);
						//printf("After: vel = %f,%f,%f \n",vel[0],vel[1],vel[2]);
						//printf("\n");

						cutNumber += 1;
						//msg(STATUS,"iStop = %li",iStop);
						iStop--;



						// Unit testing

						// Integer parts of position
						int j = (int) pos[0];
						int k = (int) pos[1];
						int l = (int) pos[2];

						long int p = j + k*sizeProd[2] + l*sizeProd[3];
						for (long int a=0; a<obj->nObjects; a++) {
							for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
								if ((obj->lookupInterior[b])==p) {
									msg(WARNING,"Neutral particle-object did not collide correct");
									printf("After: pos = %f,%f,%f \n",pos[0],pos[1],pos[2]);
									printf("After: vel = %f,%f,%f \n",vel[0],vel[1],vel[2]);
									printf("\n");
								}
							}
						}// Unit testing
					}
				}

			}

		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &cutNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	msg(STATUS,"collided Number = %i \n",cutNumber);
	cutNumber = 0;


}


void neVelAssertMax(const NeutralPopulation *pop, double max){

	double *vel = pop->vel;

	int nSpecies = pop->nSpeciesNeutral;
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

void nePosAssertInLocalFrame(const NeutralPopulation *pop, const Grid *grid){

	int *size = grid->size;
	double *pos = pop->pos;

	int nSpecies = pop->nSpeciesNeutral;

	int nDims = pop->nDims;

	for(int s=0; s<nSpecies; s++){

		long int iStart = pop->iStart[s];
		long int iStop  = pop->iStop[s];
		for(int i=iStart; i<iStop; i++){

			for(int d=0; d<nDims; d++){

				if(pos[i*nDims+d]>size[d+1]+1 || pos[i*nDims+d]<-1){
					printf("Particle i=%i (of neutral specie %i) is out of bounds"
					"in dimension %i: %f>%i",
					i, s, d, pos[i*nDims+d], size[d+1]-1);
					exit(1);

				}
			}
		}
	}
}
