/**
 * @file		object.c
 * @author		Jan Deca <jandeca@gmail.com>
 * @brief		All object-related functions are here.
 * @date		19.10.16
 */

#include "core.h"
#include "object.h"
#include "pusher.h"
#include "multigrid.h"
#include "spectral.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

/******************************************************************************
 *  LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/
// /**
//  * @brief   Object debug functions
//  * @param	obj		Object
//  * @param	ini		input settings
//  * @return	void
//  */
// static void print_gsl_mat(gsl_matrix_view A);
/**
 * @brief   Count the number of objects and fills the lookup tables.
 * @param	obj		Object
 * @return	void
 */
static void oFillLookupTables(PincObject *obj);

/**
 * @brief   Find all the object nodes which are part of the object surface.
 * @param	obj		Object
 * @param	ini		input settings
 * @return	void
 */
static void oFindObjectSurfaceNodes(PincObject *obj);



/******************************************************************************
 *  LOCAL FUNCTION DEFINITIONS
 *****************************************************************************/
//
// static void print_gsl_mat(const gsl_matrix_view A){
//
//     FILE *f;
//     f = fopen("matrix.txt", "w");
//
//     const gsl_matrix *mat = &A.matrix;
//     double element;
//
//     for(size_t i=0; i<mat->size1; i++){
//
//         for(size_t j=0; j<mat->size2; j++){
//             element = gsl_matrix_get(mat, i, j);
//             fprintf(f, "%.2g\t", element);
//         }
//         fprintf(f, "\n");
//     }
//
//     fclose(f);
// }


//delete this
static void oGhost(long int node, const int *nGhostLayersBefore,
            const int *nGhostLayersAfter, const int *trueSize,
            const long int *sizeProd, bool *ghost) {
    if (*sizeProd==1) {
        if (node>*trueSize || node<1) {
            *ghost = true;
        }
    } else {
        long int help = *(sizeProd) * (*(trueSize)+1);
        if (node < *(sizeProd) || node > help  ) {
            *ghost = true;
        }
        node = node % *(sizeProd);
        oGhost(node,nGhostLayersBefore-1,
               nGhostLayersAfter-1,trueSize-1,sizeProd-1,ghost);
    }
}
//Check whether a certain node is a ghost node.
static bool oIsGhostNode(Grid *grid, long int node) {

    long int *sizeProd = grid->sizeProd;
    int *trueSize = grid->trueSize;
    int *nGhostLayers = grid->nGhostLayers;
    int rank = grid->rank;

    bool ghost = false;

    oGhost(node,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1],&ghost);

  return ghost;
}



// Count the number of objects and fill the lookup tables
static void oFillLookupTables(PincObject *obj) {

    //printf("oFillLookupTables \n");
    int nObjects = obj->nObjects;
    // Initialise and compute the array that stores the offsets of the objects in the lookup table.
    long int *lookupInteriorOffset = malloc((nObjects+1)*sizeof(*lookupInteriorOffset));


    for (long int i=0; i<nObjects+1; i++) {
        lookupInteriorOffset[i] = 0;
    }

    for (long int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>0.5 ){ //&& !oIsGhostNode(obj->domain, i)
            lookupInteriorOffset[(int)(obj->domain->val[i]+0.5)]++;
        }
    }

    alCumSum(lookupInteriorOffset+1,lookupInteriorOffset,nObjects);
    // Initialise and compute the lookup table.
    long int *lookupInterior = malloc( (lookupInteriorOffset[nObjects])*sizeof(*lookupInterior) );
    alSetAll(lookupInterior,lookupInteriorOffset[nObjects],0);

    long int *index = malloc((nObjects)*sizeof(*index));
    for (long int i=0; i<nObjects; i++) {
        index[i]=lookupInteriorOffset[i];
    }

    for (long int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>0.5 ){ //&& !oIsGhostNode(obj->domain, i)
            lookupInterior[(index[(int)(obj->domain->val[i]+0.5)-1])] = i;
            (index[(int)(obj->domain->val[i]+0.5)-1])++;
        }
    }
    //alPrint(lookupInterior,(lookupInteriorOffset[nObjects]) );
    // Add to the object.
    //obj->nObjects = nObjects;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;

    free(index);

}


static long int oGatherSurfaceNodes(PincObject *obj, long int *nodCorLoc, \
  long int *nodCorGlob,long int *lookupSurfOff, const MpiInfo *mpiInfo){


  int size = mpiInfo->mpiSize;

  for (long int a=0; a<obj->nObjects; a++) {

      long int nodesThisCore = lookupSurfOff[a+1] - lookupSurfOff[a];

      //printf("rank = %i nodesThisCore = %li \n",mpiInfo->mpiRank,nodesThisCore);
      // Let every core know how many surface nodes everybody has.
      MPI_Allgather(&nodesThisCore, 1, MPI_LONG, nodCorLoc, 1, MPI_LONG, MPI_COMM_WORLD);

      for(long int i=size-1;i>-1;i--) nodCorLoc[i+1]=nodCorLoc[i];
      nodCorLoc[0] = 0;
      alCumSum(nodCorLoc+1,nodCorLoc,size);

      for (long int b=0; b<size+1; b++) {
        nodCorGlob[a*(size+1)+b] = nodCorLoc[b];

      }
  }

  // Find the size and initialise the array holding the capacitance matrices for all objects.
  long int capMatrixAllSize = 0;
  for (long int a=0; a<obj->nObjects; a++) {
      capMatrixAllSize +=nodCorGlob[a*(size+1)+size];
  }
  return capMatrixAllSize;
}


// Compute the capacitance matrix for each object.
void oComputeCapacitanceMatrix(PincObject *obj, dictionary *ini, const MpiInfo *mpiInfo) {

    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurf = obj->lookupSurface;
    long int *lookupSurfOff = obj->lookupSurfaceOffset;

    double *capMatrixAll = obj->capMatrixAll;
    long int *nodCorGlob = obj->capMatrixAllOffsets;
    double *capMatrixSum = obj->capMatrixSum;

    // Allocate and initialise the structures to run the potential solver.
    void (*solverInterface)() = select(ini, "methods:poisson", mgSolver_set, sSolver_set);
    void (*solve)() = NULL;
    void *(*solverAlloc)() = NULL;
    void (*solverFree)() = NULL;
    solverInterface(&solve, &solverAlloc, &solverFree);

    Grid *rhoCap = gAlloc(ini, SCALAR,mpiInfo);
    Grid *phiCap = gAlloc(ini, SCALAR,mpiInfo);

	double realTol = iniGetDouble(ini, "multigrid:tol");
	double objTol = iniGetDouble(ini, "multigrid:objTol");
	iniSetDouble(ini, "multigrid:tol", objTol);
    void *solver = solverAlloc(ini, rhoCap, phiCap, mpiInfo);
	iniSetDouble(ini, "multigrid:tol", realTol);//1e-6;
	//setting tol for residual only for cap matrix

    //msg(STATUS,"in oComputeCapacitanceMatrix");
    //exit(0);
    // for(int r=0; r<2*phiCap->rank; r++){
    //   rhoCap->bnd[r] = DIRICHLET;
    //   phiCap->bnd[r] = DIRICHLET;
    // }
    //aiPrint(phiCap->bnd,2*phiCap->rank);

    // Set Rho to zero.
    gZero(rhoCap);
	gZero(phiCap);

	//gSetBndSlices(ini, phiCap, mpiInfo);

    // Find the number of surface nodes for each object.




    // Compute the actual capacitance matrix for each object.
    for (long int a=0; a<obj->nObjects; a++) {
        long int j = 0; // Keep track of the rank
        long int inode = 0; // Keep track of the node number


        long int totSNGlob = nodCorGlob[a*(size+1)+size];
        long int beginIndex = nodCorGlob[a*(size+1)+rank];
        long int endIndex = nodCorGlob[a*(size+1)+rank+1];


        // Initialise the matrix and its inverse.
        double *capMatrix = malloc( (totSNGlob*totSNGlob) * sizeof(*capMatrix));
        double *invCapMatrix = malloc( (totSNGlob*totSNGlob) * sizeof(*invCapMatrix));
        adSetAll(capMatrix,totSNGlob*totSNGlob,0);
        adSetAll(invCapMatrix,totSNGlob*totSNGlob,0);

        // Loop over the nodes and fill the matrix
        for (long int i=0; i<totSNGlob; i++) {
			//gZero(phiCap);
            msg(STATUS,"Solving capacitance matrix for node %ld of %ld for object %ld of %ld.", \
                i+1,totSNGlob,a+1,obj->nObjects);

            // Don't loop over cores who do not have any surface nodes.
            while ((nodCorGlob[a*(size+1)+j+1]-nodCorGlob[a*(size+1)+j])==0) {
                j++;
            }

            // Set the surface node to 1 charge.
            if (rank==j) {
                rhoCap->val[lookupSurf[lookupSurfOff[a] + inode]] = 1;
                //printf("adding 1 rho to node %li \n",lookupSurf[lookupSurfOff[a] + inode]);
            }

            // Solve for the potential.
            solve(solver, rhoCap, phiCap, mpiInfo);

            //msg(STATUS,"phi size = %i",phiCap->sizeProd[4]);
        		//for (long int q = 0; q<phiCap->rank;q++){
        		//	adPrint(phiCap->val,phiCap->sizeProd[4] );
        		//	}
            // Set the surface node back to zero.
            if (rank==j) {
                rhoCap->val[lookupSurf[inode]] = 0;
            }
            // Fill column i of the capacitance matrix.
            for (int k=beginIndex; k<endIndex; k++) {
                capMatrix[totSNGlob*k + i] = phiCap->val[lookupSurf[lookupSurfOff[a] + k-beginIndex]];
                //printf("indexing to %i in cap matrix, max: %li \n",(lookupSurf[lookupSurfOff[a] + k-beginIndex]),(phiCap->sizeProd[4]));
            }

            // Increase the counters. If you looped over all nodes on this core, increase the rank and reset inode.
            inode++;
            while (inode>(nodCorGlob[a*(size+1)+j+1]-nodCorGlob[a*(size+1)+j]-1)) {
                j++;
                inode=0;
            }

        }

        // Make sure every codes has the complete matrix (needed for BLAS).
        long int mpiSendNr = (totSNGlob*totSNGlob);
        MPI_Allreduce(MPI_IN_PLACE, capMatrix, mpiSendNr, \
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //printf("(totSNGlob*totSNGlob) = %li \n", mpiSendNr);
        //adPrint(capMatrix, (totSNGlob*totSNGlob));
        // Compute the inverse of the capacitance matrix.
        // Actually, the inverse is the capacitance matrix. Probably have to rethink the variable names.
        gsl_matrix_view A = gsl_matrix_view_array(capMatrix, totSNGlob, totSNGlob);
        gsl_matrix_view invA = gsl_matrix_view_array(invCapMatrix, totSNGlob, totSNGlob);

        //debug 290619
        //print_gsl_mat(A);

        int s;
        gsl_permutation *p = gsl_permutation_alloc(totSNGlob);
        gsl_linalg_LU_decomp(&A.matrix, p, &s);
        gsl_linalg_LU_invert(&A.matrix, p, &invA.matrix);
        //print_gsl_mat(invA);

        // Add the invCapMatrix for object a to the big array.
        for (long int l=0; l<totSNGlob*totSNGlob; l++) {
            capMatrixAll[a*totSNGlob*totSNGlob+l] = invCapMatrix[l];
        }

        // Compute here to total sum of elements in the capacitance matrix (needed later).
        capMatrixSum[a] = adSum(invCapMatrix,totSNGlob*totSNGlob);
        // We need the inverse later on.
        capMatrixSum[a] = 1/capMatrixSum[a];


        free(capMatrix),
        free(invCapMatrix);
        gsl_permutation_free(p);
    }

    //long int *capMatrixAllOffsets = nodCorGlob;

    //adPrint(capMatrixAll,capMatrixAllSize*capMatrixAllSize);
    // Add to object
    //obj->capMatrixAll = capMatrixAll;
    //obj->capMatrixAllOffsets = nodCorGlob;
    //obj->capMatrixSum = capMatrixSum;

    gFree(rhoCap);
    gFree(phiCap);
    solverFree(solver);

    //free(nodCorLoc);

    //free(nodCorGlob);

}

// Construct and solve equation 5 in Miyake_Usui_PoP_2009
void oApplyCapacitanceMatrix(Grid *rho, const Grid *phi, const PincObject *obj, const MpiInfo *mpiInfo,Units *units){

    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurf = obj->lookupSurface;
    long int *lookupSurfOff = obj->lookupSurfaceOffset;

    double *capMatrixAll = obj->capMatrixAll;
    long int *capMatrixAllOffsets = obj->capMatrixAllOffsets;
    double *capMatrixSum = obj->capMatrixSum;
	double *deltaPhi = obj->deltaPhi;
	double *rhoCorr = obj->rhoCorr;

  double *bias = obj->bias;
  int biasOn = obj->biasOn;


    // Loop over the objects
    for (long int a=0; a<obj->nObjects; a++) {

        // This number is in fact the correct potential of the object.
        double capMatrixPhiSum = 0;

        // total number of surface nodes
        long int totSNGlob = capMatrixAllOffsets[a*(size+1)+size];
        long int beginIndex = capMatrixAllOffsets[a*(size+1)+rank];
        long int endIndex = capMatrixAllOffsets[a*(size+1)+rank+1];

        //double *deltaPhi = malloc(totSNGlob*sizeof(*deltaPhi));
        adSetAll(deltaPhi,totSNGlob,0);
        //double *rhoCorr = malloc(totSNGlob*sizeof(*rhoCorr));
        adSetAll(rhoCorr,totSNGlob,0);


        // Compute eq. 7.
        for (long int i=0; i<totSNGlob; i++) {
            // Make sure that each core loops only over the matrix elements/parts of the grid it has
            for (long int j=beginIndex; j<endIndex; j++) {
                capMatrixPhiSum += capMatrixAll[a*totSNGlob*totSNGlob+totSNGlob*j+i] \
                * (phi->val[lookupSurf[lookupSurfOff[a] + j-beginIndex]]);
            }
        }

		// Debug test
		// capMatrixSum[a] = 0;
		// for (long int i=0; i<totSNGlob; i++) {
        //     // Make sure that each core loops only over the matrix elements/parts of the grid it has
        //     for (long int j=beginIndex; j<endIndex; j++) {
        //         capMatrixSum[a] += capMatrixAll[a*totSNGlob*totSNGlob+totSNGlob*j+i];
        //     }
        // }
		// //MPI_Allgather(&capMatrixSum[a], 1, MPI_DOUBLE, &capMatrixSum[a], 1, MPI_DOUBLE, MPI_COMM_WORLD);
		// MPI_Allreduce(MPI_IN_PLACE, &capMatrixSum[a], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // This is phi_c for each object.
        if (biasOn==1){
          //printf("AKSJDHSKAJDHAK");
          capMatrixPhiSum = bias[a];
        } else{
          capMatrixPhiSum = capMatrixPhiSum*capMatrixSum[a];
          MPI_Allreduce(MPI_IN_PLACE, &capMatrixPhiSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

		//printf("capMatrixSum[a] = %f\n",capMatrixSum[a] );
		//printf("capMatrixPhiSum = %f\n",capMatrixPhiSum );
        msg(STATUS,"Potential-check for object %ld : %f",a,(units->potential*capMatrixPhiSum));
        //capMatrixPhiSum=0.03;

        for (long int j=beginIndex; j<endIndex; j++) {
            deltaPhi[j] = capMatrixPhiSum - phi->val[lookupSurf[lookupSurfOff[a] + j-beginIndex]];
            //printf("adding correction to node %f \n",deltaPhi[j] );
        }

        MPI_Allreduce(MPI_IN_PLACE, deltaPhi, totSNGlob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Eq. 5
        for (long int i=0; i<totSNGlob; i++) {
            for (long int j=beginIndex; j<endIndex; j++) {
                rhoCorr[i] += capMatrixAll[a*totSNGlob*totSNGlob+totSNGlob*j+i]*deltaPhi[j];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, rhoCorr, totSNGlob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Add the charge corrections.
        for (long int j=beginIndex; j<endIndex; j++) {
            rho->val[lookupSurf[lookupSurfOff[a] + j-beginIndex]] += rhoCorr[j];

        }

        //free(deltaPhi);
        //free(rhoCorr);
    }


}

//Find all the object nodes which are part of the object surface.
static void oFindObjectSurfaceNodes(PincObject *obj) {

    //printf("in oFindObjSurf \n");
    long int *sizeProd = obj->domain->sizeProd;
    double *val = obj->domain->val;

    // Initialise the array storing the offsets for the surface nodes in the lookup table.
    long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*sizeof(*lookupSurfaceOffset));
    alSetAll(lookupSurfaceOffset,obj->nObjects+1,0);

    //printf("finding offsets \n");
    // Find the 8 neighbour cells of each non-ghost node.
    long int *myNB = malloc(10*sizeof(*myNB));
    // Find the ofsetts first.
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<sizeProd[obj->domain->rank]; b++) {
            if (!oIsGhostNode(obj->domain, b)) {
                myNB[0] = b;    // me on node i,j,k
                myNB[1] = myNB[0];                  // cell i,j,k
                myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
                myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
                myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
                myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];      // cell i-1,j-1,k
                myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];   // cell i-1,j-1,k-1
                //myNB[9] = myNB[0] - sizeProd[1] - sizeProd[2];      // cell i-1,j-1,k

                int d=0;
                if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
                if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
                if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
                if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
                if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
                if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
                if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
                if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;
                //if (val[myNB[9]]>(a+0.5) && val[myNB[9]]<(a+1.5)) d++;

                // double x = pos[p];
          			// double y = pos[p+1];
          			// double z = pos[p+2];
          			// int nx = - (x<lx) + (x>=ux);
          			// int ny = - (y<ly) + (y>=uy);
          			// int nz = - (z<lz) + (z>=uz);
          			// int ne = neighborhoodCenter + nx + 3*ny + 9*nz;

                // Check if on surface.
                if (d<7.5 && d>0) { //val[myNB[0]]>(a+0.5) &&
                    lookupSurfaceOffset[a+1]++;

                }
            }
        }
        //MPI_Allreduce(MPI_IN_PLACE, &lookupSurfaceOffset[a+1], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //aiPrint(&lookupSurfaceOffset[a+1],1);
    }
    //printf("offsets done \n");
    alCumSum(lookupSurfaceOffset+1,lookupSurfaceOffset,obj->nObjects);

    // Initialise and compute the lookup table.
    long int *lookupSurface = malloc((lookupSurfaceOffset[obj->nObjects])*sizeof(*lookupSurface));
    alSetAll(lookupSurface,lookupSurfaceOffset[obj->nObjects],0);

    long int *index = malloc((obj->nObjects)*sizeof(*index));
    for (long int i=0; i<obj->nObjects; i++) {
        index[i]=lookupSurfaceOffset[i];

    }

    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
            if (!oIsGhostNode(obj->domain, b)) {
                myNB[0] = b;    // me on node i,j,k
                myNB[1] = myNB[0];                  // cell i,j,k
                myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
                myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
                myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
                myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];                  // cell i-1,j-1,k;
                myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];    // cell i-1,j-1,k-1
                //myNB[9] = myNB[0] - sizeProd[1] - sizeProd[2];      // cell i-1,j-1,k

                int d=0;
                if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
                if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
                if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
                if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
                if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
                if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
                if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
                if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;
                //if (val[myNB[9]]>(a+0.5) && val[myNB[9]]<(a+1.5)) d++;

                // Check if on surface.
                if (d<7.5 && d>0) { //val[myNB[0]]>(a+0.5) &&
                    lookupSurface[index[a]] = myNB[0];
                    index[a]++;
                    //printf("indexing %li, MAX: %li \n",index[a],lookupSurfaceOffset[obj->nObjects]+1);
                }


            }
        }
    }
    //printf("lookup surface done \n");
    // Add to object.
    //alPrint(lookupSurface,(lookupSurfaceOffset[obj->nObjects]+1));
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;

    free(myNB);
    free(index);
}
// //Find all the object nodes which are part of the object surface.
// void oFindObjectSurfaceNodes(PincObject *obj, const MpiInfo *mpiInfo) {
//
//     long int *sizeProd = obj->domain->sizeProd;
//     double *val = obj->domain->val;
//
//     // Initialise the array storing the offsets for the surface nodes in the lookup table.
//     long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*sizeof(*lookupSurfaceOffset));
//     alSetAll(lookupSurfaceOffset,obj->nObjects+1,0);
//
//     // Find the 8 neighbour cells of each non-ghost node.
//     long int *myNB = malloc(9*sizeof(*myNB));
//     // Find the ofsetts first.
//     for (long int a=0; a<obj->nObjects; a++) {
//         for (long int b=0; b<sizeProd[obj->domain->rank]; b++) {
//             if (!oIsGhostNode(obj->domain, b)) {
//                 myNB[0] = b;    // me on node i,j,k
//                 myNB[1] = myNB[0];                  // cell i,j,k
//                 myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
//                 myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
//                 myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
//                 myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
//                 myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
//                 myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];      // cell i-1,j-1,k
//                 myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];   // cell i-1,j-1,k-1
//
//                 int d=0;
//                 if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
//                 if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
//                 if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
//                 if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
//                 if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
//                 if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
//                 if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
//                 if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;
//                 printf("d = %i \n",d);
//                 // Check if on surface.
//                 if ( val[myNB[0]]>(a+0.5) && d<8.5 && d>0) {
//                     lookupSurfaceOffset[a+1]++;
//                     //msg(STATUS,"lookupSurfaceOffset[a+1] = %li",lookupSurfaceOffset[a+1]);
//                 }
//             }
//         }lookupSurfaceOffset[a+1]++;
//     }
//     alCumSum(lookupSurfaceOffset+1,lookupSurfaceOffset,obj->nObjects);
//
//     //msg(STATUS,"lookupSurfaceOffset[a+1] = %li",lookupSurfaceOffset[0+1]);
//     // Initialise and compute the lookup table.
//     long int *lookupSurface = malloc((lookupSurfaceOffset[obj->nObjects]+1)*sizeof(*lookupSurface));
//     alSetAll(lookupSurface,lookupSurfaceOffset[obj->nObjects]+1,0);
//
//     long int *index = malloc((obj->nObjects)*sizeof(*index));
//     for (long int i=0; i<obj->nObjects; i++) {
//         index[i]=lookupSurfaceOffset[i];
//     }
//
//     for (long int a=0; a<obj->nObjects; a++) {
//         for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
//             if (!oIsGhostNode(obj->domain, b)) {
//                 myNB[0] = b;    // me on node i,j,k
//                 myNB[1] = myNB[0];                  // cell i,j,k
//                 myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
//                 myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
//                 myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
//                 myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
//                 myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
//                 myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];                  // cell i-1,j-1,k;
//                 myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];    // cell i-1,j-1,k-1
//
//                 //alPrint(myNB,8);
//                 //adPrint(val,sizeProd[obj->domain->rank]);
//                 //msg(STATUS,"a = %li, val[myNB[1]] = %f",a,val[myNB[1]]);
//                 int d=0;
//                 if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
//                 if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
//                 if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
//                 if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
//                 if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
//                 if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
//                 if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
//                 if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;
//
//                 // Check if on surface
//                 //msg(STATUS,"val[myNB[1]] = %f, myNB[1] = %li, d = %i",val[myNB[1]],myNB[1],d);
//                 if ( val[myNB[1]]>(a+0.5) && d<8.5 && d>0) {
//                     lookupSurface[index[a]] = myNB[1];
//                     index[a]++;
//                     //msg(STATUS,"index[a] = %li",index[a]);
//                 }
//             }
//         }
//     }
//
//     // Add to object.
//     //alPrint(lookupSurface,1);
//     obj->lookupSurface = lookupSurface;
//     obj->lookupSurfaceOffset = lookupSurfaceOffset;
//     //msg(STATUS,"EXITING surface lookup");
//     //free(index);
//     //free(myNB);
// }


// Collect the charge inside each object.
void oCollectObjectCharge(Population *pop, Grid *rhoObj, PincObject *obj, const MpiInfo *mpiInfo) {


    //int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;

    double *val = rhoObj->val;
    long int *sizeProd = rhoObj->sizeProd;
    long int nDims = pop->nDims;

    int nSpecies = pop->nSpecies;
    double *charge = pop->charge;

    long int *lookupIntOff = obj->lookupInteriorOffset;
    long int *lookupSurfOff = obj->lookupSurfaceOffset;
    double *objectCurrent = obj->objectCurrent;

    // We might add this to the Object, although probably better  to store the rhoObj for restarts and insulators later on.
    double *chargeCounter = malloc(obj->nObjects*sizeof(*chargeCounter));
    //double *objectCurrent = malloc(obj->nObjects*(nSpecies)*sizeof(*chargeCounter));

    adSetAll(chargeCounter,obj->nObjects,0);//sets charge counter=0 for all objects
    adSetAll(objectCurrent,nSpecies*obj->nObjects,0);

    //long int *nodCorLoc = malloc((size+1)*sizeof(*nodCorLoc));
    long int *nodCorGlob = obj->capMatrixAllOffsets;//malloc(obj->nObjects*(size+1)*sizeof(*nodCorGlob));

    // for (long int a=0; a<obj->nObjects; a++) {
	//
    //     long int nodesThisCore = lookupSurfOff[a+1] - lookupSurfOff[a];
	//
    //     // Let every core know how many surface nodes everybody has.
    //     MPI_Allgather(&nodesThisCore, 1, MPI_LONG, nodCorLoc, 1, MPI_LONG, MPI_COMM_WORLD);
	//
    //     for(long int i=size-1;i>-1;i--) nodCorLoc[i+1]=nodCorLoc[i];
    //     nodCorLoc[0] = 0;
    //     alCumSum(nodCorLoc+1,nodCorLoc,size);
	//
    //     for (long int b=0; b<size+1; b++) nodCorGlob[a*(size+1)+b] = nodCorLoc[b];
    // }
    //printf("obj->nObjects*(size+1) = %li \n",obj->nObjects*(size+1));
    //alPrint(nodCorGlob,obj->nObjects*(size+1));
    //alPrint(nodCorLoc,(size+1));

    //double invNrSurfNod = 1.0/(obj->lookupSurfaceOffset[obj->nObjects]);
    double *invNrSurfNod = obj->invNrSurfNod;//malloc(obj->nObjects*sizeof(*invNrSurfNod));

	adSetAll(invNrSurfNod,obj->nObjects,0);
    for (long int a=0; a<obj->nObjects; a++) {
        invNrSurfNod[a] = 1.0/(nodCorGlob[(a+1)*(size)]);
        //printf("invNrSurfNod[a] = %f, nodCorGlob[(a+1)*(size+1)] = %li",1./invNrSurfNod[a],nodCorGlob[(a+1)*(size)]);
    }

    int cutNumber = 0;
    for(int s=0;s<nSpecies;s++) {

        long int iStart = pop->iStart[s];
        long int iStop = pop->iStop[s];

        for(long int i=iStart;i<iStop;i++){

            double *pos = &pop->pos[i*nDims];
            //double *vel = &pop->vel[i*nDims];

            // Integer parts of position
            int j = (int) pos[0];
            int k = (int) pos[1];
            int l = (int) pos[2];

            long int p = j + k*sizeProd[2] + l*sizeProd[3];
            long int pIndex = i*nDims; //j + k*sizeProd[2] + l*sizeProd[3];
            //msg(STATUS,"i, pIndex: %li,%i",(i-iStart),(pIndex-iStart*nDims));
            // Check whether p is one of the object nodes and collect the charge if so.
            // for (long int a=0; a<obj->nObjects; a++) {
            //     for (long int b=lookupSurfOff[a]; b<lookupSurfOff[a+1]; b++) {
            //         if ((obj->lookupSurface[b])==p) {
            //             chargeCounter[a] += charge[s];
			// 			//printf("adding charge = %f\n",charge[s] );
            //             //msg(STATUS,"p, pIndex: %li,%li, %li",p,(pIndex-iStart*nDims),(iStop-iStart));
            //             //msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
            //             //msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);
            //             pCut(pop, s, pIndex, pop->pos, pop->vel);
            //             cutNumber += 1;
            //             //msg(STATUS,"iStop = %li",iStop);
            //             iStop--;
			//
            //         }
            //     }
			//
            // }
			for (long int a=0; a<obj->nObjects; a++) {
                for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
                    if ((obj->lookupInterior[b])==p) {
                        chargeCounter[a] += charge[s];
                        objectCurrent[a*nSpecies + s] += charge[s];
						//printf("adding charge = %f\n",charge[s] );
                        //msg(STATUS,"p, pIndex: %li,%li, %li",p,(pIndex-iStart*nDims),(iStop-iStart));
                        //msg(STATUS,"j,k,l: %i,%i, %i",j,k,l);
                        //msg(STATUS,"j,k,l: %f,%f,%f",pos[0],pos[1],pos[2]);
                        pCut(pop, s, pIndex, pop->pos, pop->vel);
                        cutNumber += 1;
                        //msg(STATUS,"iStop = %li",iStop);
                        iStop--;

                    }
                }

            }

        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &cutNumber, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    msg(STATUS,"cutNumber = %i ",cutNumber);
    cutNumber = 0;

    MPI_Allreduce(MPI_IN_PLACE, chargeCounter, obj->nObjects, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE, objectCurrent, nSpecies*obj->nObjects, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


    for (long int a=0; a<obj->nObjects; a++) {
      for (long int s=0; s<nSpecies; s++) {
        msg(STATUS,"current for objct %i, and species %i = %f",a,s,objectCurrent[a*nSpecies + s]);
      }
    }

    //MPI_Allreduce(MPI_IN_PLACE, invNrSurfNod, obj->nObjects, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    // Add the collected charge to the surface nodes on rhoObject.
    for (long int a=0; a<obj->nObjects; a++) {

	  //printf("chargeCounter[a] = %f\n",chargeCounter[a]);
      //printf("invNrSurfNod[a] = %f\n",invNrSurfNod[a]);
	  	//int testcounter = 0;
        for (long int b=lookupSurfOff[a]; b<lookupSurfOff[a+1]; b++) {
            val[obj->lookupSurface[b]] += chargeCounter[a]*invNrSurfNod[a];
			//testcounter += 1;
        }
		// MPI_Allreduce(MPI_IN_PLACE, &testcounter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		// printf("added to %i nodes, should be %f\n",testcounter,(1./invNrSurfNod[a]) );
    }
    //free(invNrSurfNod);
    free(chargeCounter);
    //free(nodCorLoc);
    //free(nodCorGlob);

}


//
// //stores index of particles that are close to an object at the current timestep
// void oVicinityParticles(Population *pop, PincObject *obj){
//
// 	double *val = obj->domain->val;
// 	int nSpecies = pop->nSpecies;
// 	int counter = 0;
// 	long int *sizeProd = obj->domain->sizeProd;
//
// 	for(int s=0; s < nSpecies; s++) {
//
// 		long int iStart = pop->iStart[s];
// 		long int iStop = pop->iStop[s];
//
// 		for(int i=iStart;i<iStop;i++){
//
// 			double *pos = &pop->pos[3*i];
//
// 			// Integer parts of position
// 			int j = (int) pos[0];
// 			int k = (int) pos[1];
// 			int l = (int) pos[2];
//
// 			// Index of nodes surrounding particle i
// 			long int p 		= j*3 + k*sizeProd[2] + l*sizeProd[3];
// 			long int pj 	= p + 3; //sizeProd[1];
// 			long int pk 	= p + sizeProd[2];
// 			long int pjk 	= pk + 3; //sizeProd[1];
// 			long int pl 	= p + sizeProd[3];
// 			long int pjl 	= pl + 3; //sizeProd[1];
// 			long int pkl 	= pl + sizeProd[2];
// 			long int pjkl 	= pkl + 3; //sizeProd[1];
//
// 			// All neighbours must be part of bounding box, value of obj->domain->val
// 			// at a vicinity node is set to 2, obj->domain->val set to 1 for object itself.
//             // When particle does not have an object node as a neighbour, this sum will 16
// 			int sum = val[p]+val[pj]+val[pk]+val[pjk]+val[pl]+val[pjl]+val[pkl]+val[pjkl];
//
// 			if(sum < 16 && sum > 0){
// 				pop->objVicinity[counter] = i;
// 				counter++;
// 			}
// 		}
// 	}
// }

//Relies on a courant number < 1 (otherwise particle might be inside object)
//checks which particles in object vicinity will collide => overwrites pop->collisions
// void oFindParticleCollisions(Population *pop, PincObject *obj){
//
//     long int *lookupIntOff = obj->lookupInteriorOffset;
//     long int *sizeProd = obj->domain->sizeProd;
//
//     oVicinityParticles(pop, obj);
//     long int *vicinity = pop->objVicinity;
//     long int nCloseParticles = sizeof(vicinity) / sizeof(vicinity[0]);
//     long int counter = 0;
//     alSetAll(pop->collisions, nCloseParticles, 0);
//
//
//     for(int i=0;i<nCloseParticles;i++){
//
//         double *pos = &pop->pos[3*i];
//         double *vel = &pop->vel[3*i];
//         double *nextPos;
//         adAdd(pos,vel,nextPos,3);
//
//         // Integer parts of position in next time step
//         int j = (int) nextPos[0];
//         int k = (int) nextPos[1];
//         int l = (int) nextPos[2];
//
//         long int p = j + k*sizeProd[2] + l*sizeProd[3];
//
//         // Check whether p is one of the object nodes
//         for (long int a=0; a<obj->nObjects; a++) {
//             for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
//                 if ((obj->lookupInterior[b])==p) {
//                     pop->collisions[counter] = p; // CHECK THIS: p is particle position, not index. The index is 3*i, or nDims*i
//                     counter++;
//                 }
//             }
//         }
//     }
// }

//Moves a particle according to the type of collision, also creates and removes new particles
// void oParticleCollision(Population *pop, PincObject *obj, long int i){
//
//     void (*collisionType)(Population *);
//
//     pFindCollisionType(pop, obj, i, collisionType);
//
//     //collisionType();
// }


//Finds nearest 3 object surface nodes to a specific particle of index p
//3 object surface nodes needed to compute normal from cross product of surface vectors
// double *oFindNearestSurfaceNodes(Population *pop, long int particleId, PincObject *obj){
//
//     double *pos = NULL;
//     for(int i=0; i<3; i++){
//         pos[i] = pop->pos[3*particleId + i];
//     }
//
//
//     return pos;
//
// }
//
// bool oParticleIntersection(Population *pop, long int particleId, PincObject *obj){
//
//     //find nearest nodes
//     double *nearest = oFindNearestSurfaceNodes(pop, particleId, obj);
//
//     nearest += 0;
//     return false;
// }
//
// //pos_new = pos_old + vel*delta_t
// //try http://geomalgorithms.com/a05-_intersect-1.html algorithm
// //implementation based on https://rosettacode.org/wiki/Find_the_intersection_of_a_line_with_a_plane#C
// void oFindIntersectPoint(const Population *pop, long int id, double *surfNormal,
//      double *surfPoint, double *intersect){
//
//         double epsilon = 1e-6;
//         double *pos = &pop->pos[3*id];
//         double *vel = &pop->vel[3*id];
//         double *w = NULL;
//         double *Psi = vel;
//         int ndotu = adDotProd(vel,surfNormal,3);
//
//         if(ndotu < epsilon){
//             msg(ERROR, "Particle %i, will not collide with any object next timestep!", id);
//         }
//         adSub(pos, surfPoint, w, 3);
//
//         //Compute intersection
//         double si = -1.*(double) adDotProd(surfNormal,w,3) / (double) ndotu;
//         adScale(Psi,3,si);
//         adAdd(w,Psi,Psi,3);
//         adAdd(surfPoint,Psi,Psi,3);
//
//         intersect[0]=Psi[0], intersect[1]=Psi[1], intersect[1]=Psi[1];
// }

//void oParticleCollision(Population *pop, PincObject *obj, long int i){

    //msg(WARNING, "Collision types not yet implemented!");
//}

/*****************************************************************************
 *  ALLOC/DESTRUCTORS
 ****************************************************************************/

PincObject *objoAlloc(const dictionary *ini, const MpiInfo *mpiInfo, Units *units){

    int size = mpiInfo->mpiSize;
    //int mpiRank = mpiInfo->mpiRank;
    Grid *domain = gAlloc(ini, SCALAR,mpiInfo);
    //int rank = domain->rank;
    gZero(domain);

    PincObject *obj = malloc(sizeof(*obj));
    obj->domain = domain;

    oOpenH5(ini, obj, mpiInfo, units, units->chargeDensity, "object");          // for capMatrix - objects
	oReadH5(obj);
    //oCloseH5(obj);
    //Communicate the boundary nodes
		gHaloOp(setSlice, obj->domain, mpiInfo, TOHALO);


    //obj->nObjects
    //obj->lookupInterior
    //obj->lookupInteriorOffset

    // Find the number of objects in the input file
    int nObjects = 0;
    for (int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>nObjects) {
            nObjects = (int)(obj->domain->val[i]+0.5); // Note, this is not necessarily
                //the number of objects, but rather the identifier of the object with the highest number.
                //Feel free to implement something more fancy here...
        }
    }
    // Make sure each process knows the total number of objects.
    MPI_Allreduce(MPI_IN_PLACE, &nObjects, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    //msg(WARNING|ALL,"nObjects: %i",nObjects);

    obj->nObjects = nObjects;

    oFillLookupTables(obj);

    oFindObjectSurfaceNodes(obj);

    long int *nodCorLoc = malloc((size+1)*sizeof(*nodCorLoc));
    long int *nodCorGlob = malloc(obj->nObjects*(size+1)*sizeof(*nodCorGlob));

	for (int i=0;i<size+1;i++){ // AD-HOC sol, initialize manually for safety.
		nodCorLoc[i] = 0;
	}

    double *capMatrixSum = malloc(obj->nObjects*sizeof(*capMatrixSum));
    //long int *capMatrixAllOffsets = malloc(obj->nObjects*(size+1)*sizeof(*capMatrixAllOffsets));

    long int capMatrixAllSize = oGatherSurfaceNodes(obj,nodCorLoc,nodCorGlob,obj->lookupSurfaceOffset,mpiInfo);

    double *capMatrixAll = malloc( (capMatrixAllSize*(capMatrixAllSize))*sizeof(*capMatrixAll) );

	//this is an unneccessary large array, because
	// we only evaluate one obj at a time.
	double *deltaPhi = malloc(capMatrixAllSize*sizeof(*deltaPhi));
	double *rhoCorr = malloc(capMatrixAllSize*sizeof(*rhoCorr));

	double *invNrSurfNod = malloc(obj->nObjects*sizeof(*invNrSurfNod));
    // for (long int a=0; a<obj->nObjects; a++) {
    //     invNrSurfNod[a] = 1.0/(nodCorGlob[(a+1)*(size)]);
    //     //printf("invNrSurfNod[a] = %f, nodCorGlob[(a+1)*(size+1)] = %li",invNrSurfNod[a],nodCorGlob[(a+1)*(size)]);
    // }



  int nSpecies = iniGetInt(ini,"population:nSpecies");
  double *objectCurrent= malloc(nSpecies*nObjects*sizeof(*objectCurrent));
  bool biasOn = iniGetInt(ini,"object:biasOn");
  double *bias = iniGetDoubleArr(ini,"object:bias",nObjects);

  adScale(bias, nObjects, 1./units->potential);
  //printf("bias=%f\n",bias[0] );
  //exit(0);
  obj->biasOn = biasOn;
  obj->bias = bias;
    obj->capMatrixAll = capMatrixAll;
    obj->capMatrixAllOffsets = nodCorGlob;
    obj->capMatrixSum = capMatrixSum;
	obj->deltaPhi = deltaPhi;
	obj->rhoCorr = rhoCorr;
	obj->invNrSurfNod = invNrSurfNod;
  obj->objectCurrent= objectCurrent;


    free(nodCorLoc);

    return obj;
}

void oFree(PincObject *obj){

    gFree(obj->domain);

    free(obj->lookupInterior);
    free(obj->lookupInteriorOffset);
    free(obj->lookupSurface);
    free(obj->lookupSurfaceOffset);
    free(obj->capMatrixAll);
    free(obj->capMatrixAllOffsets);
    free(obj->capMatrixSum);
	free(obj->rhoCorr);
	free(obj->deltaPhi);
	free(obj->invNrSurfNod);
  free(obj->objectCurrent);
  free(obj->bias);
    free(obj);

}

void oCloseH5(PincObject *obj){

    gCloseH5(obj->domain);
}

void oOpenH5(const dictionary *ini, PincObject *obj, const MpiInfo *mpiInfo,
             const Units *units, double denorm, const char *fName){

    gOpenH5(ini, obj->domain,   mpiInfo, units, denorm, fName);
}

void oReadH5(PincObject *obj){

    // Identical to gReadH5()
    hid_t fileSpace = obj->domain->h5FileSpace;
    hid_t memSpace = obj->domain->h5MemSpace;
    hid_t file = obj->domain->h5;
    double *val = obj->domain->val;

    hid_t pList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);

    char name[64];
    sprintf(name,"Object"); //Only line which is different from gReadH5().

    hid_t dataset = H5Dopen(file,name,H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);

    H5Dclose(dataset);
    H5Pclose(pList);


}

// void oReadH5(Grid *grid, const MpiInfo *mpiInfo, const char *name){
//
// 	hid_t fileSpace = grid->h5FileSpace;
// 	hid_t memSpace = grid->h5MemSpace;
// 	hid_t file = grid->h5;
// 	double *val = grid->val;
//
// 	// Enable collective datawriting
// 	hid_t pList = H5Pcreate(H5P_DATASET_XFER);
//     H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);
//
// 	//char name[64];
// 	//sprintf(name,name);
//
// 	hid_t dataset = H5Dopen(file,name,H5P_DEFAULT);
// 	H5Dread(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);
//
// 	H5Dclose(dataset);
// 	H5Pclose(pList);
//
// }


/******************************************************************************
 *  DEPRECIATED FUNCTION DEFINITIONS
 *****************************************************************************/
// Compute the capacitance matrix. (one big matrix containing all objects)
/* void oComputeCapacitanceMatrix_v1(PincObject *obj, const dictionary *ini, const MpiInfo *mpiInfo) {

    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurface = obj->lookupSurface;
    long int *lookupSurfaceOffset = obj->lookupSurfaceOffset;

    // Allocate and initialise the structures to run the potential solver.
    void (*solverInterface)() = select(ini, "methods:poisson", mgSolver_set, sSolver_set);
    void (*solve)() = NULL;
    void *(*solverAlloc)() = NULL;
    void (*solverFree)() = NULL;
    solverInterface(&solve, &solverAlloc, &solverFree);

    Grid *rho = gAlloc(ini, SCALAR);
    Grid *phi = gAlloc(ini, SCALAR);

    void *solver = solverAlloc(ini, rho, phi);

    // Set Rho to zero.
    gZero(rho);

    // Find the number of surface nodes for each object.
    long int *nodesCoreLocal = malloc((size+1)*sizeof(*nodesCoreLocal));
    long int *nodesCoreGlobal = malloc(obj->nObjects*(size+1)*sizeof(*nodesCoreGlobal));
    long int nodesThisCore;
    for (long int a=0; a<obj->nObjects; a++) {

        nodesThisCore = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];

        // Let every core know how many surface nodes everybody has.
        MPI_Allgather(&nodesThisCore, 1, MPI_LONG, nodesCoreLocal, 1, MPI_LONG, MPI_COMM_WORLD);

        for(long int i=size-1;i>-1;i--) nodesCoreLocal[i+1]=nodesCoreLocal[i];
        nodesCoreLocal[0] = 0;
        alCumSum(nodesCoreLocal+1,nodesCoreLocal,size);

        for (long int b=0; b<size+1; b++) nodesCoreGlobal[a*(size+1)+b] = nodesCoreLocal[b];
    }

    // Find the size and initialise the array holding the capacitance matrix.
    long int capMatrixSize = 0;
    for (long int a=0; a<obj->nObjects; a++) {
        capMatrixSize +=nodesCoreGlobal[a*(size+1)+size];
    }

    // Initialise the capacitance matrix and its inverse.
    double *capMatrix = malloc( (capMatrixSize*capMatrixSize) * sizeof(*capMatrix));
    double *invCapMatrix = malloc( (capMatrixSize*capMatrixSize) * sizeof(*capMatrix));
    adSetAll(capMatrix,capMatrixSize*capMatrixSize,0);
    adSetAll(invCapMatrix,capMatrixSize*capMatrixSize,0);

    // Compute the actual capacitance matrix. (one big matrix containing all objects)
    long int meuh;
    long int moo = 0;
    long int *boo = malloc( (size+1) * sizeof(*boo));

    boo[0] = 0;
    MPI_Allgather(&(lookupSurfaceOffset[obj->nObjects]), 1, MPI_LONG, boo+1, 1, MPI_LONG, MPI_COMM_WORLD);
    alCumSum(boo+1,boo,size);

    for (long int r=0; r<size; r++) {
        if (r==rank) {
            msg(STATUS|ALL, "Computing the capacitance matrix components on core %i out of %i.", r+1, size);

            meuh = lookupSurfaceOffset[obj->nObjects];
        }
        MPI_Bcast(&meuh, 1, MPI_LONG, r, MPI_COMM_WORLD);

        for ( long int i=0; i<meuh; i++) {
            if (r==rank) {
                // Set the surface node charge to 1.
                rho->val[lookupSurface[i]] = 1;
            }
            // Solve for the potential.
            solve(solver, rho, phi, mpiInfo);

            if (r==rank) {
                // Set the surface node back to zero.
                rho->val[lookupSurface[i]] = 0;
            }

            // Fill column i of the capacitance matrix.
            for (long int k=0; k<lookupSurfaceOffset[obj->nObjects]; k++) {
                capMatrix[(boo[rank]+k) * capMatrixSize + moo + i] = phi->val[lookupSurface[k]];
            }
        }
        moo += meuh;
    }

    // Make sure every core has the complete matrix (needed for BLAS).
    MPI_Allreduce(MPI_IN_PLACE, capMatrix, (capMatrixSize*capMatrixSize), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (rank==0) {
        for (int lll=0; lll<3; lll++){
            adPrint(&capMatrix[lll*capMatrixSize],capMatrixSize);
        }
    }

    // Compute the inverse of the capacitance matrix.
    // Actually, the inverse is the capacitance matrix. Probably have to rethink the variable names.
    gsl_matrix_view A = gsl_matrix_view_array(capMatrix, capMatrixSize, capMatrixSize);
    gsl_matrix_view invA = gsl_matrix_view_array(invCapMatrix, capMatrixSize, capMatrixSize);

    int s;
    gsl_permutation *p = gsl_permutation_alloc(capMatrixSize);
    gsl_linalg_LU_decomp(&A.matrix, p, &s);
    gsl_linalg_LU_invert(&A.matrix, p, &invA.matrix);


    // Compute here the inverse total sum of elements in the capacitance matrix (needed later).
    double capMatrixInvSum = 1/adSum(invCapMatrix,capMatrixSize*capMatrixSize);

    // Add to object
    obj->capMatrix = invCapMatrix;
    obj->capMatrixInvSum = capMatrixInvSum;

    //alPrint(&capMatrixSize,1);
    //for (int i=0; i<capMatrixSize; i++) {
    //    for (int j=0; j<capMatrixSize; j++) {
    //        if (i>217 && j>217) {
    //            invCapMatrix[i*capMatrixSize+j] = 1;
    //        } else {
    //            invCapMatrix[i*capMatrixSize+j] = 0;
    //        }
    //    }
    //}

    // Need to compute here the bits needed for the mutual impedance stuff...
    // Remember, every core has the complete capacitance matrix.

    double *needCoffeeMatrix = malloc( (obj->nObjects*obj->nObjects) * sizeof(*needCoffeeMatrix));
    double *invNeedCoffeeMatrix = malloc( (obj->nObjects*obj->nObjects) * sizeof(*invNeedCoffeeMatrix));
    adSetAll(needCoffeeMatrix,obj->nObjects*obj->nObjects,0);
    adSetAll(invNeedCoffeeMatrix,obj->nObjects*obj->nObjects,0);

    long int totalNodesObject;
    long int totalNodesObjectCumA = 0;
    long int totalNodesObjectCumB = 0;

    for (int a = 0; a<obj->nObjects; a++) {
        totalNodesObject = nodesCoreGlobal[a*(size+1)+size];
        nodesThisCore = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];
        msg(WARNING|ALL,"Object %i, total %i, this core %i", a, totalNodesObject, nodesThisCore);

        totalNodesObjectCumB = 0;
        for (int b = 0; b<obj->nObjects; b++) {
            for (long int i=0;i<totalNodesObject;i++) {
                for (long int j=0;j<nodesThisCore;j++) {
                    needCoffeeMatrix[a*obj->nObjects + b] += invCapMatrix[(capMatrixSize*totalNodesObjectCumA) + j*capMatrixSize + totalNodesObjectCumB + i];
                }
            }
            totalNodesObjectCumB+=totalNodesObject;
        }
        totalNodesObjectCumA+=totalNodesObject;
    }

    //for (int lll=0; lll<obj->nObjects; lll++) {
    //    adPrint(&needCoffeeMatrix[lll*obj->nObjects],obj->nObjects);
    //}

    // Make sure every core has the complete matrix (needed for BLAS).
    MPI_Allreduce(MPI_IN_PLACE, needCoffeeMatrix, (obj->nObjects*obj->nObjects), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int lll=0; lll<obj->nObjects; lll++) {
        adPrint(&needCoffeeMatrix[lll*obj->nObjects],obj->nObjects);
    }

    gsl_matrix_view B = gsl_matrix_view_array(needCoffeeMatrix, obj->nObjects, obj->nObjects);
    gsl_matrix_view invB = gsl_matrix_view_array(invNeedCoffeeMatrix, obj->nObjects, obj->nObjects);

    int t;
    gsl_permutation *q = gsl_permutation_alloc(obj->nObjects);
    gsl_linalg_LU_decomp(&B.matrix, q, &t);
    gsl_linalg_LU_invert(&B.matrix, q, &invB.matrix);


    // Add to object
    obj->invNeedCoffeeMatrix = invNeedCoffeeMatrix;

    solverFree(solver);
} */

//Find all the object nodes which are part of the object surface.
/* void oFindObjectSurfaceNodes(PincObject *obj, const MpiInfo *mpiInfo) {

    long int *sizeProd = obj->domain->sizeProd;
    double *val = obj->domain->val;


    // Initialise the array storing the offsets for the surface nodes in the lookup table.
    long int *lookupSurfOff = malloc((obj->nObjects+1)*sizeof(*lookupSurfOff));
    alSetAll(lookupSurfOff,obj->nObjects+1,0);

    // Find the 8 neighbour cells of each non-ghost node.
    long int *myNB = malloc(9*sizeof(*myNB));
    // Find the ofsetts first.
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<sizeProd[obj->domain->rank]; b++) {
            if (!oIsGhostNode(obj->domain, b)) {
                myNB[0] = b;    // me on node i,j,k
                myNB[1] = myNB[0];                  // cell i,j,k
                myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
                myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
                myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
                myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];      // cell i-1,j-1,k
                myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];   // cell i-1,j-1,k-1

                int d=0;
                if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
                if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
                if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
                if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
                if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
                if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
                if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
                if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;

                // Check if on surface.
                if (d<7.5 && d>0) {
                    lookupSurfOff[a+1]++;
                }
            }
        }
    }
    alCumSum(lookupSurfOff+1,lookupSurfOff,obj->nObjects);

    // Initialise and compute the lookup table.
    long int *lookupSurf = malloc((lookupSurfOff[obj->nObjects])*sizeof(*lookupSurf));
    alSetAll(lookupSurf,lookupSurfOff[obj->nObjects]+1,0);

    long int *index = malloc((obj->nObjects)*sizeof(*index));
    for (long int i=0; i<obj->nObjects; i++) {
        index[i]=lookupSurfOff[i];
    }

    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
            if (!oIsGhostNode(obj->domain, b)) {
                myNB[0] = b;    // me on node i,j,k
                myNB[1] = myNB[0];                  // cell i,j,k
                myNB[2] = myNB[0] - sizeProd[3];    // cell i,j,k-1
                myNB[3] = myNB[0] - sizeProd[1];                    // cell i-1,j,k
                myNB[4] = myNB[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNB[5] = myNB[0] - sizeProd[2];                    // cell i,j-1,k
                myNB[6] = myNB[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNB[7] = myNB[0] - sizeProd[2] - sizeProd[1];                  // cell i-1,j-1,k;
                myNB[8] = myNB[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];    // cell i-1,j-1,k-1

                int d=0;
                if (val[myNB[1]]>(a+0.5) && val[myNB[1]]<(a+1.5)) d++;
                if (val[myNB[2]]>(a+0.5) && val[myNB[2]]<(a+1.5)) d++;
                if (val[myNB[3]]>(a+0.5) && val[myNB[3]]<(a+1.5)) d++;
                if (val[myNB[4]]>(a+0.5) && val[myNB[4]]<(a+1.5)) d++;
                if (val[myNB[5]]>(a+0.5) && val[myNB[5]]<(a+1.5)) d++;
                if (val[myNB[6]]>(a+0.5) && val[myNB[6]]<(a+1.5)) d++;
                if (val[myNB[7]]>(a+0.5) && val[myNB[7]]<(a+1.5)) d++;
                if (val[myNB[8]]>(a+0.5) && val[myNB[8]]<(a+1.5)) d++;

                // Check if on surface.
                if (d<7.5 && d>0) {
                    lookupSurf[index[a]] = myNB[0];
                    index[a]++;
                }
            }
        }
    }

    // Add to object.
    obj->lookupSurface = lookupSurf;
    obj->lookupSurfaceOffset = lookupSurfOff;
}
 */
//
// //depreciated...
// void oFindObjectSurfaceNodes_v1(PincObject *obj) {
//
//     long int *sizeProd = obj->domain->sizeProd;
//
//     //if(oIsGhostNode(obj->domain, 10439)) printf("test\n");
//
//     long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*\//
//                                            sizeof(*lookupSurfaceOffset));
//     for (long int i=0; i<obj->nObjects+1; i++) {
//         lookupSurfaceOffset[i] = 0;
//     }
//
//     long int *meNeighbours = malloc(7*sizeof(*meNeighbours));
//     for (long int a=0; a<obj->nObjects; a++) {
//         for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
//             if (obj->domain->val[b]>0.5 && !oIsGhostNode(obj->domain, b)) {
//                 // My neighbours
//                 meNeighbours[0] = b; // me
//                 meNeighbours[1] = meNeighbours[0] + sizeProd[1]; // right
//                 meNeighbours[2] = meNeighbours[0] - sizeProd[1]; // left
//                 meNeighbours[3] = meNeighbours[0] + sizeProd[2]; // up
//                 meNeighbours[4] = meNeighbours[0] - sizeProd[2]; // down
//                 meNeighbours[5] = meNeighbours[0] + sizeProd[3]; // back
//                 meNeighbours[6] = meNeighbours[0] - sizeProd[3]; // front
//
//                 int d=0;
//                 if (obj->domain->val[meNeighbours[1]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[2]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[3]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[4]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[5]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[6]]>(a+0.5)) d++;
//
//                 // Check if surface
//                 if (d<5.5) {
//                     lookupSurfaceOffset[a+1]++;
//                 }
//             }
//         }
//         //MPI_Allreduce(MPI_IN_PLACE, &lookupSurfaceOffset[a+1], 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//     }
//
//     alCumSum(lookupSurfaceOffset+1,lookupSurfaceOffset,obj->nObjects);
//
//
//     //Second go through to fill the table.
//     long int *lookupSurface = malloc((lookupSurfaceOffset[obj->nObjects])*\//
//                                      sizeof(*lookupSurface));
//     for (long int i=0; i<lookupSurfaceOffset[obj->nObjects]+1; i++) {
//         lookupSurface[i]=0;
//     }
//
//     long int *index = malloc((obj->nObjects)*sizeof(*index));
//     for (long int i=0; i<obj->nObjects; i++) {
//         index[i]=lookupSurfaceOffset[i];
//     }
//
//     for (long int a=0; a<obj->nObjects; a++) {
//         for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
//             if (obj->domain->val[b]>0.5 && !oIsGhostNode(obj->domain, b)) {
//                 // My neighbours
//                 meNeighbours[0] = b; // me
//                 meNeighbours[1] = meNeighbours[0] + sizeProd[1]; // right
//                 meNeighbours[2] = meNeighbours[0] - sizeProd[1]; // left
//                 meNeighbours[3] = meNeighbours[0] + sizeProd[2]; // up
//                 meNeighbours[4] = meNeighbours[0] - sizeProd[2]; // down
//                 meNeighbours[5] = meNeighbours[0] + sizeProd[3]; // back
//                 meNeighbours[6] = meNeighbours[0] - sizeProd[3]; // front
//
//                 int d=0;
//                 if (obj->domain->val[meNeighbours[1]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[2]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[3]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[4]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[5]]>(a+0.5)) d++;
//                 if (obj->domain->val[meNeighbours[6]]>(a+0.5)) d++;
//
//                 // Check if boundary
//                 if (d<5.5) {
//                     lookupSurface[index[a]] = meNeighbours[0];
//                     index[a]++;
//                 }
//             }
//         }
//     }
//
//     // Add to object
//     obj->lookupSurface = lookupSurface;
//     obj->lookupSurfaceOffset = lookupSurfaceOffset;
// }

// Construct and solve equation 5 in Miyake_Usui_PoP_2009
/* void oApplyCapacitanceMatrixoApplyCapacitanceMatrix_v1(Grid *rho, const Grid *phi, const PincObject *obj, const MpiInfo *mpiInfo){

    //int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurface = obj->lookupSurface;
    long int *lookupSurfaceOffset = obj->lookupSurfaceOffset;
    double *capMatrix = obj->capMatrix;
    double *invNeedCoffeeMatrix = obj->invNeedCoffeeMatrix;

    // Compute the righthand components for the needCoffee bit.
    double *needCoffeeRight = malloc( (obj->nObjects) * sizeof(*needCoffeeRight));
    adSetAll(needCoffeeRight,obj->nObjects,0);
    double *dummy = malloc( (obj->nObjects) * sizeof(*dummy));
    adSetAll(dummy,obj->nObjects,0);



    // Find the number of surface nodes for each object.
    long int *nodesCoreLocal = malloc((size+1)*sizeof(*nodesCoreLocal));
    long int *nodesCoreGlobal = malloc(obj->nObjects*(size+1)*sizeof(*nodesCoreGlobal));
    long int nodesThisCore;

    for (long int a=0; a<obj->nObjects; a++) {

        nodesThisCore = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];

        // Let every core know how many surface nodes everybody has.
        MPI_Allgather(&nodesThisCore, 1, MPI_LONG, nodesCoreLocal, 1, MPI_LONG, MPI_COMM_WORLD);

        for(long int i=size-1;i>-1;i--) nodesCoreLocal[i+1]=nodesCoreLocal[i];
        nodesCoreLocal[0] = 0;
        alCumSum(nodesCoreLocal+1,nodesCoreLocal,size);

        for (long int b=0; b<size+1; b++) nodesCoreGlobal[a*(size+1)+b] = nodesCoreLocal[b];
    }

    // Find the size and initialise the array holding the capacitance matrix.
    long int capMatrixSize = 0;
    for (long int a=0; a<obj->nObjects; a++) {
        capMatrixSize +=nodesCoreGlobal[a*(size+1)+size];
    }



    long int totalNodesObject = 0;
    long int totalNodesObjectCumA = 0;
    long int totalNodesObjectCumB = 0;

    for (int a = 0; a<obj->nObjects; a++) {
        totalNodesObject = nodesCoreGlobal[a*(size+1)+size];
        nodesThisCore = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];

        totalNodesObjectCumB = 0;
        for (int b = 0; b<obj->nObjects; b++) {

            for (long int i=0;i<nodesCoreGlobal[b*(size+1)+size];i++) {
                for (long int j=0;j<(lookupSurfaceOffset[b+1] - lookupSurfaceOffset[b]);j++) {

                    // No rank needed here, that's suspicious...
                    dummy[b] += capMatrix[(capMatrixSize*totalNodesObjectCumA) +j*capMatrixSize + totalNodesObjectCumB + i] * phi->val[lookupSurface[lookupSurfaceOffset[b] + j]];

                }
            }
            totalNodesObjectCumB+=totalNodesObject;
        }
        totalNodesObjectCumA+=totalNodesObject;
        for (int b = 0; b<obj->nObjects; b++) {
            needCoffeeRight[a] += dummy[b];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, needCoffeeRight, (obj->nObjects), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Now compute the potentials on each object.
    double *phiC = malloc( (obj->nObjects) * sizeof(*phiC));
    adSetAll(phiC,obj->nObjects,0);

    for (int i=0; i<obj->nObjects; i++) {
        for (int j=0; j<obj->nObjects; j++) {
            phiC[i] = invNeedCoffeeMatrix[i*obj->nObjects+j]*needCoffeeRight[j];
        }
    }

    // for testing purposes:
    for (int i=0; i<obj->nObjects; i++) { phiC[i] = 0.5; }

    // compute deltaPhi
    // compute eq. 5
    // add charge corrections.



} */



static void oMode(dictionary *ini){

	/*
	 * SELECT METHODS
	 */
	void (*acc)()   			= select(ini,	"methods:acc",
												puAcc3D1_set,
												puAcc3D1KE_set,
												puAccND1_set,
												puAccND1KE_set,
												puAccND0_set,
												puAccND0KE_set,
                        puBoris3D1KETEST_set);

	void (*distr)() 			= select(ini,	"methods:distr",
												puDistr3D1split_set,
												puDistr3D1_set,
												puDistrND1_set,
												puDistrND0_set);

	void (*extractEmigrants)()	= select(ini,	"methods:migrate",
												puExtractEmigrants3D_set,
												puExtractEmigrantsND_set,
                        						puExtractEmigrants3DOpen_set);

	void (*solverInterface)()	= select(ini,	"methods:poisson",
												mgSolver_set,
												sSolver_set);

	void (*solve)() = NULL;
	void *(*solverAlloc)() = NULL;
	void (*solverFree)() = NULL;
	solverInterface(&solve, &solverAlloc, &solverFree);

	/*
	 * INITIALIZE PINC VARIABLES
	 */
	Units *units=uAlloc(ini);
	uNormalize(ini, units);

	MpiInfo *mpiInfo = gAllocMpi(ini);
	Population *pop = pAlloc(ini,mpiInfo);
	Grid *E   = gAlloc(ini, VECTOR,mpiInfo);
	Grid *rho = gAlloc(ini, SCALAR,mpiInfo);
	Grid *rho_e = gAlloc(ini, SCALAR, mpiInfo);
	Grid *rho_i = gAlloc(ini, SCALAR, mpiInfo);
    Grid *rhoObj = gAlloc(ini, SCALAR,mpiInfo);     // for capMatrix - objects
	Grid *phi = gAlloc(ini, SCALAR,mpiInfo);



	void *solver = solverAlloc(ini, rho, phi, mpiInfo);

    PincObject *obj = objoAlloc(ini,mpiInfo,units);              // for capMatrix - objects
//TODO: look into multigrid E,rho,rhoObj

	// Creating a neighbourhood in the rho to handle migrants
	gCreateNeighborhood(ini, mpiInfo, rho);

  	// Setting Boundary slices
  	gSetBndSlices(ini, phi, mpiInfo);
	//gSetBndSlices(ini, rho, mpiInfo);
	gSetBndSlicesE(ini, E, mpiInfo);

	// Random number seeds
	gsl_rng *rngSync = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng,mpiInfo->mpiRank+1); // Seed needs to be >=1

	/*
	 * PREPARE FILES FOR WRITING
	 */

	pOpenH5(ini, pop, units, "pop");
	//double denorm = units->potential;
	gOpenH5(ini, rho, mpiInfo, units, units->chargeDensity, "rho");
	gOpenH5(ini, rho_e, mpiInfo, units, units->chargeDensity, "rho_e");
	gOpenH5(ini, rho_i, mpiInfo, units, units->chargeDensity, "rho_i");
	gOpenH5(ini, phi, mpiInfo, units, units->potential, "phi");
	gOpenH5(ini, E,   mpiInfo, units, units->eField, "E");
  // oOpenH5(ini, obj, mpiInfo, units, 1, "test");
  // oReadH5(obj, mpiInfo);


    //msg(STATUS,"opening obj file");
		gOpenH5(ini, rhoObj, mpiInfo, units, units->chargeDensity, "rhoObj");        // for capMatrix - objects
		//oOpenH5(ini, obj, mpiInfo, units, units->chargeDensity, "object");          // for capMatrix - objects
		//oReadH5(obj->domain, mpiInfo, "Object");

    //msg(STATUS,"done");


		//Count the number of objects and fill the lookup tables.
		// This is done in oAlloc now....

		//oFillLookupTables(obj,mpiInfo);
		// Find all the object nodes which are part of the object surface.
		//oFindObjectSurfaceNodes(obj, mpiInfo);


	hid_t history = xyOpenH5(ini,"history");
	pCreateEnergyDatasets(history,pop);
  xyCreateDataset(history,"/current/electrons/dataset");
  xyCreateDataset(history,"/current/ions/dataset");

	// Add more time series to history if you want
	// xyCreateDataset(history,"/group/group/dataset");

	/*
	 * INITIAL CONDITIONS
	 */

    //Compute capacitance matrix
    oComputeCapacitanceMatrix(obj, ini, mpiInfo);

	// Initalize particles
	//pPosUniform(ini, pop, mpiInfo, rngSync);
	//pPosLattice(ini, pop, mpiInfo);
	pPosUniformCell(ini,rho,pop,rng);
	//pVelZero(pop);
	//pVelMaxwell(ini, pop, rng);
	double maxVel = iniGetDouble(ini,"population:maxVel");




	// Perturb particles
	//pPosPerturb(ini, pop, mpiInfo);


	//add influx of new particles on boundary
	pPurgeGhost(pop, rho);

	// Migrate those out-of-bounds due to perturbation
	extractEmigrants(pop, mpiInfo);
	puMigrate(pop, mpiInfo, rho);

	pFillGhost(ini,rho,pop,rng);




	/*
	 * INITIALIZATION (E.g. half-step)
	 */

    // Clean objects from any charge first.
    gZero(rhoObj);                                          // for capMatrix - objects
    oCollectObjectCharge(pop, rhoObj, obj, mpiInfo);        // for capMatrix - objects
    gZero(rhoObj);                                          // for capMatrix - objects


	// Get initial charge density
	distr(pop, rho,rho_e,rho_i);
	gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
	gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);
    //gWriteH5(rho, mpiInfo, (double) 0);

	// Get initial E-field


  //gBnd(phi, mpiInfo);
	solve(solver, rho, phi, mpiInfo);
	//gNeutralizeGrid(phi, mpiInfo);
	//gBnd(phi, mpiInfo);
    //gWriteH5(phi, mpiInfo, (double) 0);
    //pWriteH5(pop, mpiInfo, (double) 0, (double)0+0.5);

	gFinDiff1st(phi, E);
	gHaloOp(setSlice, E, mpiInfo, TOHALO);
	gMul(E, -1.);
	gBnd(E, mpiInfo);


  //Boris parameters
  int nSpecies = pop->nSpecies;
	double *S = (double*)malloc((3)*(nSpecies)*sizeof(double));
	double *T = (double*)malloc((3)*(nSpecies)*sizeof(double));

  // add External E
	//gZero(E); // for testing Boris
	//gAddTo(Ext); //needs grid definition of Eext
  	puAddEext(ini, pop, E); // adds same value to whole grid


  	gMul(E, 0.5);
	puGet3DRotationParameters(ini, T, S, 0.5);
	acc(pop, E, T, S);
	gMul(E, 2.0);
	puGet3DRotationParameters(ini, T, S, 1.0);

	/*
	 * TIME LOOP
	 */

	Timer *t = tAlloc(mpiInfo->mpiRank);

	// n should start at 1 since that's the timestep we have after the first
	// iteration (i.e. when storing H5-files).
	int nTimeSteps = iniGetInt(ini,"time:nTimeSteps");
	for(int n = 1; n <= nTimeSteps; n++){

		long int totPs0 = (pop->iStop[0]- pop->iStart[0]); //debug
		long int totPs1 = (pop->iStop[1]- pop->iStart[1]);
		MPI_Allreduce(MPI_IN_PLACE, &totPs0, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &totPs1, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		msg(STATUS,"Computing time-step %i",n);
        msg(STATUS, "Nr. of particles s=0 %i: ",totPs0);
		msg(STATUS, "Nr. of particles s=1 %i: ",totPs1);



		// Check that no particle moves beyond a cell (mostly for debugging)
		pVelAssertMax(pop,maxVel);

		tStart(t);

		// Move particles
		// oRayTrace(pop, obj, deltaRho); <- do we need this still???

		puMove(pop); //puMove(pop, obj); Do not change functions such that PINC does
    // not work in other run modes!

		//add influx of new particles on boundary
		pPurgeGhost(pop, rho);

		// Migrate particles (periodic boundaries)
		extractEmigrants(pop, mpiInfo);
		puMigrate(pop, mpiInfo, rho);

		pFillGhost(ini,rho,pop,rng);


		// Check that no particle resides out-of-bounds (just for debugging)
		//pPosAssertInLocalFrame(pop, rho); //gives error with open boundary

        // Collect the charges on the objects.
        oCollectObjectCharge(pop, rhoObj, obj, mpiInfo);    // for capMatrix - objects


		// Compute charge density
		distr(pop, rho,rho_e,rho_i);
		gHaloOp(addSlice, rho, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_e, mpiInfo, FROMHALO);
		gHaloOp(addSlice, rho_i, mpiInfo, FROMHALO);


        // Keep writing Rho here.

        // Add object charge to rho.
		gAddTo(rho, rhoObj);

        //gBnd(phi, mpiInfo);
        solve(solver, rho, phi, mpiInfo);                   // for capMatrix - objects
		//gNeutralizeGrid(phi, mpiInfo);
		//gBnd(phi, mpiInfo);
        // Second run with solver to account for charges
        oApplyCapacitanceMatrix(rho, phi, obj, mpiInfo, units);    // for capMatrix - objects

		//gBnd(phi, mpiInfo);
		solve(solver, rho, phi, mpiInfo);
		//gNeutralizeGrid(phi, mpiInfo);
		//gBnd(phi, mpiInfo);
		//gHaloOp(setSlice, phi, mpiInfo, TOHALO); // Needed by sSolve but not mgSolve

		double rhoSum = gSumTruegrid(rho);
		double rhoObjSum = gSumTruegrid(rhoObj);
		MPI_Allreduce(MPI_IN_PLACE, &rhoObjSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &rhoSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		msg(STATUS,"total charge = %f PINC values",(rhoSum));
		msg(STATUS,"Object charge = %.6e C",(units->charge*rhoObjSum));
		// Compute E-field
		gFinDiff1st(phi, E);
		gHaloOp(setSlice, E, mpiInfo, TOHALO);
		gMul(E, -1.);
		gBnd(E, mpiInfo);
		//gBndE(E, mpiInfo); // always neumann cond

		//gAssertNeutralGrid(E, mpiInfo);
		// Apply external E
		//gZero(E);
		//gAddTo(Ext); //needs grid definition of Eext
		puAddEext(ini, pop, E); // adds same value to whole grid

		// Accelerate particle and compute kinetic energy for step n
		//acc(pop, E);
		acc(pop, E, T, S);

		tStop(t);

		// Sum energy for all species
		pSumKinEnergy(pop);

		// Compute potential energy for step n
		gPotEnergy(rho,phi,pop);

		// Example of writing another dataset to history.xy.h5
		// xyWrite(history,"/group/group/dataset",(double)n,value,MPI_SUM);

		if(n%1000 == 0 || n>29500){//n>122700){//50614
		//Write h5 files
			//gWriteH5(E, mpiInfo, (double) n);
			gWriteH5(rho, mpiInfo, (double) n);
			gWriteH5(rho_e, mpiInfo, (double) n);
			gWriteH5(rho_i, mpiInfo, (double) n);

			gWriteH5(phi, mpiInfo, (double) n);
			//pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
			//gWriteH5(rhoObj, mpiInfo, (double) n);
		}
		// if(n%1 == 0){
		// 	pWriteH5(pop, mpiInfo, (double) n, (double)n+0.5);
		// }

		pWriteEnergy(history,pop,(double)n,units);
    xyWrite(history,"/current/electrons/dataset",(double)n,units->current*obj->objectCurrent[0],MPI_SUM);
    xyWrite(history,"/current/ions/dataset",(double)n,units->current*obj->objectCurrent[1],MPI_SUM);
	}

	if(mpiInfo->mpiRank==0) {
    tMsg(t->total, "Time spent: ");
}

	/*
	 * FINALIZE PINC VARIABLES
	 */
	gFreeMpi(mpiInfo);


	// Close h5 files
	pCloseH5(pop);
	gCloseH5(rho);
	gCloseH5(rho_e);
	gCloseH5(rho_i);

	gCloseH5(phi);
	gCloseH5(E);
    gCloseH5(rhoObj);       // for capMatrix - objects
    oCloseH5(obj);          // for capMatrix - objects
    // 11.10.19 segfault seems to link to oClose(), as calling this
    // alters the segfault.

	xyCloseH5(history);

  // Free memory
  // sFree(solver);
  // mgFreeSolver(solver);
  solverFree(solver);
  gFree(rho);
  gFree(rho_e);
  gFree(rho_i);
  gFree(phi);
  free(S);
  free(T);

  gFree(E);
  pFree(pop);
  uFree(units);
    gFree(rhoObj);          // for capMatrix - objects
    oFree(obj);             // for capMatrix - objects




	gsl_rng_free(rngSync);
	gsl_rng_free(rng);


}

funPtr oMode_set(){ // dictionary *ini
	  // TODO: sanity
	return oMode;
}
