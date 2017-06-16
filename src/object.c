/**
 * @file		object.c
 * @author		Jan Deca <jandeca@gmail.com>
 * @brief		All object-related functions are here.
 * @date		19.10.16
 */

#include "core.h"
#include "object.h"
#include "multigrid.h"
#include "spectral.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdbool.h>

/******************************************************************************
 *  LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/

/**
 * @brief   Count the number of objects and fills the lookup tables.
 * @param	obj		Object
 * @param	ini		input settings
 * @return	void
 */
void oFillLookupTables(Object *obj, const MpiInfo *mpiInfo);

/**
 * @brief   Check whether a certain node is a ghost node.
 * @param	grid	Grid
 * @param	node	long int
 * @return	bool
 */
bool isGhostNode(Grid *grid, long int node);
void oGhost(long int node, const int *nGhostLayersBefore,
            const int *nGhostLayersAfter, const int *trueSize,
            const long int *sizeProd, bool *ghost);

/******************************************************************************
 *  LOCAL FUNCTION DEFINITIONS
 *****************************************************************************/

//Check whether a certain node is a ghost node.
bool isGhostNode(Grid *grid, long int node) {
    
    long int *sizeProd = grid->sizeProd;
    int *trueSize = grid->trueSize;
    int *nGhostLayers = grid->nGhostLayers;
    int rank = grid->rank;
    
    bool ghost = false;
    
    oGhost(node,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1],&ghost);
    
    return ghost;
}
void oGhost(long int node, const int *nGhostLayersBefore,
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

// Count the number of objects and fill the lookup tables
void oFillLookupTables(Object *obj, const MpiInfo *mpiInfo) {
    
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
    
    
    // Initialise and compute the array that stores the offsets of the objects in the lookup table.
    long int *lookupInteriorOffset = malloc((nObjects+1)*sizeof(*lookupInteriorOffset));
    for (long int i=0; i<nObjects+1; i++) {
        lookupInteriorOffset[i] = 0;
    }
   
    for (long int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>0.5){
            lookupInteriorOffset[(int)(obj->domain->val[i]+0.5)]++;
        }
    }
    alCumSum(lookupInteriorOffset+1,lookupInteriorOffset,nObjects);

    // Initialise and compute the lookup table.
    long int *lookupInterior = malloc((lookupInteriorOffset[nObjects])*sizeof(*lookupInterior));
    for (long int i=0; i<lookupInteriorOffset[nObjects]; i++) {
        lookupInterior[i]=0;
    }
   
    long int *index = malloc((nObjects)*sizeof(*index));
    for (long int i=0; i<nObjects; i++) {
        index[i]=lookupInteriorOffset[i];
    }
   
    for (long int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>0.5 && !isGhostNode(obj->domain, i)){
            lookupInterior[(index[(int)(obj->domain->val[i]+0.5)-1])] = i;
            (index[(int)(obj->domain->val[i]+0.5)-1])++;
        }
    }
    
    // Add to the obj Object.
    obj->nObjects = nObjects;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;
}

//Find all the object nodes which are part of the object surface.
void oFindObjectSurfaceNodes(Object *obj, const MpiInfo *mpiInfo) {
    
    long int *sizeProd = obj->domain->sizeProd;
    double *val = obj->domain->val;
    
    // Initialise the array storing the offsets for the surface nodes in the lookup table.
    long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*sizeof(*lookupSurfaceOffset));
    alSetAll(lookupSurfaceOffset,obj->nObjects+1,0);
    
    // Find the 8 neighbour cells of each non-ghost node.
    long int *myNB = malloc(9*sizeof(*myNB));
    // Find the ofsetts first.
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<sizeProd[obj->domain->rank]; b++) {
            if (!isGhostNode(obj->domain, b)) {
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
                    lookupSurfaceOffset[a+1]++;
                }
            }
        }
    }
    alCumSum(lookupSurfaceOffset+1,lookupSurfaceOffset,obj->nObjects);
    
    // Initialise and compute the lookup table.
    long int *lookupSurface = malloc((lookupSurfaceOffset[obj->nObjects])*sizeof(*lookupSurface));
    alSetAll(lookupSurface,lookupSurfaceOffset[obj->nObjects]+1,0);
    
    long int *index = malloc((obj->nObjects)*sizeof(*index));
    for (long int i=0; i<obj->nObjects; i++) {
        index[i]=lookupSurfaceOffset[i];
    }
    
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
            if (!isGhostNode(obj->domain, b)) {
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
                    lookupSurface[index[a]] = myNB[0];
                    index[a]++;
                }
            }
        }
    }
    
    // Add to object.
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;
}

// Compute the capacitance matrix. (one big matrix containing all objects)
void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini, const MpiInfo *mpiInfo) {
    
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
    MPI_Barrier(MPI_COMM_WORLD);
    
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
}
// Construct and solve equation 5 in Miyake_Usui_PoP_2009
void oApplyCapacitanceMatrix(Grid *rho, const Grid *phi, const Object *obj, const MpiInfo *mpiInfo){
    
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
    
    
    
}
// Collect the charge inside each object.
void oCollectObjectCharge(Population *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo) {
    
    double *val = rhoObj->val;
    long int *sizeProd = rhoObj->sizeProd;
    
    int nSpecies = pop->nSpecies;
    double *charge = pop->charge;
    
    long int *lookupIntOff = obj->lookupInteriorOffset;
    long int *lookupSurfOff = obj->lookupSurfaceOffset;
    
    // We might add this to the Object, although probably better  to store the rhoObj for restarts and insulators later on.
    double *chargeCounter = malloc(obj->nObjects*sizeof(*chargeCounter));
    adSetAll(chargeCounter,obj->nObjects,0);
    
    //double invNrSurfNod = 1.0/(obj->lookupSurfaceOffset[obj->nObjects]);
    double *invNrSurfNod = malloc(obj->nObjects*sizeof(*invNrSurfNod));
    for (long int a=0; a<obj->nObjects; a++) {
        invNrSurfNod[a] = 1.0/(obj->lookupSurfaceOffset[a+1]);
    }
    
    for(int s=0;s<nSpecies;s++) {
        
        long int iStart = pop->iStart[s];
        long int iStop = pop->iStop[s];
        
        for(int i=iStart;i<iStop;i++){
            
            double *pos = &pop->pos[3*i];
            double *vel = &pop->vel[3*i];
            
            // Integer parts of position
            int j = (int) pos[0];
            int k = (int) pos[1];
            int l = (int) pos[2];
            
            long int p = j + k*sizeProd[2] + l*sizeProd[3];
            
            // Check whether p is one of the object nodes and collect the charge if so.
            for (long int a=0; a<obj->nObjects; a++) {
                for (long int b=lookupIntOff[a]; b<lookupIntOff[a+1]; b++) {
                    if ((obj->lookupInterior[b])==p) {
                        chargeCounter[a] += charge[s];
                        pCut(pop, s, p, pos, vel);
                    }
                }
            }
        }
    }
    // Add the collected charge to the surface nodes on rhoObject.
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=lookupSurfOff[a]; b<lookupSurfOff[a+1]; b++) {
            val[obj->lookupSurface[b]] += chargeCounter[a]*invNrSurfNod[a];
        }
    }
}


/*****************************************************************************
 *  ALLOC/DESTRUCTORS
 ****************************************************************************/

Object *oAlloc(const dictionary *ini){
    
    Grid *domain = gAlloc(ini, SCALAR);
    
    long int *lookupInterior = NULL;
    long int *lookupInteriorOffset = NULL;
    long int *lookupSurface = NULL;
    long int *lookupSurfaceOffset = NULL;
    int nObjects = 0;
    
    Object *obj = malloc(sizeof(*obj));
    
    obj->domain = domain;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;
    obj->nObjects = nObjects;
    
    return obj;
}

void oFree(Object *obj){
    
    gFree(obj->domain);
    
    free(obj->lookupInterior);
    free(obj->lookupInteriorOffset);
    free(obj->lookupSurface);
    free(obj->lookupSurfaceOffset);
    free(obj);
    
}

void oCloseH5(Object *obj){
    
    gCloseH5(obj->domain);
}

void oOpenH5(const dictionary *ini, Object *obj, const MpiInfo *mpiInfo,
             const double *denorm, const double *dimen, const char *fName){
   
    gOpenH5(ini, obj->domain,   mpiInfo, denorm, dimen, "object");
}

void oReadH5(Object *obj, const MpiInfo *mpiInfo){
    
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
    
    //Communicate the boundary nodes -> DON'T DO THIS HERE!
    gHaloOp(setSlice, obj->domain, mpiInfo, TOHALO);
    
    //Count the number of objects and fill the lookup tables.
    oFillLookupTables(obj,mpiInfo);
    
    // Find all the object nodes which are part of the object surface.
    oFindObjectSurfaceNodes(obj, mpiInfo);
}

/******************************************************************************
 *  GLOBAL FUNCTION DEFINITIONS
 *****************************************************************************/


