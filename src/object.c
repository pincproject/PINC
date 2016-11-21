/**
 * @file		object.c
 * @author		Jan Deca <jandeca@gmail.com>
 * @brief		All object-related functions are here.
 * @date		19.10.16
 */

#include "core.h"
#include "object.h"
#include "multigrid.h"
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
 * @brief   Find all the object nodes which are part of the object surface.
 * @param	obj		Object
 * @param	ini		input settings
 * @return	void
 */
void oFindObjectSurfaceNodes(Object *obj, const MpiInfo *mpiInfo);

/**
 * @brief   Compute the capacitance matrix for each object.
 * @param	obj		Object
 * @param	ini		input settings
 * @return	void
 */
void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini, const MpiInfo *mpiInfo);


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

//Count the number of objects and fills the lookup tables.
void oFillLookupTables(Object *obj, const MpiInfo *mpiInfo) {
    
    // Find the number of objects in the input file.
    int nObjects = 0;
    for (int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i] > nObjects) {
            nObjects = (int)(obj->domain->val[i]+0.5); // Note, this is not necessarily
                // the number of objects, but rather the identifier of the object with the highest number.
        }
    }
    // Make sure each process knows about the total number of objects.
    MPI_Allreduce(MPI_IN_PLACE, &nObjects, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    // Initialise and compute the array storing the offsets of the objects in the lookup table.
    long int *lookupInteriorOffset = malloc((nObjects+1)*sizeof(*lookupInteriorOffset));
    alSetAll(lookupInteriorOffset, nObjects+1, 0);
   
    for (long int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>0.5 && !isGhostNode(obj->domain, i)){
            lookupInteriorOffset[(int)(obj->domain->val[i]+0.5)]++;
        }
    }
    alCumSum(lookupInteriorOffset+1,lookupInteriorOffset,nObjects);
    
    // Initialise and compute the lookup table.
    long int *lookupInterior = malloc((lookupInteriorOffset[nObjects])*sizeof(*lookupInterior));
    alSetAll(lookupInterior,lookupInteriorOffset[nObjects],0);
   
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
    
    // Add to the object.
    obj->nObjects = nObjects;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;

}

//Find all the object nodes which are part of the object surface.
void oFindObjectSurfaceNodes(Object *obj, const MpiInfo *mpiInfo) {
    
    long int *sizeProd = obj->domain->sizeProd;
    
    // Initialise the array storing the offsets for the surface nodes in the lookup table.
    long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*sizeof(*lookupSurfaceOffset));
    alSetAll(lookupSurfaceOffset,obj->nObjects+1,0);
    
    // Find the 8 neighbour cells of each non-ghost node.
    long int *myNeighbours = malloc(9*sizeof(*myNeighbours));
    // Find the ofsetts first.
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<sizeProd[obj->domain->rank]; b++) {
            if (!isGhostNode(obj->domain, b)) {
                myNeighbours[0] = b;    // me on node i,j,k
                myNeighbours[1] = myNeighbours[0];                  // cell i,j,k
                myNeighbours[2] = myNeighbours[0] - sizeProd[3];    // cell i,j,k-1
                myNeighbours[3] = myNeighbours[0] - sizeProd[1];                    // cell i-1,j,k
                myNeighbours[4] = myNeighbours[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNeighbours[5] = myNeighbours[0] - sizeProd[2];                    // cell i,j-1,k
                myNeighbours[6] = myNeighbours[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNeighbours[7] = myNeighbours[0] - sizeProd[2] - sizeProd[1];                  // cell i-1,j-1,k;
                myNeighbours[8] = myNeighbours[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];    // cell i-1,j-1,k-1
                
                int d=0;
                if (obj->domain->val[myNeighbours[1]]>(a+0.5) && obj->domain->val[myNeighbours[1]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[2]]>(a+0.5) && obj->domain->val[myNeighbours[2]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[3]]>(a+0.5) && obj->domain->val[myNeighbours[3]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[4]]>(a+0.5) && obj->domain->val[myNeighbours[4]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[5]]>(a+0.5) && obj->domain->val[myNeighbours[5]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[6]]>(a+0.5) && obj->domain->val[myNeighbours[6]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[7]]>(a+0.5) && obj->domain->val[myNeighbours[7]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[8]]>(a+0.5) && obj->domain->val[myNeighbours[8]]<(a+1.5)) d++;

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
                myNeighbours[0] = b;    // me on node i,j,k
                myNeighbours[1] = myNeighbours[0];                  // cell i,j,k
                myNeighbours[2] = myNeighbours[0] - sizeProd[3];    // cell i,j,k-1
                myNeighbours[3] = myNeighbours[0] - sizeProd[1];                    // cell i-1,j,k
                myNeighbours[4] = myNeighbours[0] - sizeProd[1] - sizeProd[3];      // cell i-1,j,k-1
                myNeighbours[5] = myNeighbours[0] - sizeProd[2];                    // cell i,j-1,k
                myNeighbours[6] = myNeighbours[0] - sizeProd[2] - sizeProd[3];      // cell i,j-1,k-1
                myNeighbours[7] = myNeighbours[0] - sizeProd[2] - sizeProd[1];                  // cell i-1,j-1,k;
                myNeighbours[8] = myNeighbours[0] - sizeProd[2] - sizeProd[1] - sizeProd[3];    // cell i-1,j-1,k-1
                
                int d=0;
                if (obj->domain->val[myNeighbours[1]]>(a+0.5) && obj->domain->val[myNeighbours[1]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[2]]>(a+0.5) && obj->domain->val[myNeighbours[2]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[3]]>(a+0.5) && obj->domain->val[myNeighbours[3]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[4]]>(a+0.5) && obj->domain->val[myNeighbours[4]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[5]]>(a+0.5) && obj->domain->val[myNeighbours[5]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[6]]>(a+0.5) && obj->domain->val[myNeighbours[6]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[7]]>(a+0.5) && obj->domain->val[myNeighbours[7]]<(a+1.5)) d++;
                if (obj->domain->val[myNeighbours[8]]>(a+0.5) && obj->domain->val[myNeighbours[8]]<(a+1.5)) d++;
                
                // Check if on surface.
                if (d<7.5 && d>0) {
                    lookupSurface[index[a]] = myNeighbours[0];
                    index[a]++;
                }
            }
        }
    }
    
    // Add to object.
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;
}

// Compute the capacitance matrix for each object.
void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini, const MpiInfo *mpiInfo) {
    
    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurface = obj->lookupSurface;
    long int *lookupSurfaceOffset = obj->lookupSurfaceOffset;
    
    // Allocate and initialise the structures to run the potential solver.
    MgAlgo mgAlgo = getMgAlgo(ini);
    Grid *rho = gAlloc(ini, SCALAR);
    Grid *res = gAlloc(ini, SCALAR);
    Grid *phi = gAlloc(ini, SCALAR);
    Multigrid *mgPhi = mgAlloc(ini, phi);
    Multigrid *mgRho = mgAlloc(ini, rho);
    Multigrid *mgRes = mgAlloc(ini, res);

    gZero(rho);
    
    // Find the number of surface nodes for each object.
    long int *nodesCoreLocal = malloc((size+1)*sizeof(*nodesCoreLocal));
    long int *nodesCoreGlobal = malloc(obj->nObjects*(size+1)*sizeof(*nodesCoreGlobal));
    
    double *capMatrixSum = malloc(obj->nObjects*sizeof(*capMatrixSum));

    for (long int a=0; a<obj->nObjects; a++) {
        
        long int nodesThisCore = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];

        // Let every core know how many surface nodes everybody has.
        MPI_Allgather(&nodesThisCore, 1, MPI_LONG, nodesCoreLocal, 1, MPI_LONG, MPI_COMM_WORLD);
        
        for(long int i=size-1;i>-1;i--) nodesCoreLocal[i+1]=nodesCoreLocal[i];
        nodesCoreLocal[0] = 0;
        alCumSum(nodesCoreLocal+1,nodesCoreLocal,size);

        for (long int b=0; b<size+1; b++) nodesCoreGlobal[a*(size+1)+b] = nodesCoreLocal[b];
    }
    
    // Find the size and initialise the array holding the capacitance matrices for all objects.
    long int capMatrixAllSize = 0;
    for (long int a=0; a<obj->nObjects; a++) {
        capMatrixAllSize +=nodesCoreGlobal[a*(size+1)+size];
    }
    double *capMatrixAll = malloc(capMatrixAllSize*capMatrixAllSize*sizeof(*capMatrixAll));
    //long int *capMatrixAllOffsets = malloc(obj->nObjects*(size+1)*sizeof(*capMatrixAllOffsets));

    // Compute the actual capacitance matrix for each object.
    for (long int a=0; a<obj->nObjects; a++) {
        long int j = 0; // Keep track of the rank
        long int inode = 0; // Keep track of the node number

        // Initialise the matrix and its inverse.
        double *capMatrix = malloc((nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]) * sizeof(*capMatrix));
        double *invCapMatrix = malloc((nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]) * sizeof(*invCapMatrix));
        adSetAll(capMatrix,nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size],0);
        adSetAll(invCapMatrix,nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size],0);

        // Loop over the nodes and fill the matrix
        for (long int i=0; i<nodesCoreGlobal[a*(size+1)+size]; i++) {
            msg(STATUS,"Solving capacitance matrix for node %ld of %ld for object %ld of %ld.", \
                i+1,nodesCoreGlobal[a*(size+1)+size],a+1,obj->nObjects);
            // Set the surface node to 1 charge.
            if (rank==j) {
                rho->val[lookupSurface[inode]] = 1;
            }
            // Solve for the potential.
            mgSolve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
            // Set the surface node back to zero.
            if (rank==j) {
                rho->val[lookupSurface[inode]] = 0;
            }
            // Fill column i of the capacitance matrix.
            for (int k=nodesCoreGlobal[a*(size+1)+rank]; k<nodesCoreGlobal[a*(size+1)+rank+1]; k++) {
                capMatrix[i+nodesCoreGlobal[a*(size+1)+size]*k] = \
                    phi->val[lookupSurface[k-nodesCoreGlobal[a*(size+1)+rank]]];
            }
            
            // Increase the counters. If you looped over all nodes on this core, increase the rank and reset inode.
            inode++;
            if (inode>(nodesCoreGlobal[a*(size+1)+j+1]-nodesCoreGlobal[a*(size+1)+j]-1)) {
                j++;
                inode=0;
            }

        }
    
        // Make sure every codes has the complete matrix (needed for BLAS).
        MPI_Allreduce(MPI_IN_PLACE, capMatrix, (nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]), \
                      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 
        // Compute the inverse of the capacitance matrix.
        // Actually, the inverse is the capacitance matrix. Probably have to rethink the variable names.
        gsl_matrix_view A = gsl_matrix_view_array(capMatrix, nodesCoreGlobal[a*(size+1)+size], nodesCoreGlobal[a*(size+1)+size]);
        gsl_matrix_view invA = gsl_matrix_view_array(invCapMatrix, nodesCoreGlobal[a*(size+1)+size], nodesCoreGlobal[a*(size+1)+size]);
        
        int s;
        gsl_permutation *p = gsl_permutation_alloc (nodesCoreGlobal[a*(size+1)+size]);
        gsl_linalg_LU_decomp(&A.matrix, p, &s);
        gsl_linalg_LU_invert(&A.matrix, p, &invA.matrix);
        
        // Add the invCapMatrix for object a to the big array.
        for (long int l=0; l<nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]; l++) {
            capMatrixAll[a*nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]+l] = \
                invCapMatrix[l];
        }
        
        // Compute here to total sum of elements in the capacitance matrix (needed later).
        capMatrixSum[a] = adSum(invCapMatrix,nodesCoreGlobal[a*(size+1)+size]*nodesCoreGlobal[a*(size+1)+size]);
        // We need the inverse later on.
        capMatrixSum[a] = 1/capMatrixSum[a];

    }

    long int *capMatrixAllOffsets = nodesCoreGlobal;
    
    // Add to object
    obj->capMatrixAll = capMatrixAll;
    obj->capMatrixAllOffsets = capMatrixAllOffsets;
    obj->capMatrixSum = capMatrixSum;
}

// Construct and solve equation 5 in Miyake_Usui_PoP_2009
void oApplyCapacitanceMatrix(Grid *rho, const Grid *phi, const Object *obj, const MpiInfo *mpiInfo){

    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;
    long int *lookupSurface = obj->lookupSurface;
    
    double *capMatrixAll = obj->capMatrixAll;
    long int *capMatrixAllOffsets = obj->capMatrixAllOffsets;
    double *capMatrixSum = obj->capMatrixSum;
    
    double capMatrixPhiSum = 0;
    
    // Loop over the objects
    for (long int a=0; a<obj->nObjects; a++) {
        
        // total number of surface nodes
        long int totSurfNodGlobal = capMatrixAllOffsets[a*(size+1)+size];
        long int beginIndex = capMatrixAllOffsets[a*(size+1)+rank];
        long int endIndex = capMatrixAllOffsets[a*(size+1)+rank+1];

        
        double *deltaPhi = malloc(totSurfNodGlobal*sizeof(*deltaPhi));
        adSetAll(deltaPhi,totSurfNodGlobal,0);
        double *rhoCorrection = malloc(totSurfNodGlobal*sizeof(*rhoCorrection));
        adSetAll(rhoCorrection,totSurfNodGlobal,0);
        
        // Compute eq. 7.
        for (long int i=0; i<totSurfNodGlobal; i++) {
            // Make sure that each core loops only over the matrix elements/parts of the grid it has
            for (long int j=beginIndex; j<endIndex; j++) {
                capMatrixPhiSum += capMatrixAll[i+totSurfNodGlobal*j] * (phi->val[lookupSurface[j-beginIndex]]);
            }
        }
        // This is phi_c for each object.
        capMatrixPhiSum = capMatrixPhiSum*capMatrixSum[a];
        MPI_Allreduce(MPI_IN_PLACE, &capMatrixPhiSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        msg(STATUS,"Potential-check for object %ld : %f",a,capMatrixPhiSum);
        //capMatrixPhiSum=0.03;

        for (long int j=beginIndex; j<endIndex; j++) {
            deltaPhi[j] = capMatrixPhiSum - phi->val[lookupSurface[j-beginIndex]];
        }
        MPI_Allreduce(MPI_IN_PLACE, deltaPhi, totSurfNodGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        // Eq. 5
        for (long int i=0; i<totSurfNodGlobal; i++) {
            for (long int j=beginIndex; j<endIndex; j++) {
                rhoCorrection[i] += capMatrixAll[i+totSurfNodGlobal*j]*deltaPhi[j];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, rhoCorrection, totSurfNodGlobal, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        for (long int j=beginIndex; j<endIndex; j++) {
            
            if (rhoCorrection[j]==0) {
                printf("trouble");
            }
            rho->val[lookupSurface[j-beginIndex]] += rhoCorrection[j];
        }
    }
}

void oCollectObjectCharge(Population *pop, Grid *rhoObj, Object *obj, const MpiInfo *mpiInfo) {

    double *val = rhoObj->val;
    long int *sizeProd = rhoObj->sizeProd;
    
    int nSpecies = pop->nSpecies;
    double *charge = pop->charge;
    
    // We might add this to the Object, although probably better  to store the rhoObj for restarts and insulators later on.
    double *chargeCounter = malloc(obj->nObjects*sizeof(*chargeCounter));
    adSetAll(chargeCounter,obj->nObjects,0);

    //double invNrSurfaceNodes = 1.0/(obj->lookupSurfaceOffset[obj->nObjects]);
    double *invNrSurfaceNodes = malloc(obj->nObjects*sizeof(*invNrSurfaceNodes));
    for (long int a=0; a<obj->nObjects; a++) {
        invNrSurfaceNodes[a] = 1.0/(obj->lookupSurfaceOffset[a+1]);
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

            // Check wether p is one of the object nodes and collect the charge if so.
            for (long int a=0; a<obj->nObjects; a++) {
                for (long int b=obj->lookupInteriorOffset[a]; b<obj->lookupInteriorOffset[a+1]; b++) {
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
        for (long int b=obj->lookupSurfaceOffset[a]; b<obj->lookupSurfaceOffset[a+1]; b++) {
            val[obj->lookupSurface[b]] += chargeCounter[a]*invNrSurfaceNodes[a];
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
    
    double *capMatrixAll = NULL;
    long int *capMatrixAllOffsets = NULL;
    double *capMatrixSum = NULL;
    
    Object *obj = malloc(sizeof(*obj));
    
    obj->domain = domain;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;
    obj->nObjects = nObjects;
    
    obj->capMatrixAll = capMatrixAll;
    obj->capMatrixAllOffsets =capMatrixAllOffsets;
    obj->capMatrixSum = capMatrixSum;
    
    return obj;
}

void oFree(Object *obj){
    
    gFree(obj->domain);
    
    free(obj->lookupInterior);
    free(obj->lookupInteriorOffset);
    free(obj->lookupSurface);
    free(obj->lookupSurfaceOffset);
    free(obj->capMatrixAll);
    free(obj->capMatrixAllOffsets);
    free(obj->capMatrixSum);
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
    
    // Count the number of objects and fills the lookup tables.
    oFillLookupTables(obj,mpiInfo);
    
    // Find all the object nodes which are part of the object surface.
    oFindObjectSurfaceNodes(obj, mpiInfo);
}
















/******************************************************************************
 *  DEPRECIATED FUNCTION DEFINITIONS
 *****************************************************************************/

//depreciated...
void oFindObjectSurfaceNodes_v1(Object *obj, const MpiInfo *mpiInfo) {
    
    long int *sizeProd = obj->domain->sizeProd;
    
    //if(isGhostNode(obj->domain, 10439)) printf("test\n");
    
    long int *lookupSurfaceOffset = malloc((obj->nObjects+1)*sizeof(*lookupSurfaceOffset));
    for (long int i=0; i<obj->nObjects+1; i++) {
        lookupSurfaceOffset[i] = 0;
    }
    
    long int *meNeighbours = malloc(7*sizeof(*meNeighbours));
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
            if (obj->domain->val[b]>0.5 && !isGhostNode(obj->domain, b)) {
                // My neighbours
                meNeighbours[0] = b; // me
                meNeighbours[1] = meNeighbours[0] + sizeProd[1]; // right
                meNeighbours[2] = meNeighbours[0] - sizeProd[1]; // left
                meNeighbours[3] = meNeighbours[0] + sizeProd[2]; // up
                meNeighbours[4] = meNeighbours[0] - sizeProd[2]; // down
                meNeighbours[5] = meNeighbours[0] + sizeProd[3]; // back
                meNeighbours[6] = meNeighbours[0] - sizeProd[3]; // front
                
                int d=0;
                if (obj->domain->val[meNeighbours[1]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[2]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[3]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[4]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[5]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[6]]>(a+0.5)) d++;
                
                // Check if surface
                if (d<5.5) {
                    lookupSurfaceOffset[a+1]++;
                }
            }
        }
        //MPI_Allreduce(MPI_IN_PLACE, &lookupSurfaceOffset[a+1], 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    
    alCumSum(lookupSurfaceOffset+1,lookupSurfaceOffset,obj->nObjects);
    
    
    //Second go through to fill the table.
    long int *lookupSurface = malloc((lookupSurfaceOffset[obj->nObjects])*sizeof(*lookupSurface));
    for (long int i=0; i<lookupSurfaceOffset[obj->nObjects]+1; i++) {
        lookupSurface[i]=0;
    }
    
    long int *index = malloc((obj->nObjects)*sizeof(*index));
    for (long int i=0; i<obj->nObjects; i++) {
        index[i]=lookupSurfaceOffset[i];
    }
    
    for (long int a=0; a<obj->nObjects; a++) {
        for (long int b=0; b<obj->domain->sizeProd[obj->domain->rank]; b++) {
            if (obj->domain->val[b]>0.5 && !isGhostNode(obj->domain, b)) {
                // My neighbours
                meNeighbours[0] = b; // me
                meNeighbours[1] = meNeighbours[0] + sizeProd[1]; // right
                meNeighbours[2] = meNeighbours[0] - sizeProd[1]; // left
                meNeighbours[3] = meNeighbours[0] + sizeProd[2]; // up
                meNeighbours[4] = meNeighbours[0] - sizeProd[2]; // down
                meNeighbours[5] = meNeighbours[0] + sizeProd[3]; // back
                meNeighbours[6] = meNeighbours[0] - sizeProd[3]; // front
                
                int d=0;
                if (obj->domain->val[meNeighbours[1]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[2]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[3]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[4]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[5]]>(a+0.5)) d++;
                if (obj->domain->val[meNeighbours[6]]>(a+0.5)) d++;
                
                // Check if boundary
                if (d<5.5) {
                    lookupSurface[index[a]] = meNeighbours[0];
                    index[a]++;
                }
            }
        }
    }
    
    // Add to object
    obj->lookupSurface = lookupSurface;
    obj->lookupSurfaceOffset = lookupSurfaceOffset;
}

