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
 * @brief   Compute the capacitance matrix.
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
void oWoop(long int node, const int *nGhostLayersBefore,
           const int *nGhostLayersAfter, const int *trueSize,
           const long int *sizeProd, bool *woo);

/******************************************************************************
 *  LOCAL FUNCTION DEFINITIONS
 *****************************************************************************/


bool isGhostNode(Grid *grid, long int node) {
    
    //const double *val = grid->val;
    long int *sizeProd = grid->sizeProd;
    int *trueSize = grid->trueSize;
    int *nGhostLayers = grid->nGhostLayers;
    int rank = grid->rank;
    
    bool woo = false;
    
    oWoop(node,&nGhostLayers[rank-1],&nGhostLayers[2*rank-1],&trueSize[rank-1],&sizeProd[rank-1],&woo);

    return woo;
}

void oWoop(long int node, const int *nGhostLayersBefore,
           const int *nGhostLayersAfter, const int *trueSize,
           const long int *sizeProd, bool *woo) {
    
    if (*sizeProd==1) {
        if (node>*trueSize || node<1) {
            *woo = true;
        }
    } else {
        long int help = *(sizeProd) * (*(trueSize)+1);
        if (node < *(sizeProd) || node > help  ) {
            *woo = true;
        }
        node = node % *(sizeProd);
        oWoop(node,nGhostLayersBefore-1,
                    nGhostLayersAfter-1,trueSize-1,sizeProd-1,woo);
    }
}


void oFillLookupTables(Object *obj, const MpiInfo *mpiInfo) {
    
    int nObjects = 0;
    for (int i=0; i<obj->domain->sizeProd[obj->domain->rank]; i++) {
        if (obj->domain->val[i]>nObjects) {
            nObjects = (int)(obj->domain->val[i]+0.5); // Note, this is not necessary
                //the number of objects, but rather the identifier of the object with the highest number.
        }
    }
    // Make sure each process knows about the total number of objects.
    MPI_Allreduce(MPI_IN_PLACE, &nObjects, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    // Initialise and compute the array storing the offsets of the objects in the lookup table.
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
        if (obj->domain->val[i]>0.5){
            lookupInterior[(index[(int)(obj->domain->val[i]+0.5)-1])] = i;
            (index[(int)(obj->domain->val[i]+0.5)-1])++;
        }
    }
    
    // Add to the obj Object.
    obj->nObjects = nObjects;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;
}

// Put comments and documentation still
void oFindObjectSurfaceNodes(Object *obj, const MpiInfo *mpiInfo) {
    
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

void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini, const MpiInfo *mpiInfo) {
    
    long int *lookupSurface = obj->lookupSurface;
    long int *lookupSurfaceOffset = obj->lookupSurfaceOffset;
    
    // Allocate structures to run the potential solver
    MgAlgo mgAlgo = getMgAlgo(ini);
    Grid *rho = gAlloc(ini, SCALAR);
    Grid *res = gAlloc(ini, SCALAR);
    Grid *phi = gAlloc(ini, SCALAR);
    double *valrho = rho->val;
    double *valphi = phi->val;
    Multigrid *mgPhi = mgAlloc(ini, phi);
    Multigrid *mgRho = mgAlloc(ini, rho);
    Multigrid *mgRes = mgAlloc(ini, res);

    gZero(rho);
    
    int rank = mpiInfo->mpiRank;
    int size = mpiInfo->mpiSize;

    for (long int a=0; a<obj->nObjects; a++) {
        int total_global = 0;
        int total_local = lookupSurfaceOffset[a+1] - lookupSurfaceOffset[a];
        
        // MPI_LONG_INT sum not defined apparently....
        MPI_Allreduce(&total_local, &total_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        //printf("%ld, %ld\n",total_local, total_global);
        // total_global is the size of the capacitance matrix! -> so we should have done things globally after all...
        //total_global=5;
        double *capMatrix = malloc((total_global*total_global) * sizeof(*capMatrix));
        for (long int i=0; i<(total_global*total_global); i++) {
            capMatrix[i] = 0;
        }
        double *invCapMatrix = malloc((total_global*total_global) * sizeof(*invCapMatrix));
        for (long int i=0; i<(total_global*total_global); i++) {
            invCapMatrix[i] = 0;
        }
        int *totalPerNode = malloc((size+1)*sizeof(*totalPerNode));
        MPI_Allgather(&total_local, 1, MPI_INT, totalPerNode, 1, MPI_INT, MPI_COMM_WORLD);
        //aiPrint(totalPerNode,size+1);
        for(long int i=size;i>-1;i--) totalPerNode[i+1]=totalPerNode[i];
        totalPerNode[0] = 0;
        aiCumSum(totalPerNode+1,totalPerNode,size);
        aiPrint(totalPerNode,size+1);
        
        int j = 0;
        int inode = 0;
        
        for (int i=0; i<total_global; i++) {
            msg(STATUS,"Solving capacitance matrix for node %ld of %ld for object %ld of %ld.",i+1,total_global,a+1,obj->nObjects);
            //printf("%d %d %ld, a",j,inode,lookupSurface[inode]);
            if (rank==j) {
                valrho[lookupSurface[inode]] = 1;
            }
            //printf("b");
            mgSolve(mgAlgo, mgRho, mgPhi, mgRes, mpiInfo);
            //printf("c");
            if (rank==j) {
                valrho[lookupSurface[inode]] = 0;
            }
            //printf("d");
            
            //printf("e\n");
            //Now add the result to the matrix in the ith column -> this needs to be done globally!!!
         
            //if (rank==j) {
            //    printf("%f\n",phi->val[lookupSurface[inode]]);
               // for (long int k=0; k<total_global; k++) {
                  //  printf("%d %d, %f\n",inode,k, phi->val[lookupSurface[k]]);
                //    capMatrix[i+total_global*k] = phi->val[lookupSurface[k]];
            //}
            //}
            //printf("%d %d %d\n",rank, totalPerNode[rank],totalPerNode[rank+1]);
            for (int k=totalPerNode[rank]; k<totalPerNode[rank+1]; k++) {
                //printf("%d\n",k);
                // lookUpSurface is not the same on all nodes, that's why it hangs... NOT WORKING FOR MORE THAN ONE OBJECT!!!
                capMatrix[i+total_global*k] = phi->val[lookupSurface[k-totalPerNode[rank]]];
                //printf("%f %d\n",phi->val[lookupSurface[k]],k);
            }
 /*
            for (int k=totalPerNode[rank]; k<totalPerNode[rank+1]; k++) {
              //  //capMatrix[i+total_global*k] = phi->val[lookupSurface[k]];
                printf("%d %d %f\n",k,rank,phi->val[lookupSurface[k]]);
            }
*/
            inode++;
            if (inode>(totalPerNode[j+1]-totalPerNode[j]-1)) {
                j++;
                inode=0;
            }

        }
    
        MPI_Allreduce(MPI_IN_PLACE, capMatrix, (total_global*total_global), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //adPrint(capMatrix,(total_global*total_global));
        
        gsl_matrix_view A = gsl_matrix_view_array(capMatrix, total_global, total_global);
        gsl_matrix_view invA = gsl_matrix_view_array(invCapMatrix, total_global, total_global);
        
        int s;
        gsl_permutation *p = gsl_permutation_alloc (total_global);
        
        gsl_linalg_LU_decomp(&A.matrix, p, &s);
        
        gsl_linalg_LU_invert(&A.matrix, p, &invA.matrix);
    
    }

}



/*****************************************************************************
 *  ALLOC/DESTRUCTORS
 ****************************************************************************/

Object *oAlloc(const dictionary *ini){
    
    Grid *domain = gAlloc(ini, SCALAR);
    
    long int *lookupInterior = NULL;
    long int *lookupInteriorOffset = NULL;
    int nObjects = 0;
    
    Object *obj = malloc(sizeof(*obj));
    
    obj->domain = domain;
    obj->lookupInterior = lookupInterior;
    obj->lookupInteriorOffset = lookupInteriorOffset;
    obj->nObjects = nObjects;
    
    return obj;
}

void oFree(Object *obj){
    
    gFree(obj->domain);
    
    free(obj->lookupInterior);
    free(obj->lookupInteriorOffset);
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
    
    //Communicate the boundary nodes
    gHaloOp(setSlice, obj->domain, mpiInfo, TOHALO);
    
    // Count the number of objects and fills the lookup tables.
    oFillLookupTables(obj,mpiInfo);
    
    // Find all the object nodes which are part of the object surface.
    oFindObjectSurfaceNodes(obj, mpiInfo);
}

/******************************************************************************
 *  GLOBAL FUNCTION DEFINITIONS
 *****************************************************************************/


