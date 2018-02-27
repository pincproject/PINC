/**
 * @file		object.c
 * @author		Jan Deca <jandeca@gmail.com>
 * @brief		All object-related functions are here.
 * @date		19.10.16
 */

#include "core.h"
#include "object.h"

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

/******************************************************************************
 *  LOCAL FUNCTION DEFINITIONS
 *****************************************************************************/

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
             const Scales *scales, double denorm, const char *fName){
   
    gOpenH5(ini, obj->domain,   mpiInfo, scales, denorm, "object");
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
    
    //Count the number of objects and fills the lookup tables.
    oFillLookupTables(obj,mpiInfo);
}

/******************************************************************************
 *  GLOBAL FUNCTION DEFINITIONS
 *****************************************************************************/


