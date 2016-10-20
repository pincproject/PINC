/**
 * @file		object.c
 * @author		Jan Deca <...>
 * @brief		Implementing objects
 * @date		19.10.16
 */

#include "core.h"
#include "object.h"

/******************************************************************************
 * LOCAL FUNCTION DECLARATIONS
 *****************************************************************************/

/*****************************************************************************
 *		ALLOC/DESTRUCTORS
 ****************************************************************************/

Object *oAlloc(const dictionary *ini){
    
    //Call gAlloc
    Grid *domain = gAlloc(ini, SCALAR);
    
    long int *lookupInterior = NULL;
    long int *lookupInteriorOffset = NULL;
    int nObjects = 0;
    
    /* Store in Object */
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
    
    hid_t fileSpace = obj->domain->h5FileSpace;
    hid_t memSpace = obj->domain->h5MemSpace;
    hid_t file = obj->domain->h5;
    double *val = obj->domain->val;
    
    // Enable collective datawriting
    hid_t pList = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(pList, H5FD_MPIO_COLLECTIVE);
    
    char name[64];
    sprintf(name,"Object");
    
    hid_t dataset = H5Dopen(file,name,H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memSpace, fileSpace, pList, val);
    
    H5Dclose(dataset);
    H5Pclose(pList);
}

/******************************************************************************
 * FUNCTION DEFINITIONS
 *****************************************************************************/


