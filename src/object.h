/**
 * @file		object.h
 * @author		Jan Deca <jandeca@gmail.com>
 * @brief		All object-related functions are here.
 * @date		19.10.16
 */

#ifndef OBJECT_H
#define OBJECT_H

/**
 * @brief Represents an object
 */
typedef struct{
	Grid *domain;					///< Represents precense of objects
	long int *lookupInterior;		///< Indices of the interior of the objects
	long int *lookupInteriorOffset;	///< Offset in the above per object (nObjects+1 elements)
    long int *lookupSurface;        ///< Indices of the surface nodes of the objects
    long int *lookupSurfaceOffset;  ///< Offset in the above per object (nObjects+1 elements)
    double *capMatrixAll;           ///< Array holding the capacitance matrices for each object
    long int *capMatrixAllOffsets;  ///< Offset in the above per object (nObjects*(size+1) elements)
    double *capMatrixSum;           ///< Array holding the total sum of capMatrix elements (nObjects elements)
	int nObjects;					///< Number of objects
} Object;

/**
 * @brief Allocates an Object object
 * @param	ini			Input file
 * @return				Pointer to Object
 *
 * Remember to free using oFree().
 *
 * NB! Assumes 1 ghost point on all edges for now.
 */
Object *oAlloc(const dictionary *ini);

/**
 * @brief Frees allocated object
 * @param	obj     Object
 * @return	void
 */
void oFree(Object *obj);

/**
 * @brief	Opens .grid.h5-file to read in object
 * @param	ini				Dictionary to input file
 * @param	obj             Object
 * @param	mpiInfo			MpiInfo
 * @param	denorm			Quantity denormalization factors
 * @param	dimen			Quantity dimensionalizing factors
 * @param	fName			Filename
 * @return	void
 * @see gOpenH5()
 *
 * Remember to call oCloseH5().
 *
 */
void oOpenH5(const dictionary *ini, Object *obj, const MpiInfo *mpiInfo,
             const double *denorm, const double *dimen, const char *fName);

/**
 * @brief	Closes .grid.h5-file
 * @param	obj     Object
 * @return	void
 */
void oCloseH5(Object *obj);

/**
 * @brief	Read values from .grid.h5-field to Object
 * @param	obj             Object
 * @param	mpiInfo			MpiInfo
 * @return	void
 * @see gReadH5()
 *
 * Reads the input objects and creates the various lookup tables needed.
 */
void oReadH5(Object *obj, const MpiInfo *mpiInfo);

/**
 * TO IMPLEMENT!
 *
 */
void oRayTrace(Population *pop, const Object *obj);

void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini, const MpiInfo *mpiInfo);

void oApplyCapacitanceMatrix(Grid *rho, const Grid *phi, const Object *obj, const MpiInfo *mpiInfo);

#endif // OBJECT_H







