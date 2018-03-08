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
 * @param	axisDenorm		Axis denormalization factor
 * @param	quantityDenorms	Quantity denormalization factor
 * @param	fName			Filename
 * @return	void
 * @see gOpenH5()
 *
 * Remember to call oCloseH5().
 *
 */
void oOpenH5(const dictionary *ini, Object *obj, const MpiInfo *mpiInfo,
             const Units *units, double denorm, const char *fName);
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

#endif // OBJECT_H
