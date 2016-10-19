/**
 * @file		object.h
 * @author		Jan Deca <...>
 * @brief		Implementing objects
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

Object *oAlloc(const dictionary *ini);
void oFree(Object *obj);
void oOpenH5(	const dictionary *ini, Object *obj, const MpiInfo *mpiInfo,
				const char *fName);
void oCloseH5(Grid *grid);

void oRayTrace(Population *pop, const Object *obj);

#endif // OBJECT_H
