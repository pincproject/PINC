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
	Grid *domain;					///< Represents presence of objects
	long int *lookupInterior;		///< Indices of the interior of the objects
	long int *lookupInteriorOffset;	///< Offset in the above per object (nObjects+1 elements)
  long int *lookupSurface;        ///< Indices of the surface nodes of the objects
  long int *lookupSurfaceOffset;  ///< Offset in the above per object (nObjects+1 elements)
  double *capMatrix;              ///< Array holding the capacitance matrices for each object
  double capMatrixInvSum;         ///< Array holding the total sum of capMatrix elements (nObjects elements)
  double *invNeedCoffeeMatrix;    ///< Better name wanted...
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
 * @brief   Compute the capacitance matrix. (one big matrix containing all objects)
 * @param	obj		Object
 * @param	ini		input settings
 * @return	void
 *
 * Compute the capacitance matrix.
 */
void oComputeCapacitanceMatrix(Object *obj, const dictionary *ini,
                               const MpiInfo *mpiInfo);

/**
 * @brief	Apply the capacitance matrix
 * @param   rho         Grid
 * @param   phi         Grid
 * @param	obj         Object
 * @param	mpiInfo		MpiInfo
 * @return	void
 *
 * Construct and solve equation 5 in Miyake_Usui_PoP_2009.
 */
void oApplyCapacitanceMatrix(Grid *rho, const Grid *phi, const Object *obj,
                             const MpiInfo *mpiInfo);

/**
 * @brief	Collect the charge inside each object
 * @param   pop         Population
 * @param   rhoObj      Grid
 * @param	obj         Object
 * @param	mpiInfo		MpiInfo
 * @return	void
 *
 * Collect the charge inside each object.
 */
void oCollectObjectCharge(Population *pop, Grid *rhoObj, Object *obj,
                          const MpiInfo *mpiInfo);

/**
 * TO IMPLEMENT!
 *
 */
void oRayTrace(Population *pop, const Object *obj);

#endif // OBJECT_H







