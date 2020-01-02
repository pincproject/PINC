/**
 * @file		grid.h
 * @brief		Input/output functions.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 * 				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 */

#ifndef GRID_H
#define GRID_H

/**
 * @brief Defines valid values to use with gAlloc()
 */
enum{
	SCALAR = 1,
	VECTOR = -1
};

/**
 * @brief Defines direction of operation in gHaloOp(), gHaloOpDim()
 */
typedef enum{
	TOHALO = 0,
	FROMHALO = 1
} opDirection;

/**
 * @brief Allocates a Grid object as specified in the input file
 * @param	ini			Input file
 * @param	nValues		Number of values per grid point (use SCALAR or VECTOR)
 * @return				Pointer to Grid
 *
 * nValues=1 for scalars or the number of dimensions for vectors. For
 * convenience, use SCALAR and VECTOR rather than specifying the numbers
 * manually.
 *
 * Remember to free using gFree().
 *
 * NB! Assumes 1 ghost point on all edges for now.
 */

Grid *gAlloc(const dictionary *ini, int nValues, const MpiInfo *mpiInfo);

/**
 * @brief Frees allocated grid
 * @param	grid	Grid
 * @return	void
 */
void gFree(Grid *grid);

/**
 * @brief Set boundary slices
 * @param   grid    Grid
 * @param   mpiInfo Mpi Info
 * @return  void
 *
 * Set the boundary slices. With Dirichlet and Neumann bnd conditions the
 * values which corresponds to phi and the gradient of phi is set here.
 *
 * For now it only has the capability of constant boundaries, but it could
 * later be expanded without too much trouble.
 *
 *
 */

void gSetBndSlices(const dictionary *ini, Grid *grid,const MpiInfo *mpiInfo);

/**
 * @brief Allocates the memory for an MpiInfo struct according to input file
 * @param	ini		Input file dictionary
 * @return	Pointer to MpiInfo
 */
MpiInfo *gAllocMpi(const dictionary *ini);

/**
 * @brief Frees the memory of an MpiInfo struct
 * @param	mpiInfo		MpiInfo
 * @return	void
 */
void gFreeMpi(MpiInfo *mpiInfo);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param sliceOp			Slicing operation
 * @param *Grid	Grid struct
 * @param *mpiInfo		MpiInfo struct
 * @param d				Along which dimension it should exhange ghost cells
 *
 * Since each true subdomain is surrounded by ghost layers, both to implement
 * boundary conditions and to facility communication between different
 * processes each taking care of a spesified subdomain, this function fills the
 * ghost layers with the values from the surrounding subdomains. When called it
 * has the option to either set the ghost layer as neighboring subdomain, or
 * add to the existing ghost layers. This is specified by the first argument
 * which can either be setSlice or addSlice.
 *
 *
 * Example with 2 subdomains and a 2D 2x2 true grid:
 * @code
  	1	1	1	1			2	2	2	2
 	1	1	1	1			2	2	2	2
 	1	1	1	1			2	2	2	2
 	1	1	1	1			2	2	2	2
 * @endcode
 *
 * If we then want to the ghost layers in the subdomains with the values from
 * the surrounding subdomain we use the following syntax.
 *
 * @code
	gHaloOpDim(setSlice, grid, mpiInfo, 1)
 * @endcode
 * The result should look like:
 *
 * @code
  	2	1	1	2			1	2	2	1
 	2	1	1	2			1	2	2	1
 	2	1	1	2			1	2	2	1
 	2	1	1	2			1	2	2	1
 * @endcode
 *
 * If we go back to the original grids and we instead want to add to the
 * ghostlayer we can use addSlice instead.
 * @code
	gHaloOpDim(addSlice, grid, mpiInfo, 1)
 * @endcode
 * This should produce:
 * @code
  	3	1	1	3			3	2	2	3
 	3	1	1	3			3	2	2	3
 	3	1	1	3			3	2	2	3
 	3	1	1	3			3	2	2	3
 * @endcode
 *
 * If needed it should be quick to facilitate for more slice operations, in
 * addition to set and add.
 *
 * NB! Only works with 1 ghost layer.
 * @see gHaloOp
 */
void gHaloOpDim(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, opDirection dir);


/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param sliceOp			Slicing operation
 * @param *grid				Grid struct
 * @param *mpiInfo			MpiInfo struct
 *
 * A wrapper to the gHaloOpDim function, that is used when the user wants the
 * interaction in all the dimensions.
 *
 * NB! Only works with 1 ghost layer.
 * @see gExchangeSlice
 * @see gHaloOpDim
 */
void gHaloOp(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, opDirection dir);



/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param sliceOp			Slicing operation
 * @param *Grid	Grid struct
 * @param *mpiInfo		MpiInfo struct
 * @param d				Along which dimension it should exhange ghost cells
 * @param int boundary  Decidec upper = 1 or lower = 0
 *
 * NB! Only works with 1 ghost layer.
 * @see gHaloOp
 */
void gHaloOpDimOpen(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, opDirection dir,int boundary);


/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param sliceOp			Slicing operation
 * @param *grid				Grid struct
 * @param *mpiInfo			MpiInfo struct
 *
 * A wrapper to the gHaloOpDim function, that is used when the user wants the
 * interaction in all the dimensions.
 *
 * NB! Only works with 1 ghost layer.
 * @see gExchangeSlice
 * @see gHaloOpDim
 */
void gHaloOpOpen(funPtr sliceOp, Grid *grid, const MpiInfo *mpiInfo, opDirection dir);




/**
 * @brief Extracts a (dim-1) dimensional slice of grid values.
 * @param	slice 		Return array
 * @param	grid		Grid
 * @param	d			Perpendicular direction to slice
 * @param	offset 		Offset of slice
 * @return				Void
 *
 * This function gets extracts a slice from a N dimensional grid. The integer d
 * decides in which direction the slice is perpendicular to, and the offset
 * decides which which slice it picks out. It needs a preallocated slice array
 * where the extracted slice will be stored.
 *
 * 2D example: Here we have a 5x4 grid and we want to extract a slice
 * corresponding to the second row, where x (d=0) is a constant 1.
 *
 * @code
 * 15   16   17   18   19
 *
 * 10   11   12   13   14
 *
 *  5    6    7    8    9
 *
 *  0    1    2    3    4
 * @endcode
 *
 * @code
	 getSlice(slice, grid, 0, 1);
 * @endcode
 * After running this the slice array consists of
 * slice = \f( [1, 6, 11, 16] \f)
 *
 * @see setSlice
 * @see gHaloOpDim
 **/

void getSlice(double *slice, const Grid *grid, int d, int offset);

/**
 * @brief places a (dim-1) dimensional slice onto a selected slice on the grid.
 * @param	slice		Slice containing a layer of values
 * @param	grid		Grid
 * @param	d 			Perpendicular direction to slice
 * @param	offset 		Offset of slice
 * @return				Void
 *
 * This function places a a slice on a grid. If we have a slice and want to
 * insert it onto a grid this function is used.
 *
 * Example: We have a 1D slice consisting of 6 2s and want to insert it onto
 * the third row, from the bottom.
 * @code
 *	111111
 *	111111
 *	111111
 *	111111
 *
 *	setSlice(slice, grid, 0, 2);
 *
 *	111111
 *	222222
 *	111111
 *	111111
 * @endcode
 *
 * @see setSlice
 * @see gHaloOpDim
 */
void setSlice(const double *slice, Grid *grid, int d, int offset);

/**
 * @brief Adds a slice to a slice in a Grid
 * @param	slice		Slice of values to add into grid
 * @param	grid		Grid
 * @param	d			Perpendicular direction to slice grid
 * @param	offset		Offset of slice in grid
 * @return				void
 *
 * Similar to setSlice() but adds slice to existing values rather than replacing
 * them.
 */
void addSlice(const double *slice, Grid *grid, int d, int offset);


/**
 * @brief Set all values in grid to zero
 * @param	grid	Grid
 * @return			void
 */
void gZero(Grid *grid);

/**
 * @brief Set grid quantity to vector (or scalar) for all grid points
 * @param	grid	Grid
 * @param	value	Array (vector) of values to set
 *
 * Each grid point is set to have the vector value specified by 'value'. Hence
 * value is expected to have length grid->size[0]
 */
void gSet(Grid *grid, const double *value);

/**
 * @brief Copy a grid
 * @param	original	 Original grid
 * @param	copy	     Copied grid
 *
 * Each grid point is set to have the vector value specified by 'value'. Hence
 * value is expected to have length grid->size[0]
 */
void gCopy(const Grid *original, Grid *copy);

/**
 * @brief Multiply all values in grid by a number
 * @param	grid	Grid
 * @param	num		Number to multiply by
 * @return			void
 */
void gMul(Grid *grid, double num);

/**
 * @brief Add all values in grid by a number
 * @param	grid	Grid
 * @param	num		Number to multiply by
 * @return			void
 */
void gAdd(Grid *grid, double num);

/**
 * @brief Add all values in grid by a number
 * @param	grid	Grid
 * @param	num		Number to multiply by
 * @return			void
 */
void gSub(Grid *grid, double num);

/**
 * @brief   Square a grid
 * @param	grid	Grid
 * @return			void
 */

void gSquare(Grid *grid);

/**
 * @brief   Returns total number of true grid points for all subdomains
 * @param	grid	Grid
 * @param   mpiInfo MpiInfo
 * @return	double	Number of total true nodes in simulation
 */

long int gTotTruesize(const Grid *grid, const MpiInfo *mpiInfo);

/**
 * @brief   Returns global location of a node in a 3D grid
 * @param	grid		Grid
 * @param   mpiInfo 	MpiInfo
 * @param   p			long int
 * @return	location	i,j,k values of node location
 */
void gNodeToGrid3D(Grid *grid, const MpiInfo *mpiInfo, long int p, double *pos);
/**
 * @brief Performs a central space finite difference on a grid
 * @param 	scalar 	Value to do the finite differencing on
 * @return	field	Field returned after derivating
 */

void gFinDiff1st(const Grid *scalar, Grid *field);

/**
 * @brief Performs a 2nd order central space finite difference on a grid
 * @param 	rho 	Value to do the finite differencing on
 * @return	phi		Field returned after derivating
 *
 *
 */
void gFinDiff2nd3D(Grid *phi,const Grid *rho);

/**
 * @brief Performs a 2nd order central space finite difference on a grid
 * @param 	rho 	Value to do the finite differencing on
 * @return	phi		Field returned after derivating
 *
 *
 */
void gFinDiff2ndND(Grid *phi,const Grid *rho);

 /**
 * @brief Normalize E-field
 * @param	ini		Input file dictionary
 * @param	E		E-field
 *
 * Normalizes an non-normalized (but non-dimensional) E-field store in a Grid
 * object according to step-size, time step and mass and charge of specie 0.
 */
void gNormalizeE(const dictionary *ini, Grid *E);

void gNormalizePhi(const dictionary *ini, Grid *phi);

/**
 * @brief Set total rho to 0
 * @param rho	rho-field
 *
 * Sets total charge density to 0. This can be useful to avoid infinite
 * potential when using periodic boundary conditions.
 *
 */

void gNeutralizeGrid(Grid *rho, const MpiInfo *mpiInfo);

/**
* @brief Adds a grid to another.
* @param	result		Grid added to
* @param	addition	Grid that is added to the other
*
* Adds one grid to another. result = result + addition
*
*/
void gAddTo(Grid *result, Grid *addition);

/**
* @brief Subtracts a grid from another.
* @param	result		Grid subtracted from
* @param	subtraction	Grid that is subtracted from the other
*
* Subtracts one grid to another. result = result - addition
*
*/

void gSubFrom(Grid *result, const Grid *subtraction);

/**
* @brief Sums up the values in the true grid
* @param	grid       Grid
* @param	mpiInfo	   MpiInfo
*
* Sums the true grid, uses a NDimensional algorithm
*
*/
double gSumTruegrid(const Grid *grid);


void gAssertNeutralGrid(const Grid *rho, const MpiInfo *mpiInfo);

/**
 * @brief Applies boundary conditions to edge
 * @param 	grid		Grid to apply boundary conditions to
 * @param	mpiInfo		Info about subdomain
 *
 * @return 	grid		Returns grid with changed boundary
 *
 * Applies boundary conditions
 */
void gBnd(Grid *grid, const MpiInfo *mpiInfo);

/**
 * @brief	Assign particles artificial positions suitable for debugging
 * @param			ini				Input file dictionary
 * @param[in,out]	grid			Grid
 *
 * A quantity will be artificially assigned values depending on the position of
 * the grid points. For scalar valued quantities, node j will have value j,
 * node (j,k) will have value j+k*10 and so on for higher dimensions. For
 * instance at node (j,k,l)=(4,5,6) the value will be 456.
 *
 * For vector valued quantities the integer part of all values at a given point
 * are the same (unless there is 10 values or more), but the decimal part
 * increments by 0.1 for each value.
 *
 * In addition, 1000*mpiRank is added to all values. For instance the third
 * value (z-component) at grid point (j,k,l)=(4,5,6) of the subdomain with MPI
 * rank 2 will be 2654.3.
 */
void gValDebug(Grid *grid, const MpiInfo *mpiInfo);

/**
 * @brief	Creates .grid.h5-file to store population in
 * @param	ini				Dictionary to input file
 * @param	grid			Grid
 * @param	mpiInfo			MpiInfo
 * @param	denorm			Quantity denormalization factor
 * @param	fName			Filename
 * @return	void
 * @see gWriteH5(), gCloseH5()
 *
 * An output file is created whose filename is as explained in openH5File().
 * Remember to call gCloseH5().
 *
 * The file will have one dataset in the root group for each time-step a grid
 * quantity is stored, named "n=<timestep>" where <timestep> is signified with
 * one decimal allowing for interleaved quantities.
 *
 * The grid data is stored in normalized values. To obtain physical values, the
 * output contains two attributes: The "Quantity denormalization factor" is a
 * factor that denormalizes the values of the stored quantity upon
 * multiplication. The values are stored on a integer-based grid. If these
 * integer position of the grid nodes are multiplied by the "Axis
 * denormalization factor", physical coordinates are obtained. "Axis
 * denormalization facotr" is basically the same as the step size in physical
 * units.
 */
void gOpenH5(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo,
			 const Units *units, double denorm, const char *fName);

/**
 * @brief	Store values in Grid to .grid.h5-file
 * @param	grid			Grid
 * @param	mpiInfo			MpiInfo
 * @param	n				Timestep of quantity to be stored
 * @return	void
 *
 * n is double to allow storing quantities at half time-steps
 * (e.g. leapfrog). All but the most significant decimals are discarded.
 *
 * The function will fail ungracefully if trying to write to an existing
 * dataset (testing omitted for performance reasons).
 */
void gWriteH5(const Grid *grid, const MpiInfo *mpiInfo, double n);

/**
 * @brief	Read values from .grid.h5-file to Grid
 * @param	grid			Grid
 * @param	mpiInfo			MpiInfo
 * @param	n				Timestep to read from file
 * @return	void
 *
 * n is double to allow reading quantities from half time-steps
 * (e.g. leapfrog). All but the most significant decimals are discarded.
 *
 * The function will fail ungracefully if trying to read from a non-existing
 * dataset (testing omitted for performance reasons).
 */
void gReadH5(Grid *grid, const MpiInfo *mpiInfo, double n);

/**
 * @brief	Closes .grid.h5-file
 * @param	grid		Grid
 * @return	void
 */
void gCloseH5(Grid *grid);

/**
 * @brief Creates a neighborhood in MpiInfo
 * @param			ini		Dictionary to input file
 * @param[in,out]	mpiInfo	MpiInfo
 * @param			grid	Some grid
 *
 * Prior to creating a neighborhood with this function domain decomposition
 * functions for particles (particle migration) will not work.
 */
void gCreateNeighborhood(const dictionary *ini, MpiInfo *mpiInfo, Grid *grid);

/**
 * @brief Destroys a neighborhood
 * @param	mpiInfo	MpiInfo
 */
void gDestroyNeighborhood(MpiInfo *mpiInfo);

/**
 * @brief Computes potential energy
 * @param		rho		Charge density
 * @param		phi		Electric potential
 * @param[out]	pop		Population to store results in
 *
 * Computes the potential energy as
 * \f[
 *	U=\sum_j\rho_j\phi_j
 * \f]
 *
 * Hence the computed energy will be for the same timestep as rho and phi is in.
 * Only total energy (for all species) is calculated, and hence, the result is
 * stored in potEnergy[nSpecies] in Population. Only the energy for this
 * subdomain is computed. The energy for the whole domain is gathered during
 * storing to h5-file.
 *
 * NB: rho and phi must be same kind of Grid.
 * NB: Assumes electrostatic approximation.
 */
void gPotEnergy(const Grid *rho, const Grid *phi, Population *pop);

/**
 * @brief Returns normalized size of global domain
 * @param	ini		Input file dictionary
 * @return	Array of normalized dimensions.
 *
 * trueSize*nSubdomains-1 for all dimension.
 */
int *gGetGlobalSize(const dictionary *ini);
long int gGetGlobalVolume(const dictionary *ini);

void dumpWholeGrid(dictionary *ini, Grid *grid);
void dumpTrueGrid(dictionary *ini, Grid *grid);
void fillGridIndexes(Grid *grid);

/**
 * @brief Removes halo (ghost nodes) from grid quantity
 * @param[in,out]	grid	Grid
 * @return			void
 *
 * Removes halo from grid by in-place operations. Grid can be accessed using
 * grid->size and grid->sizeProd as usual after removal.
 *
 * Note: The memory for grid->val is not re-allocated to consume less memory.
 * This is for optimization purposes. If a halo want to be temporarily removed,
 * and later re-inserted, this can be done without reallocation operations.
 */
void gRemoveHalo(Grid *grid);

/**
 * @brief Inserts a halo (ghost nodes) on grid quantity
 * @param[in,out]	grid			Grid
 * @param			nGhostLayers	Number of ghost layers to insert
 * @return			void
 *
 * Inserts a helo to grid  by in-place operations.
 *
 * To re-insert a halo previously removed by gRemoveHalo(), make a copy of
 * nGhostLayers (e.g. using memcpy()) from grid before removing halo, and use
 * it with this function when inserting halo.
 *
 * Note: The memory for grid->val is not re-allocated to get sufficient memory
 * by this function. The developer must make sure grid->val has sufficient
 * memory allocated prior to calling this function or segmentation fault will
 * occur. This is for optimization purposes.
 *
 * This function assumes grid has no ghost layers from before (as if
 * gRemoveHalo() has been called). Could be extended in future.
 */
void gInsertHalo(Grid *grid, const int *nGhostLayers);

/**
 * @brief Initial grid configurations, temporary untill a read hdf5 method is
 * available
 */

 #define PI 3.14159265

void gFillHeavi(Grid *grid, int dim,const MpiInfo *mpiInfo);
void gFillHeaviSol(Grid *grid, int rank, const MpiInfo *mpiInfo);
void gFillPoint(Grid *grid, const MpiInfo *mpiInfo);
void gFillPointSol(Grid *grid, const MpiInfo *mpiInfo);
void gFillSin(Grid *grid, int d, const MpiInfo *mpiInfo, int norm);
void gFillSinSol(Grid *grid, int d, const MpiInfo *mpiInfo);
void gFillSinESol(Grid *grid, int d ,const MpiInfo *mpiInfo);
void gFillExp(Grid *grid, const MpiInfo *mpiInfo);
void gFillRng(Grid *grid, const MpiInfo *mpiInfo, const gsl_rng *rng);
void gFillCst(Grid *grid, const MpiInfo *mpiInfo);
void gFillPolynomial(Grid *grid , const MpiInfo *mpiInfo);


#endif // GRID_H
