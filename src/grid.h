/**
 * @file		grid.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Input/output functions.
 * @date		26.06.16
 */

#ifndef GRID_H
#define GRID_H



/**
 * @brief Allocates a Grid object as specified in the input file
 * @param	ini			Input file
 * @param	nValues		Number of values per grid point
 * @return				Pointer to Grid
 *
 * Use nValues=1 for scalar field, nValues=3 for 3D vector field and so on.
 *
 * Remember to free using gFree().
 *
 * NB! Assumes 1 ghost point on all edges for now.
 */

Grid *gAlloc(const dictionary *ini, int nValues);

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
 *  Set the boundary slices. With Dirichlet and Neumann bnd conditions the values
 *  which corresponds to phi and the gradient of phi is set here.
 *
 *  For now it only has the capability of constant boundaries, but it could
 *  later be expanded without too much trouble.
 *
 *
 */

void gSetBndSlices(Grid *grid,MpiInfo *mpiInfo);

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
 * @param sliceOp			SliceOpPointer
 * @param *Grid	Grid struct
 * @param *mpiInfo		MpiInfo struct
 * @param d				Along which dimension it should exhange ghost cells
 *
 *
 * Since each true subdomain is surrounded by ghost layers, both to implement boundary conditions
 * and to facility communication between different processes each taking care of a spesified subdomain,
 * this function fills the ghost layers with the values from the surrounding subdomains. When called
 * it has the option to either set the ghost layer as neighboring subdomain, or add to the existing
 * ghost layers. This is specified by the first argument which can either be setSlice or addSlice.
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
 *	If we then want to the ghost layers in the subdomains with the values from the surrounding subdomain
 *	we use the following syntax.
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
 *	If we go back to the original grids and we instead want to add to the ghostlayer
 *	we can use addSlice instead.
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
 *	If needed it should be quick to facilitate for more slice operations, in addition to set and add.
 *
 * NB! Only works with 1 ghost layer.
 * @see gHaloOp
 */

void gHaloOpDim(SliceOpPointer sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, int reverse);

void gHaloOpDim2(SliceOpPointer sliceOp, Grid *grid, const MpiInfo *mpiInfo, int d, int inverse);

/**
 * @brief Send and recieves the overlapping layers of the subdomains
 * @param sliceOp			SliceOpPointer
 * @param *grid				Grid struct
 * @param *mpiInfo			MpiInfo struct
 *
 * A wrapper to the gHaloOpDim function, that is used when the user wants the interaction in
 * all the dimensions.
 *
 * NB! Only works with 1 ghost layer.
 * @see gExchangeSlice
 * @see gHaloOpDim
 * @see SliceOpPointer
 */
void gHaloOp(SliceOpPointer sliceOp, Grid *grid, const MpiInfo *mpiInfo, int reverse);

/**
 * @brief Extracts a (dim-1) dimensional slice of grid values.
 * @param	slice 		Return array
 * @param	grid		Grid
 * @param	d			Perpendicular direction to slice
 * @param	offset 		Offset of slice
 * @return				Void
 *
 * This function gets extracts a slice from a N dimensional grid. The integer d
 * decides in which direction the slice is perpendicular to, and the offset decides
 * which which slice it picks out. It needs a preallocated slice array where
 * the extracted slice will be stored.
 *
 * 2D example: Here we have a 5x4 grid and we want to extract a slice corresponding
 * to the second row, where x (d=0) is a constant 1.
 *
 * @code
 * 15   16   17   18   19
 *
 * 10   11   12   13   14
 *
 *  5    6    7    8    9

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
void gFinDiff2nd2D(Grid *phi,const Grid *rho);

 /**
 * @brief Normalize E-field
 * @param	ini		Input file dictionary
 * @param	E		E-field
 *
 * Normalizes an non-normalized (but non-dimensional) E-field store in a Grid
 * object according to step-size, time step and mass and charge of specie 0.
 */
void gNormalizeE(const dictionary *ini, Grid *E);

/**
 * @brief Set total rho to 0
 * @param rho	rho-field
 *
 *	Sets total charge density to 0. This can be useful to avoid infinite potential
 *	when using periodic boundary conditions.
 *
 */

void gNeutralizeGrid(Grid *rho, const MpiInfo *mpiInfo);

/**
* @brief Adds a grid to another.
* @param	result		Grid added to
* @param	addition	Grid that is added to the other
*
*	Adds one grid to another. result = result + addition
*
*/
void gAddTo(Grid *result, Grid *addition);

/**
* @brief Subtracts a grid from another.
* @param	result		Grid subtracted from
* @param	subtraction	Grid that is subtracted from the other
*
*	Subtracts one grid to another. result = result - addition
*
*/

void gSubFrom(Grid *result, Grid *subtraction);


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
 * @param	denorm			Quantity denormalization factors
 * @param	dimen			Quantity dimensionalizing factors
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
 * In PINC it is made an distinction between _non-dimensionalizing_ and
 * _normalizing_. Input quantities are non-dimensionalized by specifying them
 * in terms of Debye lengths, plasma frequency, elementary charges and so on
 * rather than using SI-units. Further on, the program normalizes them with
 * respect to for instance cell size in order to make computations as fast as
 * possible. The data stored in .grid.h5 is non-dimensionalized _and_ normalized
 * as it is often much cheaper to just rescale the axis in the visualization
 * tool rather than re-scaling all quantities in PINC.
 *
 * The file will have four attributes of size nDims attached to the root group
 * ("/") which is useful for interpreting the data. These are:
 *	- Axis denormalization factor
 *	- Axis dimensionalizing factor
 *	- Quantity denormalization factor
 *	- Quantity dimensionalizing factor
 *
 * The axis denormalization factor can be multiplied to the integer axis to
 * convert it to be in terms of Debye lengths. Another multiplication by axis
 * dimensionalizing factor converts the axes to meters. Likewise for the
 * quantity factors. However, since the quantity factors depend on which
 * quantity it is (e.g. charge density or electric field), it must be specified
 * in the inputs denorm and dimen in this function. They are expected to be of
 * length nDims.
 */
void gOpenH5(const dictionary *ini, Grid *grid, const MpiInfo *mpiInfo, const double *denorm, const double *dimen, const char *fName);

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
 * @brief	Read values from .grid.h5-fiel to Grid
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

/**
 * @brief Initial grid configurations, temporary untill a read hdf5 method is available
 */

 #define PI 3.14159265

void fillHeaviside(Grid *grid, const MpiInfo *mpiInfo);
void fillHeaviSol(Grid *grid, const MpiInfo *mpiInfo);
void fillPointCharge(Grid *grid, const MpiInfo *mpiInfo);
void fillPointSol(Grid *grid, const MpiInfo *mpiInfo);
void fillSin(Grid *grid, const MpiInfo *mpiInfo);
void fillSinSol(Grid *grid, const MpiInfo *mpiInfo);
void fillExp(Grid *grid, const MpiInfo *mpiInfo);
void fillRng(Grid *grid, const MpiInfo *mpiInfo, const gsl_rng *rng);
void fillCst(Grid *grid, const MpiInfo *mpiInfo);
void fillPolynomial(Grid *grid , const MpiInfo *mpiInfo);


#endif // GRID_H
