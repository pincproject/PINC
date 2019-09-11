/**
 * @file		aux.h
 * @brief		Auxiliary functions.
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 *				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 */

#ifndef AUX_H
#define AUX_H

/**
 * @name Timer functions
 */
///@{

/**
 * @brief	Allocates a Timer struct
 * @return	Pointer to Timer struct
 * @see		Timer, tFree, tMsg()
 *
 * Remember to free using tFree().
 */
Timer *tAlloc();

/**
 * @brief	Frees a Timer struct allocated with tAlloc()
 * @param	timer 	Pointer to Timer struct
 * @see		Timer, tAlloc()
 */
void tFree(Timer *t);


/**
 * @brief	Starts the time counter
 * @param	timer 	Pointer to Timer struct
 * @see		tStop()
 *
 *	Sets t->start = clock().
 *
 */
void tStart(Timer *t);


/**
 * @brief	Stops the timer stopwatch
 * @param	timer 	Pointer to Timer struct
 * @see		Timer, tStart()
 *
 *	Computes the time, since tStart was called and then the result is adtimerded to the
 *	total time the stopwatch have been running.
 *
 */
void tStop(Timer *t);


/**
 * @brief	Resets the timer
 * @param	timer 	Pointer to Timer struct
 * @see		Timer, tStart()
 */
void tReset(Timer *t);

/**
 * @brief	Prints a message along with a time converted to a sensible format
 * @param	nanoSec 	nanoseconds
 * @see		string		message before the time measurement
 *
 * Prints the message along with the the time measurement, converted from
 * nanoseconds to a suitable format.
 *
 * Useful for testing execution speed of chunks of code. Each call to tMsg()
 * prints the total time the program has been running, along with the time since
 * last call to tMsg() before it resets the timer.
 *
 */
void tMsg(long long int nanoSec, const char *string);

///@}

/**
 * @brief Concatenates strings
 * @param	n	Number of strings to concatenate
 * @param	...	List of pointer to strings to concatenate
 * @return		Pointer to concatenated string
 * @see	strcat()
 *
 * Contrary to strcat() this function allocates a new string of suitable length.
 * Make sure to call free().
 */
char *strCatAlloc(int n,...);

/** @name Array functions
 * This group of functions perform simple and common array/vector operations on
 * one array/vector 'a' (unary operations) or two arrays/vectors 'a' and 'b'
 * (binary operations). Some operations may have scalar-valued results, in which
 * case the return value of the function is used, while others have
 * vector-valued results, in which case a variable 'res' must be pre-allocated
 * to hold the result. Example:
 *
 * @code
 *	int n = 5;
 *	int a[] = {1,2,3,4,5};
 *	int b[] = {2,3,4,5,6};
 *	int *res = malloc(n*sizeof(*res));
 *	int s = aiSum(a,n);	// The sum of all elements in a (the scalar 15)
 *	aiAdd(a,b,res,n);		// The sum of a and b (the vector {3,5,7,9,11})
 * @endcode
 *
 * All functions starts with the prefix 'a' (for array) followed by one or
 * two letters signifying the datatype of the arrays:
 *
 *  prefix  | datatype
 *	--------|------------------------------------
 *	ad      | double
 *	ai      | int
 *	al      | long int
 *	ail		| int but output is a long int array
 *
 * The return values of certain 'ai'-functions may still be long int since, for
 * instance, the product of many int's may end up in the long int range. If the
 * return value is stored in a plain int variable the return value will cast
 * nicely to int without emitting a warning. Likewise, the number of elements
 * of the arrays (n) is always of type long int to support large arrays, but
 * using a plain int input should cast nicely. The only place where long int
 * doesn't cast nicely to int is in arrays since the array would mis-align upon
 * casting. Therefore, some 'ai' functions returning int arrays as results has
 * an 'ail' variant which returns an array of long ints (for when the result
 * may be in the long int range).
 *
 * In-place operations are supported, i.e. the 'res' array may very well be the
 * same as one or both of the input array. Example:
 *
 * @code
 *	int n = 3;
 *	int a[] = {1,2,3}
 *	aiMul(a,a,a,n);		// a now equals {1,4,9}
 * @endcode
 *
 * Finally notice that pointer arithmetics along with adjusted values of 'n' may
 * be utilized to affect only parts of an array. Example:
 *
 * @code
 *	int n = 5;
 *	int a[] = {1,1,1,1,1};
 *	aiSetAll(a,n,2);		// a now equals {2,2,2,2,2}
 *	aiSetAll(a+1,3,5);		// a now equals {2,5,5,5,2}
 *	aiSetAll(&a[1],3,7);	// a now equals {2,7,7,7,2}
 * @endcode
 *
 * @param		a		Input array/vector (unary operations)
 * @param		b		Input array/vector (unary and binary operations)
 * @param[out]	res		Resulting array/vector (if result is vector-valued)
 * @param		n		Number of elements in a (and b and typically res)
 * @param		value	Input scalar value to use in operation (if any)
 * @return				Resulting scalar (if result is scalar-valued)
 */
///@{
///@brief Adds two arrays
void adAdd(const double *a, const double *b, double *res, long int n);
///@brief Adds two arrays
void aiAdd(const int *a, const int *b, int *res, long int n);
///@brief Adds two arrays
void alAdd(const long int *a, const long int *b, long int *res, long int n);
///@brief Multiplies two arrays element-wise (Hadamard)
void adMul(const double *a, const double *b, double *res, long int n);
///@brief Multiplies two arrays element-wise (Hadamard)
void aiMul(const int *a, const int *b, int *res, long int n);
///@brief Multiplies two arrays element-wise (Hadamard)
void alMul(const long int *a, const long int *b, long int *res, long int n);
///@brief Shifts an array by a constant value
void adShift(double *a, long int n, double value);
///@brief Shifts an array by a constant value
void aiShift(int *a, long int n, int value);
///@brief Shifts an array by a constant value
void alShift(long int *a, long int n, long int value);
///@brief Units an array by a constant value
void adScale(double *a, long int n, double value);
///@brief Units an array by a constant value
void aiScale(int *a, long int n, int value);
///@brief Units an array by a constant value
void alScale(long int *a, long int n, int value);
///@brief Returns maximum value in an array
double adMax(const double *a, long int n);
///@brief Returns unit normal vector. Assumes 3D eucledian vector
void adNormal(const double *a, const double *b, double *res, long int n);
///@brief Returns vector reflected by an arbitrary plane
void adReflect(const double *ray, const double *a, const double *b, double *res);
///@brief Returns maximum value in an array
int aiMax(const int *a, long int n);
///@brief Returns maximum value in an array
long int alMax(const long int *a, long int n);
///@brief Returns minimum value in an array
double adMin(const double *a, long int n);
///@brief Returns minimum value in an array
int aiMin(const int *a, long int n);
///@brief Returns minimum value in an array
long int alMin(const long int *a, long int n);
///@brief Returns most significant extremum (maximum or minimum) of an array.
/// For instance, if -6 is the minimum and 5 is the maximum the return value is
/// -6.
double adExt(const double *a, long int n);
///@brief Returns most significant extremum (maximum or minimum) of an array.
/// For instance, if -6 is the minimum and 5 is the maximum the return value is
/// -6.
int aiExt(const int *a, long int n);
///@brief Returns most significant extremum (maximum or minimum) of an array.
/// For instance, if -6 is the minimum and 5 is the maximum the return value is
/// -6.
long int alExt(const long int *a, long int n);
///@brief Returns sum of all elements
double adSum(const double *a, long int n);
///@brief Returns sum of all elements
long int aiSum(const int *a, long int n);
///@brief Returns sum of all elements
long int alSum(const long int *a, long int n);
///@brief Returns average of all elements
double adAvg(const double *a, long int n);
///@brief Returns average of all elements
double aiAvg(const int *a, long int n);
///@brief Returns average of all elements
double alAvg(const long int *a, long int n);
///@brief Returns product of all elements
double adProd(const double *a, long int n);
///@brief Returns product of all elements
long int aiProd(const int *a, long int n);
///@brief Returns product of all elements
long int alProd(const int *a, long int n);
///@brief Returns cross product of two eucledian 3D vectors
void adCrossProd(const double *a, const double *b, double *res); 

///@brief Returns dot product of vectors
int adDotProd(const double *a, const double *b, long int n);
///@brief Returns dot product of vectors
int aiDotProd(const int *a, const int *b, long int n);
///@brief Returns dot product of vectors
int alDotProd(const long int *a, const long int *b, long int n);
///@brief Returns 1 if arrays are equal, 0 otherwise. Every element must be
/// at maximum 'tol' apart to be considered equal (max-norm).
int adEq(const double *a, const double *b, long int n, double tol);
///@brief Returns 1 if arrays are equal, 0 otherwise
int aiEq(const int *a, const int *b, long int n);
///@brief Returns 1 if arrays are equal, 0 otherwise
int alEq(const long int *a, const long int *b, long int n);
///@brief Determine cumulative product of elements in 'a' starting at 1.
/// Hence the cumulative product of {5,4,3} is {1,5,20,60}. Notice that the
/// result is of lenght n+1 in this case.
void adCumProd(const double *a, double *res, long int n);
///@brief Determine cumulative product of elements in 'a' starting at 1.
/// Hence the cumulative product of {5,4,3} is {1,5,20,60}. Notice that the
/// result is of lenght n+1 in this case.
void aiCumProd(const int *a, int *res, long int n);
///@brief Determine cumulative product of elements in 'a' starting at 1.
/// Hence the cumulative product of {5,4,3} is {1,5,20,60}. Notice that the
/// result is of lenght n+1 in this case.
void ailCumProd(const int *a, long int *res, long int n);
///@brief Determine cumulative product of elements in 'a' starting at 1.
/// Hence the cumulative product of {5,4,3} is {1,5,20,60}. Notice that the
/// result is of lenght n+1 in this case.
void alCumProd(const long int *a, long int *res, long int n);
///@brief Determine cumulative sum of elements in 'a' starting at 0.
/// Hence the cumulative sum of {5,4,3} is {0,5,9,12}. Notice that the
/// result is of lenght n+1 in this case.
void adCumSum(const double *a, double *res, long int n);
///@brief Determine cumulative sum of elements in 'a' starting at 0.
/// Hence the cumulative sum of {5,4,3} is {0,5,9,12}. Notice that the
/// result is of lenght n+1 in this case.
void aiCumSum(const int *a, int *res, long int n);
///@brief Determine cumulative sum of elements in 'a' starting at 0.
/// Hence the cumulative sum of {5,4,3} is {0,5,9,12}. Notice that the
/// result is of lenght n+1 in this case.
void ailCumSum(const int *a, long int *res, long int n);
///@brief Determine cumulative sum of elements in 'a' starting at 0.
/// Hence the cumulative sum of {5,4,3} is {0,5,9,12}. Notice that the
/// result is of lenght n+1 in this case.
void alCumSum(const long int *a, long int *res, long int n);
///@brief Sets all elements in array to 'value'
void adSetAll(double *a, long int n, double value);
///@brief Sets all elements in array to 'value'
void aiSetAll(int *a, long int n, int value);
///@brief Sets all elements in array to 'value'
void alSetAll(long int *a, long int n, long int value);
///@brief Set n elements in array manually, e.g. adSet(a,5,1.,2.,3.,4.,5.);
void adSet(double *a, long int n, ...);
///@brief Set n elements in array manually, e.g. adSet(a,5,1.,2.,3.,4.,5.);
void aiSet(int *a, long int n, ...);
///@brief Set n elements in array manually, e.g. adSet(a,5,1.,2.,3.,4.,5.);
void alSet(long int *a, long int n, ...);
///@brief See adPrint(). varName is the name to output for the variable.
void adPrintInner(double *a, long int inc, long int end, char *varName);
///@brief See aiPrint(). varName is the name to output for the variable.
void aiPrintInner(int *a, long int inc, long int end, char *varName);
///@brief See alPrint(). varName is the name to output for the variable.
void alPrintInner(long int *a, long int inc, long int end, char *varName);
///@brief Prints an array in a nice format (for debugging only).
#define adPrint(a,n) do { adPrintInner(a,1,n,#a); } while (0)
///@brief Prints an array in a nice format (for debugging only).
#define aiPrint(a,n) do { aiPrintInner(a,1,n,#a); } while (0)
///@brief Prints an array in a nice format (for debugging only).
#define alPrint(a,n) do { alPrintInner(a,1,n,#a); } while (0)
///@brief Prints an array in a nice format (for debugging only).
#define adPrintSpaced(a,inc,end) do { adPrintInner(a,inc,end,#a); } while (0)
///@brief Prints an array in a nice format (for debugging only).
#define aiPrintSpaced(a,inc,end) do { aiPrintInner(a,inc,end,#a); } while (0)
///@brief Prints an array in a nice format (for debugging only).
#define alPrintSpaced(a,inc,end) do { alPrintInner(a,inc,end,#a); } while (0)
///@}

#endif // AUX_H
