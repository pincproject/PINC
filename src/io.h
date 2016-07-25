/**
 * @file		io.h
 * @author		Sigvald Marholm <sigvaldm@fys.uio.no>,
 				Gullik Vetvik Killie <gullikvk@fys.uio.no>
 * @copyright	University of Oslo, Norway
 * @brief		Input/output functions.
 * @date		26.06.16
 */

#ifndef IO_H
#define IO_H

/**
 * @brief	The PINC equivalent of printf().
 * @param	kind	STATUS, WARNING or ERROR depending on what to output.
 * @param	format	printf-like format specifier
 * @param	...		printf-like arguments
 * @return	void
 * @see msgKind, fMsg(), printf()
 *
 * This replaces printf() in the PINC context. Similar syntax to printf().
 * In the case of an ERROR, the program is terminated. Appends end-of-line
 * automatically at the end.
 *
 * The message will by default be printed by all nodes calling msg(), however
 * kind can be bitwise ORed with ONCE to only allow the master to display this
 * message, e.g. STATUS|ONCE.
 */
void msg(msgKind kind, const char* restrict format,...);

/**
 * @brief	Prints message to file given by msgfiles:<fNameKey> in ini-file.
 * @param	ini			Input file dictionary
 * @param	fNameKey	Name of key holding the filename
 * @param	format		printf-like format specifier
 * @param	...			printf-like arguments
 * @return	void
 * @see		msg(), printf(), fprintf()
 *
 * This replaces fprintf() in the PINC context being the preferred way to write
 * messages (not data) to files. This guarantees that no file names are hard
 * coded but configurable from the ini-file.
 *
 * Beware that this is not a high performance writing function, and very
 * frequent invocations (for instance per particle) should be avoided.
 */
void fMsg(dictionary *ini, const char* restrict fNameKey, const char* restrict format, ...);

/**
 * @name ini functions
 * Opens .ini input files as dictionary and reads from it. The base
 * functionality is provided by the iniparser library which provides the
 * "dictionary" datatype, but its features is largely extended by the
 * ini-functions (e.g. error handling and array inputs) which serves as a layer
 * on top of iniparser. Do not call iniparser directly.
 *
 * The notation "section:key" is used to indicate keys under various sections
 * in the ini-file.
 *
 * Comma-separated lists (using "," as delimeter) is interpreted as arrays.
 *
 * Remember to free returned arrays.
 *
 * Example of use:
 *
 * @code
 *	dictionary *ini = iniOpen(argc,argv);
 *	int value = iniGetInt(ini,"section:key");
 *	iniClose(ini);
 * @endcode
 *
 * @param		ini			ini-file dictionary
 * @param		key			Key to get value from ("section:key")
 * @param[out]	nElements	Number of elements in returned array
 */
///@{

/**
 * @brief	Opens PINC input ini-file specified in PINC's arguments.
 * @brief	argc	Argument count (as passed to PINC)
 * @param	argv	Arugment vector (as passed to PINC)
 * @return	Allocated iniparser dictionary holding data from ini-file
 *
 * Performs sanity check on argc and argv and opens the input file specified in
 * argv[1]. It also empties all files specified in msgfiles.
 *
 * Following arguments are used to override settings from the input file, which
 * may be useful if running parameter sweeps from external scripts. E.g. if
 * argv[2]=="grid:nSubdomains=2,2,2" the value of this key is overriden in
 * the returned dictionary. Any number of such arguments can be accepted, and
 * in any order, but the key must already exist in the input file (a
 * restriction inherited from the underlying iniparser library).
 *
 * The special argument "getnp" is used to get the number of MPI processes PINC
 * requires within an external script, for instance to know which "-np"
 * argument to pass to mpirun. When PINC stumbles upon "getnp" it immediately
 * terminates after returning the requested number, even before processing
 * further arguments. Since changing settings (i.e. "grid:nSubdomains") may
 * change the required number of processes, "getnp" should always be the last
 * argument.
 *
 * Close dictionary using iniClose() after use.
 */
dictionary* iniOpen(int argc, char *argv[]);

///@brief Close dictionary
void iniClose(dictionary *ini);

///@brief Get integer
int iniGetInt(const dictionary* ini, const char *key);
///@brief Get long int
long int iniGetLongInt(const dictionary* ini, const char *key);
///@brief Get double
double iniGetDouble(const dictionary* ini, const char *key);
///@brief Allocate and get string (remeber to free)
char* iniGetStr(const dictionary *ini, const char *key);

///@brief Allocate and get array of integers (remember to free)
int* iniGetIntArr(const dictionary *ini, const char *key, int *nElements);
///@brief Allocate and get array of long ints (remember to free)
long int* iniGetLongIntArr(const dictionary *ini, const char *key, int *nElements);
///@brief Allocate and get array of doubles (remember to free)
double* iniGetDoubleArr(const dictionary *ini, const char *key, int *nElements);

/**
 * @brief Get the array of strings associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements in returned array
 * @return			NULL-terminated array of NULL-terminated strings
 *
 * Output is similar to listToStrArr(). Remember to free resulting string array
 * using freeStrArr().
 */
char** iniGetStrArr(const dictionary *ini, const char *key, int *nElements);

///@brief Get the number of elements in an array/comma-separated list
int iniGetNElements(const dictionary* ini, const char* key);

/**
 * @brief Assert that a number of entries are arrays of equal length.
 * @param			ini			Dictionary to search
 * @param			nKey		Number of keys to search for
 * @param			...			Keys to search for
 * @return			Number of elements in arrays
 * @see				iniGetNElements()
 *
 * Example:
 * @code
 *	iniAssertEqualNElements(ini,3,"mySec:a","mySec:b","mySec:c");
 * @endcode
 *
 * This code does nothing if the specified entries have equal length. Unequal
 * lengths will lead to an error.
 */
int iniAssertEqualNElements(const dictionary *ini, int nKeys, ...);

///@}

/**
 * @brief Frees dynamically allocated NULL-terminated array of NULL-terminated strings
 * @param	strArr	Pointer to array of strings
 * @return	void
 * @see listToStrArr(), iniGetStrArr()
 */
void freeStrArr(char** strArr);

/**
 * @brief Open (or create) an .h5-file
 * @param	ini		Input file dictionary
 * @param	fName	File name
 * @param	fSubExt	File sub-extension
 * @return	HDF5 file identifier
 * @see		H5Fcreate(), H5Fclose()
 *
 * This function opens an existing .h5 file or creates it if it doesn't
 * already exist. It takes the place of H5Fopen() and H5Fcreate() in PINC. The
 * file name is <fName>.<fSubExt>.h5 or possibly <prefix>_<fName>.<fSubExt>.h5
 * (see below).
 *
 * Similarly as .h5 indicates the file type being h5, <fSubExt> indicates _what
 * kind_ of .h5-file it is, e.g. if it is a grid quantity (.grid.h5) or
 * population data (.pop.h5). These kinds are standardized within PINC and have
 * dedicated functions unifying file handling of these data types. fName is a
 * name given by the developer to further specify the _contents_ of the file.
 * E.g. is it charge density (rho.grid.h5) or electric field (E.grid.h5).
 *
 * The file will be stored in the folder specified by "files:output" in the
 * input file. Examples of valid values of "files:output":
 *	- output = data/
 *  - output = ./data/
 *	- output = ../data/
 *	- output = ~/data/
 *	- output = /home/me/data/
 *
 * However, this variable also allows the user to specify a prefix to each file.
 * This is indicated by _not_ ending the variable with '/'. For instance:
 *	- output = prefix
 *	- output = data/prefix
 *
 * This will output files such a "prefix_rho.grid.h5".
 *
 * Parent directories are created unless they already exists.
 *
 * Close return value using H5Fclose().
 */
hid_t openH5File(const dictionary* ini, const char *fName, const char *fSubExt);

/**
 * @brief Sets array of double as attributes in h5-file
 * @param	h5		.h5-file identifier
 * @param	name	Attribute name
 * @param	value	Attribute value
 * @param	size	Size of attribute
 *
 * If the attribute name already exists, the old attributes are overwritten.
 * This function is not capable of setting multi-dimensional or non-double
 * arrays as attributes (that's part of its simplification compared to the
 * functions in the HDF5 library).
 */
void setH5Attr(hid_t h5, const char *name, const double *value, int size);

/**
 * @brief Creates a group in a .h5-file recursively
 * @param	h5		.h5-file identifier
 * @param	name	Group name
 * @return	void
 *
 * createH5Group() handles creating multiple levels of groups recursively unlike
 * H5Gcreate() which will crash ungracefully. It also will not crash if a group
 * already exists. Note that the last group must have a trailing slash.
 *
 * Examples:
 *	- /group/group/dataset	- Will create /group/group
 *	- /group/group/		 	- Will create /group/group
 */
void createH5Group(hid_t h5, const char *name);

/**
 * @brief Creates a .xy.h5-file for storing (x,y) datasets
 * @param	ini		Dictionary to input file
 * @param	name	File name
 * @return	Handle for H5-file
 *
 * For conventions regarding the file name, see openH5File().
 * See xyWriteH5() for how to write (x,y) datapoits to the file.
 * Remember to close using xyCloseH5() or PINC will fail ungracefully.
 */
hid_t xyOpenH5(const dictionary *ini, const char *fName);

/**
 * @brief Closes a .xy.h5-file
 * @param	h5		Identifier to h5-file to close
 * @param			void
 */
void xyCloseH5(hid_t h5);

/**
 * @brief Writes an (x,y) datapoint to a dataset in an H5-file
 * @param	h5		.h5-file identifier
 * @param	name	Dataset name
 * @param	x		x-value
 * @param	y		y-value
 * @param	op		MPI reduction operation performed on y
 * @return	void
 *
 * A datapoint on a curve in the xy-plane is appended at the end of a dataset in
 * the h5-file. For instance, the x-axis can represent the time steps and the y
 * axis some kind of energy or error residual.
 *
 * In parallel executions, the y value is reduced across all MPI nodes using the
 * specified MPI reduction operation, for instance MPI_SUM to sum the y-value of
 * all MPI nodes before writing to file. If the x value differs amongst the
 * nodes, the x-value of rank 0 is simply used.
 *
 * The dataset must be created beforehand by calling xyCreateDataset() and the
 * file is created by xyOpenH5(). Remember to close the H5 file using
 * xyCloseH5().
 *
 * Example:
 * @code
 *	hid_t hist = xyOpenH5(ini,"timesweep");
 *
 *	xyCreateDataset(hist,"/energy/potential");
 *	xyCreateDataset(hist,"residual");
 *
 *	for(int n=0;n<N;n++){
 *		...
 *		// Energy summed across nodes using MPI_SUM and stored to file
 *		xyWrite(hist,"/energy/potential",(double)n,energy,MPI_SUM);
 *		xyWrite(hist,"residual",(double)n,res,MPI_SUM);
 *		...
 *	}
 *	xyCloseH5(hist);
 * @endcode
 */
void xyWrite(hid_t h5, const char* name, double x, double y, MPI_Op op);

/**
 * @brief Creates a dataset in a .xy.h5 file
 * @param	h5		Identifier to .h5-file to create dataset in
 * @param	name	Dataset name
 * @return			void
 *
 * Creates a dataset specified by its name and all parent groups with it. For
 * example, see xyWrite().
 */
void xyCreateDataset(hid_t h5, const char *name);

/**
 * @brief Writes grid structs to a parsefile
 * @param ini 		dictionary of the input file
 * @param grid 		grid struct
 *
 * Debug help
 */

#endif // IO_H
