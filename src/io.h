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
 * @brief	Function selector for assigning function pointers
 * @param	ini		Input file dictionary
 * @param	key		Key to method parameter ("section:key")
 * @param	...		Addresses to set-functions to select from
 * @return 	Function pointer as specified by set-function
 *
 * The select() macro enables the developer to let a method be selectable by
 * the user through the input files. This is best described through an example.
 * Consider that two different numerical methods may be used for the same thing,
 * each represented by its own function:
 *
 * @code
 *	void fun1(const char *str){
 *		msg(STATUS|ONCE,"Function 1 called with input: %s",str);
 *	}
 *
 *	void fun2(const char *str1, const char *str2){
 *		msg(STATUS|ONCE,"Function 2 called with inputs: %s, %s",str1,str2);
 *	}
 * @endcode
 *
 * Both of these functions _must_ return void, but they can have any kind of
 * input parameters. fun() may then be a function pointer, selected to point to
 * either of these two functions, and subsequently used as a normal function:
 *
 * @code
 *	void (*fun)() = select(ini,"method:fun",fun1_set,fun2_set);
 *	fun("Input 1","Input 2");
 * @endcode
 *
 * If "method:fun" in the ini-file equals "fun1", the output will read
 *
 * @code
 *	Function 1 called with input: Input 1
 * @endcode
 *
 * whereas if it reads "fun2", the output reads
 *
 * @code
 *	Function 2 called with inputs: Input 1, Input 2
 * @endcode
 *
 * Notice that the number of arguments in fun1() and fun2() may be different.
 * Unused arguments are simply ignored. However, overlapping arguments must be
 * the same kind since fun() is invoked the same way each time.
 *
 * Finally, notice that select() do not actually return the addresses to fun1()
 * and fun2() itself, but instead calls a _set()-function (internally) to do it.
 * These can be defined as follows:
 *
 * @code
 *	funPtr fun1_set(dictionary *ini){
 *		return fun1;
 *	}
 *
 *	funPtr fun2_set(dictionary *ini){
 *
 *		int nDims = iniGetInt(ini,"grid:nDims");
 *		if(nDims!=3) msg(ERROR|ONCE, "fun2() requires nDims=3");
 *
 *		return &fun2; // Ampersand optional
 *	}
 * @endcode
 *
 * Thus, when you write "fun1_set" in the select() call, fun1_set is the address
 * of the function being called when select finds a matching value "fun1" in the
 * input file. Underscore merely acts as a syntactical delimeter, and select()
 * matches only the first part with the input file.
 *
 * The purpose of having _set()-functions rather than just setting the function
 * pointers directly is that the developer of each method can (and should)
 * attach sanity checks to the methods inside the _set() functions. Thus every
 * time a new, selectable method is developed, the developer must also make a
 * _set() function specifying/checking the circumstances under which the method
 * works. In our example, it is assumed that the fun2() method only works for
 * 3-dimensional problems. It may not be beneficial to read the ini file and
 * test for 3-dimensionality inside fun2() since it may run inside a time
 * critical loop. The _set() functions, however, are only called once. The
 * _set() functions must always be made with one input, the ini-file, and with
 * return type funPtr (which is just a pointer to a function returning void).
 *
 * In some cases, the numerical methods pointed to by e.g. fun() may need to be
 * initialized, and in different ways depending on which method is being used.
 * Then, since underscore acts merely as a delimeter, one may create
 * _init()-functions which trigger on the same key in the ini-file. E.g.
 *
 * @code
 *	void (*fun)()     = select(ini,"method:fun",fun1_set ,fun2_set );
 *	void (*initFun)() = select(ini,"method:fun",fun1_init,fun2_init);
 *	initFun();
 *	fun();
 * @endcode
 */
#define select(ini,key,...) selectInner(ini,key,#__VA_ARGS__,__VA_ARGS__)

/**
 * @brief	Function selector for assigning function pointers
 * @param	ini		Input file dictionary
 * @param	key		Key to method parameter ("section:key")
 * @param	string	Variadic arguments stringified by select() macro
 * @param	...		Addresses to set-functions to select from
 * @return 	Function pointer as specified by set-function
 *
 * Merely a wrapper which stringifies the function pointers in the variadic
 * arguments to make the selectInner() function able to match function names by
 * values in the ini-file.
 */
funPtr selectInner(dictionary *ini, const char *key, const char *string,...);

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
 * Example of use:
 *
 * @code
 *	dictionary *ini = iniOpen(argc,argv);
 *	int value = iniGetInt(ini,"section:key");
 *	iniClose(ini);
 * @endcode
 *
 *
 * In the case of fetching arrays the returned array will have as many elements
 * as specified by nElements. If the entry in the ini-file has fewer elements
 * than specified they will be repeated, e.g. if 5 elements is read from the
 * entry "1,2" the elements 1,2,1,2,1 will be returned.
 *
 * Remember to free returned arrays.
 *
 * @param		ini			ini-file dictionary
 * @param		key			Key to get value from ("section:key")
 * @param		nElements	Number of elements to fetch (in case of arrays)
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
int* iniGetIntArr(const dictionary *ini, const char *key, int nElements);
///@brief Allocate and get array of long ints (remember to free)
long int* iniGetLongIntArr(const dictionary *ini, const char *key, int nElements);
///@brief Allocate and get array of doubles (remember to free)
double* iniGetDoubleArr(const dictionary *ini, const char *key, int nElements);

/**
 * @brief Get the array of strings associated to a key.
 * @param			ini			Dictionary to search
 * @param			key			Key string to look for
 * @param[out]		nElements	Number of elements to get
 * @return			NULL-terminated array of NULL-terminated strings
 *
 * Output is similar to listToStrArr(). Remember to free resulting string array
 * using freeStrArr().
 */
char** iniGetStrArr(const dictionary *ini, const char *key, int nElements);

///@brief Get the number of elements in an array/comma-separated list
int iniGetNElements(const dictionary* ini, const char* key);

/**
 * @brief Apply multiplicator to entries in ini-file with suffix.
 * @param[in,out]	ini			Dictionary to search
 * @param			key			Key string to look for
 * @param			suffix		Suffix string to look for
 * @param			mul			Multiplier(s)
 * @param			mulLen		Length of mul
 * @return			void
 *
 * Changes the dictionary for later use. Searches through elements in a list
 * specified by key and, and if suffix is present, multiplies each entry by
 * a multiplier. Entry i in the list is multiplied by mul[i%mulLen]. Notice
 * that this is purely a "parsing" feature. I.e. it does not know how to expand
 * a single element into three, that is done on a later stage by get-functions.
 */
void iniApplySuffix(dictionary *ini, const char *key, const char *suffix, const double *mul, int mulLen);

void parseIndirectInput(dictionary *ini);

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
